#ifndef _actual_eos_H_
#define _actual_eos_H_

#include <extern_parameters.H>
#include <fundamental_constants.H>
#include <prob_parameters.H>
#include <network.H>
#include <iostream>
#include <vector>
#include <actual_eos_data.H>
#include <eos_util.H>

const std::string eos_name = "csm";

inline
void actual_eos_init ()
{

    // constant ratio of specific heats
    if (eos_gamma > 0.e0_rt) {
        gamma_const = eos_gamma;
    } else {
        gamma_const = 5.0_rt / 3.0_rt;
    }

}

template <typename I>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool is_input_valid (I input)
{
  static_assert(std::is_same<I, eos_input_t>::value, "input must be an eos_input_t");

  bool valid = true;

  if (input == eos_input_rh ||
      input == eos_input_tp ||
      input == eos_input_ps ||
      input == eos_input_ph ||
      input == eos_input_th) {
    valid = false;
  }

  return valid;
}

template <typename I, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_eos (I input, T& state)
{
    static_assert(std::is_same<I, eos_input_t>::value, "input must be an eos_input_t");

    const Real R = C::k_B * C::n_A;

    Real p, e, cv, poverrho;

	Real ndens[n_elems]; // number density of each element

	state.mu = 0.0_rt; 


	// if we're not doing simple grey opacity and we are within the temperature bounds <T_cold,T_saha>
	// and we aren't in the ambient medium (aux[4]>0) then do the full saha solve to get n_e (aka rhoYe) and mu
	if(!problem::do_grey_es && state.T < problem::T_saha && state.T > problem::T_cold &&  state.aux[4] < 0.9_rt){



    brent_args args;


//	args.T = std::max(state.T,2.e3_rt); //don't want to run into issues at low temps

	Real z_a, z_b, z_avg; // bounds for Zbar <z_a,z_b>
	// we want to make sure the solution z_avg (set below) is bracketed by these values
	// technically fulll neutrality Z_avg->0 while full ionized goes to Z_avg->Z

	Real mu_n = 0.0_rt, mu_i = 0.0_rt, ne_i = 0.0_rt;
	Real n_e = 0.0_rt;
	std::vector<std::vector<Real>> ion_fracs;

	Real n_ion = 0.0_rt;
	
	for(int i = 0; i < n_elems; ++i){ //solve for some estimate values, not really used, here for debugging
		mu_n += state.xn[i]*A[i]; // mean molecular weight if fully neutral
		mu_i += state.xn[i]/A[i]*(1.0_rt+Z[i]); //mean molecular weight if fully ionized
		ne_i += Z[i]*state.xn[i]/A[i]; // electron density if fully ionized

	}

	mu_i = 1.0/mu_i;
	ne_i *= state.rho/C::m_p;

	
	for(int i = 0; i < n_elems; ++i){

		z_avg = 1.0;//this doesn't matter, it gets overridden anyways
		z_b = 1.1*Z[i]; // z_b is always <=Z[i], but make it a little larger in case of rounding error
		z_a = min(1.e-6,pow(10.0_rt,-(50.0_rt-pow(log10(state.T),2.0)))); // kind of a kluge 
		// basically what i'm trying to do is adjust the lower bound on Zbar to be temperature dependent
		// since as you go to higher temperatures, Zbar increases, but as T gets really low, Zbar->0
		// i played around with some different compositions and this seemed to work OK for CSM conditions

	
		std::vector<Real> ifrac;

		for(int j = 0; j < num_I[i]; ++j) ifrac.push_back(0.0);

		ion_fracs.push_back(ifrac);

		ndens[i] = state.xn[i]*state.rho/(A[i]*C::m_p);
		state.mu += ndens[i]*A[i];
		n_ion += ndens[i];
		args.k = i;
		args.n = ndens[i];
		args.ion_frac = ion_fracs[i];
		args.f = new Real[num_I[i]-1];

		int solve_error = find_root(z_a,z_b,&args,z_avg);

		if(solve_error == -1) std::cerr  << "Z = " << Z[i] << " root not bracketed!\n";
		else if (solve_error == 0) std::cerr << "Z = " << Z[i] << " max iter reached!\n";

		ion_fracs[i][0] = args.ion_frac[0];

		for(int j = 1; j < num_I[i]; ++j){

			n_e += j*args.ion_frac[j]*ndens[i];
			ion_fracs[i][j] = args.ion_frac[j];

		}


	}

	state.mu /= n_ion + n_e;

	state.aux[0] = C::m_p*n_e/state.rho; //this should really be Ye, but since
										// we aren't doing neutrinos
										// Ye only appears in the opacity
										// and we want the terms to work out so that kappa*rhoYe = n_e*sigma_T
	state.aux[1] = 1.0_rt/state.mu; // aux[1] stores the inverse of mu



	} // if we aren't doing the saha solve
	else if(problem::do_grey_es || state.T >= problem::T_saha || state.aux[4] > 0.9_rt){

		//if using grey opacity or T>T_saha or ambient medium we skip the saha solve
		// and just assume completely ionized gas
		//
		Real ne_i = 0.0, invmu_i = 0.0;

		for(int i = 0; i < n_elems; ++i){

			invmu_i += state.xn[i]/A[i]*(1.0_rt+Z[i]);
			ne_i += Z[i]*state.xn[i]/A[i];

		}

		ne_i *= state.rho/C::m_p;

		if(state.aux[4] > 0.9_rt) ne_i *= 1.e-2_rt; // just to make ambient medium opacity smaller
		
		state.mu = 1.0_rt/invmu_i;
		state.aux[0] = C::m_p*ne_i/state.rho;
		state.aux[1] = invmu_i;



	}
	else{ // completely neutral, this should be hit if T<T_cold

		Real mu_n = 0.0;

		for(int i =0; i < n_elems; ++i) mu_n += state.xn[i]*A[i];

		state.mu = mu_n;
		state.aux[1] = 1.0_rt/mu_n;
		state.aux[0] = C::m_p/state.rho*1.e-10_rt; //some small nonzero number
		

	}


    switch (input) {

    case eos_input_rt:

        // dens, temp and xmass are inputs
        cv = R / (state.mu * (gamma_const - 1.0_rt));
        e = cv * state.T;
        if constexpr (has_energy<T>::value) {
            state.cv = cv;
            state.e = e;
        }
        if constexpr (has_pressure<T>::value) {
            state.p = (gamma_const - 1.0_rt) * state.rho * e;
        }
        if constexpr (has_pressure<T>::value && has_energy<T>::value) {
            state.gam1 = gamma_const;
        }

        break;

    case eos_input_rh:

        // dens, enthalpy, and xmass are inputs
#ifndef AMREX_USE_GPU
        amrex::Error("EOS: eos_input_rh is not supported in this EOS.");
#endif

        break;

    case eos_input_tp:

        // temp, pres, and xmass are inputs
#ifndef AMREX_USE_GPU
        amrex::Error("EOS: eos_input_tp is not supported in this EOS.");
#endif

        break;

    case eos_input_rp:

        // dens, pres, and xmass are inputs

        if constexpr (has_pressure<T>::value) {
            poverrho = state.p / state.rho;
            state.T = poverrho * state.mu * (1.0_rt / R);
            if constexpr (has_energy<T>::value) {
                state.e = poverrho * (1.0_rt / (gamma_const - 1.0_rt));
                state.gam1 = gamma_const;
            }
        }

        break;

    case eos_input_re:

        // dens, energy, and xmass are inputs

        if constexpr (has_energy<T>::value) {
            poverrho = (gamma_const - 1.0_rt) * state.e;
            state.T = poverrho * state.mu * (1.0_rt / R);

            if constexpr (has_pressure<T>::value) {
                state.p = poverrho * state.rho;
                state.gam1 = gamma_const;

                // sound speed
                state.cs = sqrt(gamma_const * poverrho);

                state.dpdr_e = poverrho;
                state.dpde = (gamma_const - 1.0_rt) * state.rho;
            }
        }

        // Try to avoid the expensive log function.  Since we don't need entropy
        // in hydro solver, set it to an invalid but "nice" value for the plotfile.
        if constexpr (has_entropy<T>::value) {
            state.s = 1.0_rt;
        }

        break;

    case eos_input_ps:

        // pressure entropy, and xmass are inputs

#ifndef AMREX_USE_GPU
        amrex::Error("EOS: eos_input_ps is not supported in this EOS.");
#endif

        break;

    case eos_input_ph:

        // pressure, enthalpy and xmass are inputs
#ifndef AMREX_USE_GPU
        amrex::Error("EOS: eos_input_ph is not supported in this EOS.");
#endif

        break;

    case eos_input_th:

        // temperature, enthalpy and xmass are inputs

        // This system is underconstrained.
#ifndef AMREX_USE_GPU
        amrex::Error("EOS: eos_input_th is not a valid input for the gamma law EOS.");
#endif

        break;

    default:

#ifndef AMREX_USE_GPU
        amrex::Error("EOS: invalid input.");
#endif

        break;

    }

}



inline
void actual_eos_finalize ()
{
}

#endif

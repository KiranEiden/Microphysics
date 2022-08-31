#ifndef actual_opacity_H
#define actual_opacity_H

#include <AMReX.H>

#include <extern_parameters.H>
#include <prob_parameters.H>
#include <fundamental_constants.H>

AMREX_INLINE
void actual_opacity_init () {}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_opacity (Real& kp, Real& kr, Real rho, Real temp, Real rhoYe, Real nu,
                     bool get_Planck_mean, bool get_Rosseland_mean)
{




		Real kappa_es = problem::kap;

		if(get_Planck_mean) kp = kappa_es*rho*problem::eps;
		if(get_Rosseland_mean) kr = kappa_es*rho;


}

#endif

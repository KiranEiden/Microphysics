#ifndef VODE_RHS_SIMPLIFIED_SDC_H
#define VODE_RHS_SIMPLIFIED_SDC_H

#include <network.H>
#include <burn_type.H>
#include <extern_parameters.H>

#include <vode_type_simplified_sdc.H>
#include <actual_network.H>
#include <actual_rhs.H>
#include <numerical_jacobian.H>
#ifdef NONAKA_PLOT
#include <nonaka_plot.H>
#endif

// The f_rhs routine provides the right-hand-side for the DVODE solver.
// This is a generic interface that calls the specific RHS routine in the
// network you're actually using.

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void rhs(const Real time, burn_t& state, dvode_t& vode_state, RArray1D& ydot, const bool in_jacobian=false)
{

    // update rho

    fill_unevolved_variables(time, state, vode_state);

    // ensure that the mass fractions are valid -- only vode_state is
    // updated here

    clean_state(time, state, vode_state);

    // convert to the burn_t -- this does an EOS call to get T
    // and populates the (burn_t) state

    vode_to_burn(time, vode_state, state);

    // make sure that the temperature is valid

    if (state.T <= EOSData::mintemp || state.T >= MAX_TEMP) {

        for (int n = 1; n <= VODE_NEQS; ++n) {
            ydot(n) = 0.0_rt;
        }

        return;

    }

    // call the specific network to get the RHS

    actual_rhs(state, ydot);

#ifdef NONAKA_PLOT
    if (! in_jacobian) {
        nonaka_rhs(time, state, ydot);
    }
#endif

    // apply fudge factor:

    if (react_boost > 0.0_rt) {
        // the correct size here is the size of the reactive state RHS,
        // not the SDC state RHS (neqs, not VODE_NEQS)
        for (int n = 1; n <= neqs; ++n) {
            ydot(n) *= react_boost;
        }
    }

    // convert back to the vode type -- this will add the advective terms

    rhs_to_vode(time, state, vode_state, ydot);

}



// Analytical Jacobian

template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void jac (burn_t& state, dvode_t& vode_state, MatrixType& pd)
{
    // NOTE: the time at which to evaluate the Jacobian is not
    // explicitly passed. VODE always evaluates the analytic
    // Jacobian at vode_state.tn, and we pass vode_state.


    // update rho, rho*u, ... in the burn_t state

    fill_unevolved_variables(vode_state.tn, state, vode_state);

    // ensure that the mass fractions are valid -- only vode_state is
    // updated here

    clean_state(vode_state.tn, state, vode_state);

    // convert to the burn_t -- this does an EOS call to get T
    // and populates the (burn_t) state

    vode_to_burn(vode_state.tn, vode_state, state);


    if (state.T <= EOSData::mintemp || state.T >= MAX_TEMP) {

        for (int j = 1; j <= VODE_NEQS; ++j) {
            for (int i = 1; i <= VODE_NEQS; ++i) {
                pd(i,j) = 0.0_rt;
            }
        }

        return;

    }

    if (vode_state.jacobian_type == 3) {
        jac_info_t jac_info;
        jac_info.h = vode_state.H;

        // the numerical Jacobian for SDC will automatically
        // put it into the correct form (in terms of energy)
        // so we can just operate on the VODE Jacobian storage
        // directly

        numerical_jac(state, jac_info, pd);

        // apply fudge factor:

        if (react_boost > 0.0_rt) {
            pd.mul(react_boost);
        }

        // jacobian = 3 is a hybrid numerical + analytic
        // implementation for simplified SDC.  It uses a one-sided
        // difference numerical approximation for the reactive
        // part and then algebraically transforms it to the
        // correct state variables for SDC.  This means we need to
        // account for the NumSpec+1 RHS evals
        vode_state.NFE += NumSpec + 1;

    } else {

        // Call the specific network routine to get the Jacobian.

        actual_jac(state, pd);

        // The Jacobian from the nets is in terms of dYdot/dY, but we want
        // it was dXdot/dX, so convert here.
        for (int n = 1; n <= NumSpec; n++) {
            for (int m = 1; m <= neqs; m++) {
                pd(n,m) = pd(n,m) * aion[n-1];
            }
        }

        for (int m = 1; m <= neqs; m++) {
            for (int n = 1; n <= NumSpec; n++) {
                pd(m,n) = pd(m,n) * aion_inv[n-1];
            }
        }

        // apply fudge factor:

        if (react_boost > 0.0_rt) {
            pd.mul(react_boost);
        }

        // The system we integrate has the form (rho X_k, rho e)

        // pd is now of the form:
        //
        //  SFS         / d(rho X1dot)/dX1  d(rho X1dit)/dX2 ... 1/cv d(rho X1dot)/dT \
        //              | d(rho X2dot)/dX1  d(rho X2dot)/dX2 ... 1/cv d(rho X2dot)/dT |
        //  SFS-1+nspec |   ...                                                       |
        //  SEINT       \ d(rho Edot)/dX1   d(rho Edot)/dX2  ... 1/cv d(rho Edot)/dT  /
        //
        //                   SFS                                         SEINT

    }

}
#endif

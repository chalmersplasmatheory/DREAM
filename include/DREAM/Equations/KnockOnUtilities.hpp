/**
 * Implementation of calculations that are used to describe knock-on collisions.
 *
 * Handles quantities that show up in both the general as well as the approximate
 * (Chiu-Harvey, Rosenbluth-Putvinski) knock-on collision terms.
 *
 * The basis for the physics/maths implemented herein is the Møller scattering of
 * relativistic electron-electron collisions, specialized to the case where the
 * target particle is stationary.
 */
#ifndef _DREAM_EQUATIONS_KNOCK_ON_UTILITIES_HPP
#define DREAM_EQUATIONS_KNOCK_ON_UTILITIES_HPP

#include "DREAM/Constants.hpp"
#include "DREAM/DREAMException.hpp"
#include "FVM/Grid/Grid.hpp"
#include <cmath>

namespace DREAM::KnockOnUtilities {

// Helper functions for the delta integrand
real_t evaluateF(real_t y, real_t xi1, real_t xi_star);
real_t evaluateD(real_t t1, real_t t2, real_t xi1, real_t xi_star);
void evaluateT1T2(
    real_t &t1, real_t &t2, real_t xi0_j1, real_t xi0_j2, real_t BOverBmin, real_t xi0Cutoff
);
void evaluateZ1Z2(real_t &z1, real_t &z2, real_t t1, real_t t2, real_t xi_star);

// the integrand that is integrated over theta inside EvaluateDeltaContribution.
real_t deltaIntegrand(
    real_t theta, len_t ir, real_t xi0_j1, real_t xi0_j2, real_t xi01, real_t xi_star,
    const FVM::Grid *grid
);

enum orbitIntegrationQuadrature {
    QUADRATURE_MIDPOINT,
    QUADRATURE_ADAPTIVE
};


/**
 * Gyro-, bounce and FVM averaged kinematic delta function, defined as
 * the function \delta_j in the docs/notes/theory notes.
 * Important: multiplicity effects in the trapped regions are not accounted for,
 * e.g. it needs to be summed for co- and counter moving trapped particles
 * for both xi0 and xi01 when used in practice.
 */

real_t EvaluateDeltaContribution(
    len_t ir, real_t xi_star, real_t xi01, real_t xi0_f1, real_t xi0_f2, real_t Vp1, real_t theta1,
    real_t theta2, len_t n_points_integral, const FVM::Grid *grid, orbitIntegrationQuadrature quad=QUADRATURE_MIDPOINT
);

real_t EvaluateDelta(
    len_t ir, real_t xi_star, real_t xi01, real_t xi0_f1, real_t xi0_f2, real_t Vp1, real_t theta1,
    real_t theta2, len_t n_points_integral, const FVM::Grid *grid
);

inline void estimateBoundingTheta(
    len_t ir, len_t j, len_t l, real_t &theta1, real_t &theta2, const FVM::Grid *grid
) {
    real_t xi0_f2 = grid->GetMomentumGrid(ir)->GetP2_f(j + 1);
    if (xi0_f2 < 0) {
        theta1 = grid->GetThetaBounce1_f2(ir, 0, j);
        theta2 = grid->GetThetaBounce2_f2(ir, 0, j);
    } else {
        theta1 = grid->GetThetaBounce1_f2(ir, 0, j + 1);
        theta2 = grid->GetThetaBounce2_f2(ir, 0, j + 1);
    }
    theta1 = std::max(theta1, grid->GetThetaBounce1(ir, 0, l));
    theta2 = std::min(theta2, grid->GetThetaBounce2(ir, 0, l));
}

/**
 * Cosine of the angle between incident electron and knock-on,
 * as function of their respective momenta.
 *   p1: momentum magnitude of the incident electron (assumed to
 *       collide with stationary electron).
 *   p: knock-on momentum magnitude.
 */
inline real_t evaluateXiStar(real_t p, real_t p1) {
    // handle the well-behaved (and common) limit of infinitely energetic incident electron
    real_t g = sqrt(1 + p * p);
    if (std::isinf(p1)) {
        return sqrt((g - 1) / (g + 1));
    }
    real_t g1 = sqrt(1 + p1 * p1);
    return sqrt((g1 + 1) * (g - 1) / ((g1 - 1) * (g + 1)));
}

/**
 * Moller differential cross section, between collisions (p1, 0) -> (p, p')
 * where p1 is the momentum of the incident particle (on a stationary one),
 * and p is interpreted as the knock-on momentum (when gamma1 > 2*gamma - 1)
 */
inline real_t evaluateMollerDifferentialCrossSection(real_t p, real_t p1) {
    real_t prefactor = 2 * M_PI * Constants::r0 * Constants::r0 * Constants::c;
    real_t g = sqrt(1 + p * p);
    real_t g1 = sqrt(1 + p1 * p1);
    real_t term1 = (g1 - 1) * (g1 - 1);
    real_t term2 =
        -(g - 1) * (g1 - g) / (g1 * g1) * (2 * g1 * g1 + 2 * g1 - 1 - (g - 1) * (g1 - g));
    return prefactor * (term1 + term2) * g1 * g1 /
           (p1 * p1 * (g - 1) * (g - 1) * (g1 - g) * (g1 - g));
}

/**
 * Returns the function 'S' obtained by writing the
 * Mœller differential cross section Sigma on the form
 *   Sigma = gamma/p * dS/dp.
 * Is obtained when exactly integrating the knock-on source
 * over a FVM momentum grid cell (in the knock-on momenum p).
 */

inline real_t evaluateMollerFlux(real_t p, real_t p1) {
    real_t prefactor = 2 * M_PI * Constants::r0 * Constants::r0 * Constants::c;

    real_t g = sqrt(1 + p * p);
    if (std::isinf(p1)) {
        return -prefactor / (g - 1);
    }

    real_t g1 = sqrt(1 + p1 * p1);
    real_t term1 = g1 * g1 * (1 / (g1 - g) - 1 / (g - 1));
    real_t term2 = g;
    real_t term3 = -(2 * g1 - 1) / (g1 - 1) * log((g - 1) / (g1 - g));

    return prefactor / (p1 * p1) * (term1 + term2 + term3);
}

} // namespace DREAM::KnockOnUtilities
#endif /*_DREAM_EQUATIONS_KNOCK_ON_UTILITIES_HPP*/

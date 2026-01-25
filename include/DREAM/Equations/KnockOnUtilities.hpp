/**
 * Implementation of calculations that are used to describe knock-on collisions.
 *
 * This header provides:
 *  - Relativistic kinematics and Møller cross-section utilities (namespace Kinematics)
 *  - Analytic utilities for the gyro- and pitch-averaged kinematic delta function
 *    (namespace AnalyticDelta)
 *  - Orbit-averaged delta-function integrals
 *  - Grid-aware helpers for assembling delta matrix elements on DREAM grids
 *
 * Handles quantities that show up in both the general as well as the approximate
 * (Chiu-Harvey, Rosenbluth-Putvinski) knock-on collision terms.
 *
 * The basis for the physics/maths implemented herein is the Møller scattering of
 * relativistic electron-electron collisions, specialized to the case where the
 * target particle is stationary.
 */
#ifndef _DREAM_EQUATIONS_KNOCK_ON_UTILITIES_HPP
#define _DREAM_EQUATIONS_KNOCK_ON_UTILITIES_HPP

#include <cmath>

#include "DREAM/Constants.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/RadialGrid.hpp"

namespace DREAM::KnockOnUtilities {

// ============================================================================
// Low-level analytic building blocks (no grid dependence)
// ============================================================================
namespace Kinematics {
// Functions related to relativistic collision kinematics and Møller cross sections.

/**
 * Cosine of the angle between incident electron and knock-on,
 * as function of their respective momenta.
 *   p1: momentum magnitude of the incident electron (assumed to
 *       collide with stationary electron).
 *   p: knock-on momentum magnitude.
 */
inline real_t EvaluateXiStar(real_t p, real_t p1) {
    // handle the well-behaved (and common) limit of infinitely energetic incident electron
    real_t g = sqrt(1 + p * p);
    if (std::isinf(p1)) {
        return sqrt((g - 1) / (g + 1));
    }
    real_t g1 = sqrt(1 + p1 * p1);
    return sqrt((g1 + 1) * (g - 1) / ((g1 - 1) * (g + 1)));
}

/**
 * Møller differential cross section, between collisions (p1, 0) -> (p, p')
 * where p1 is the momentum of the incident particle (on a stationary one),
 * and p is interpreted as the knock-on momentum (when gamma1 > 2*gamma - 1)
 */
inline real_t EvaluateMollerDifferentialCrossSection(real_t p, real_t p1) {
    real_t prefactor = 2 * M_PI * Constants::r0 * Constants::r0 * Constants::c;
    real_t g = sqrt(1 + p * p);
    if (std::isinf(p1)) {
        return prefactor / ((g-1)*(g-1));
    }
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
 *   Sigma = gamma/p * dS/dp = dS/dgamma.
 * Is obtained when exactly integrating the knock-on source
 * over a FVM momentum grid cell (in the knock-on momentum p).
 */

inline real_t EvaluateMollerFlux(real_t p, real_t p1) {
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
}  // namespace Kinematics

namespace AnalyticDelta {
/**
 * Analytic, theta-local utilities for evaluating the kinematic delta function
 * after gyro- and pitch-averaging, with no geometry or grid dependence.
 */
real_t EvaluateXiAntiderivative(real_t y, real_t xi1, real_t xi_star);
real_t EvaluateXiIntegral(real_t t1, real_t t2, real_t xi1, real_t xi_star);
void ComputeXiIntegrationBounds(
    real_t &t1, real_t &t2, real_t xi0_f1, real_t xi0_f2, real_t BOverBmin, real_t xi0Cutoff
);
void ComputeXi1Bounds(real_t &z1, real_t &z2, real_t t1, real_t t2, real_t xi_star);
real_t EvaluateLocalContribution(
    real_t xi0_f1, real_t xi0_f2, real_t xi01, real_t BOverBmin, real_t xi_star, real_t xi0Cutoff,
    real_t &xi1_over_xi01
);
}  // namespace AnalyticDelta

}  // namespace DREAM::KnockOnUtilities
#endif /*_DREAM_EQUATIONS_KNOCK_ON_UTILITIES_HPP*/

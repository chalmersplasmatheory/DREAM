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

#include "DREAM/Equations/KnockOnQuadrature.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/RadialGrid.hpp"

namespace DREAM::KnockOnUtilities {

// ============================================================================
// Low-level analytic building blocks (no grid dependence)
// ============================================================================
namespace Kinematics {
real_t EvaluateXiStar(real_t p, real_t p1);
real_t EvaluateMollerDifferentialCrossSection(real_t p, real_t p1);
real_t EvaluateMollerFlux(real_t p, real_t p1);
real_t MaximumKnockOnMomentum(real_t p1);
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

// ============================================================================
// Orbit-averaging routines based on a RadialGrid for local-geometry evaluation
// ============================================================================

/**
 * Evaluate local geometric quantities and return true if this poloidal location
 * theta can be reached both by some xi0 in [xi0_f1, xi0_f2] and by xi01.
 */
bool CheckIfReachableAndSetGeometricQuantities(
    len_t ir, real_t theta, real_t xi0_f1, real_t xi0_f2, real_t xi01, real_t &BOverBmin,
    real_t &Jacobian, real_t &xi0Cutoff, FVM::FluxSurfaceAverager *fsa
);

// Evaluate the theta-local integrand of the orbit-averaged \delta_j function.
real_t OrbitDeltaIntegrand(
    real_t theta, len_t ir, real_t xi0_f1, real_t xi0_f2, real_t xi01, real_t xi_star,
    const FVM::RadialGrid *rg
);

// Unnormalized primary-volume-weighted orbit numerator. Its theta integral
// equals (Vp1/p1^2)*delta_j before trapped-sign/cell corrections. The existing
// integrand already divides by the outgoing cell width; do not divide again.
// This explicit-pitch interface introduces no grid-centre substitution.
real_t OrbitNumeratorIntegrand(
    real_t theta, len_t ir, real_t xi0_f1, real_t xi0_f2, real_t xi01, real_t xi_star,
    const FVM::RadialGrid *rg
);

// Integrate the raw numerator using prepared geometry/rule. When mirror=true,
// both scattering signs share each evaluation of the geometric quantities.
real_t EvaluateOrbitNumerator(
    len_t ir, real_t xi_star, real_t xi01, real_t xi0_f1, real_t xi0_f2,
    const FVM::RadialGrid *rg, Quadrature::OrbitGeometry& geometry,
    const Quadrature::GaussLegendreRule& rule, bool mirror=false
);

// Continuum integral of the unnormalized numerator over represented primary
// pitch. Includes trapped signs, but no native void-column or boundary-row
// transfer. The secondary interval is averaged, the primary interval integrated.
real_t EvaluatePitchCellNumerator(
    len_t ir, real_t xi_star, Quadrature::Interval primary, Quadrature::Interval secondary,
    const FVM::RadialGrid *rg, Quadrature::OrbitGeometry& geometry,
    const Quadrature::GaussLegendreRule& poloidalRule,
    const Quadrature::GaussLegendreRule& pitchRule
);

/**
 * Fundamental gyro-, bounce-, and FVM-averaged kinematic delta function.
 *
 * Corresponds to the quantity \delta_j in the theory notes.
 * Note: trapping multiplicity effects are NOT included; callers must
 * explicitly account for co- and counter-moving trapped contributions.
 */
real_t EvaluateOrbitAveragedDelta(
    len_t ir, real_t xi_star, real_t xi01, real_t xi0_f1, real_t xi0_f2, real_t Vp1,
    const FVM::RadialGrid *rg, len_t poloidalOrder
);

/**
 * Orbit-averaged delta function \delta_j including DREAM trapping conventions.
 *
 * This wraps EvaluateOrbitAveragedDelta() and adds mirrored contributions
 * for trapped regions in xi0 and/or xi01.
 */
real_t EvaluateOrbitAveragedDeltaWithTrappingCorrection(
    len_t ir, real_t xi_star, real_t xi01, real_t xi0_f1, real_t xi0_f2, real_t Vp1,
    const FVM::RadialGrid *rg, len_t poloidalOrder
);

// ============================================================================
// Grid-aware helpers for constructing matrix elements in equation terms
// ============================================================================

/**
 * Unnormalized cell-averaged matrix element A_jl/Hgrid_l on a DREAM grid.
 * This is a diagnostic; complete columns are normalized by the helper below.
 */
real_t EvaluateDeltaMatrixElementOnGrid(
    len_t ir, real_t xi_star, len_t j, len_t l, const FVM::Grid *grid_knockon,
    const FVM::Grid *grid_primary, len_t poloidalOrder, len_t pitchOrder=4
);

/**
 * Set the angular column for primary pitch cell l. Normalize the column
 * so its outgoing pitch integral is one, and return the raw integral
 * Hquad_l / Hgrid_l. A void primary cell gives a zero column and integral.
 */
real_t SetDeltaMatrixColumnOnGrid(
    len_t ir, real_t xi_star, len_t l, const FVM::Grid *grid_knockon, const FVM::Grid *grid_primary,
    real_t *deltaCol, len_t poloidalOrder, len_t pitchOrder=4
);

// Reuse the prepared geometry and rules across columns on this surface.
real_t SetDeltaMatrixColumnOnGrid(
    len_t ir, real_t xi_star, len_t l, const FVM::Grid *grid_knockon, const FVM::Grid *grid_primary,
    real_t *deltaCol, Quadrature::OrbitGeometry &geometry,
    const Quadrature::GaussLegendreRule &poloidalRule, const Quadrature::GaussLegendreRule &pitchRule
);

// Set a single momentum matrix element S_{ik} on a DREAM grid.
real_t EvaluateMollerFluxMatrixElementOnGrid(
    len_t i, len_t k, const FVM::Grid *grid_knockon, const FVM::Grid *grid_primary, real_t pCutoff
);

// Evaluate the total Møller flux on the entire knock-on grid,
// identical to sum_i dp_i * EvaluateMollerFluxMatrixElementOnGrid(i,k,...).
real_t EvaluateMollerFluxIntegratedOverKnockonGrid(
    len_t k, const FVM::Grid *grid_knockon, const FVM::Grid *grid_primary, real_t pCutoff
);

}  // namespace DREAM::KnockOnUtilities
#endif /*_DREAM_EQUATIONS_KNOCK_ON_UTILITIES_HPP*/

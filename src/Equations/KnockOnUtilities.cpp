/**
 * Implementation of calculations that show up in kinetic models of
 * knock-on collisions of relativistic electrons (e.g. Møller scattering).
 *
 * The theory is described in detail in doc/notes/theory.tex,
 * and variable names throughout this file mirrors the corresponding
 * ones in the notes.
 *
 * There are a few recurring definitions throughout this implementation.
 * xi01 - is the incident electron pitch
 * [xi0_f1, xi0_f2] - is the knock-on pitch cell range over which we average the source.
 * xi_star - is the cosine of the true scattering angle (which is smeared out
 *   due to gyro-, orbit- and FVM averaging).
 */
#include "DREAM/Equations/KnockOnUtilities.hpp"

#include <gsl/gsl_integration.h>

#include <cmath>
#include <limits>

#include "DREAM/DREAMException.hpp"

using namespace DREAM;

namespace {

constexpr real_t RELAXED_EPS = 5 * std::numeric_limits<real_t>::epsilon();

constexpr real_t QUAD_ADAPTIVE_ATOL = 1e-6;

/*
 * RELAXED_EPS:
 * Used to guard against accumulated floating-point roundoff in chained
 * geometric and momentum-space transformations (e.g. xi0 → xi → xi1),
 * particularly near trapped–passing and endpoint boundaries.
 *
 * Chosen to be safely above machine epsilon while remaining negligible
 * on physical scales.
 *
 * QUAD_ADAPTIVE_ATOL:
 * Tolerance for terminating the adaptive integration algorithm.
 * The value was chosen empirically to balance computational speed and accuracy,
 * and is shown in unit testing to be sufficient to meet tolerances.
 */

}  // anonymous namespace

/**
 * Analytic antiderivative of the gyro-averaged kinematic delta function.
 * Evaluates the variable "F" defined in the theory notes.
 */
real_t KnockOnUtilities::AnalyticDelta::EvaluateXiAntiderivative(
    real_t y, real_t xi1, real_t xi_star
) {
    return asin((y - xi1 * xi_star) / sqrt((1 - xi1 * xi1) * (1 - xi_star * xi_star))) / M_PI;
}

/**
 * Analytic pitch integration of the gyro-averaged kinematic delta function,
 * from xi=t1 to xi=t2 (clipped to the kinematically available region [a-b, a+b],
 * in which case the values of F are limited to +/- 0.5).
 * Evaluates the variable "D" defined in the theory notes.
 */
real_t KnockOnUtilities::AnalyticDelta::EvaluateXiIntegral(
    real_t t1, real_t t2, real_t xi1, real_t xi_star
) {
    real_t a = xi_star * xi1;
    real_t b = sqrt(std::max(0.0, (1 - xi_star * xi_star) * (1 - xi1 * xi1)));
    real_t F1 = -0.5;
    real_t F2 = 0.5;
    if (t2 < a + b) {
        F2 = EvaluateXiAntiderivative(t2, xi1, xi_star);
    }
    if (t1 > a - b) {
        F1 = EvaluateXiAntiderivative(t1, xi1, xi_star);
    }

    return F2 - F1;
}

/**
 * Compute the local bounds [t1, t2] obtained when changing from low-field side pitch [xi0_f1,
 * xi0_f2] in the FVM cell average to their local values xi(xi0, theta). The derivation of t1 and t2
 * is given in the theory notes.
 *
 * xi0Cutoff is the minimum low-field side pitch that can reach
 * this poloidal location (i.e. this value of BOverBmin).
 */
void KnockOnUtilities::AnalyticDelta::ComputeXiIntegrationBounds(
    real_t &t1, real_t &t2, real_t xi0_f1, real_t xi0_f2, real_t BOverBmin, real_t xi0Cutoff
) {
    // Clamp local pitch to zero where xi0 cannot reach this theta
    if (-RELAXED_EPS <= xi0_f1 && xi0_f1 <= xi0Cutoff) {
        t1 = 0;
    } else {
        t1 = FVM::MomentumGrid::XiFromXi0(xi0_f1, BOverBmin);
    }
    if (-xi0Cutoff <= xi0_f2 && xi0_f2 <= RELAXED_EPS) {
        t2 = 0;
    } else {
        t2 = FVM::MomentumGrid::XiFromXi0(xi0_f2, BOverBmin);
    }
}

/**
 * Evaluates the local limits [z1, z2] that xi1 spans when xi is swept between [t1, t2]
 * (e.g. when xi0 is cell-averaged over [xi0_f1, xi0_f2]), at fixed scattering angle cosine xi_star.
 * These limits are derived in the theory notes.
 */
void KnockOnUtilities::AnalyticDelta::ComputeXi1Bounds(
    real_t &z1, real_t &z2, real_t t1, real_t t2, real_t xi_star
) {
    if (t2 < -xi_star) {
        z1 = t2 * xi_star - sqrt(1 - t2 * t2) * sqrt(1 - xi_star * xi_star);
    } else if (t1 <= -xi_star && -xi_star <= t2) {
        z1 = -1;
    } else {
        z1 = t1 * xi_star - sqrt(1 - t1 * t1) * sqrt(1 - xi_star * xi_star);
    }

    if (t2 < xi_star) {
        z2 = t2 * xi_star + sqrt(1 - t2 * t2) * sqrt(1 - xi_star * xi_star);
    } else if (t1 <= xi_star && xi_star <= t2) {
        z2 = 1;
    } else {
        z2 = t1 * xi_star + sqrt(1 - t1 * t1) * sqrt(1 - xi_star * xi_star);
    }
}

/**
 * Evaluates the analytically gyro- and xi-cell-averaged kinematic
 * delta function at a fixed poloidal location, characterized
 * by its BOverBmin and xi0Cutoff values.
 *
 * In the theory notes, the returned value corresponds to the D_j function.
 *
 * Returns 0 if no contribution exists at this theta, and
 * handles the xi01 = +/-1 (CH/RP) limits carefully.
 */

real_t KnockOnUtilities::AnalyticDelta::EvaluateLocalContribution(
    real_t xi0_f1, real_t xi0_f2, real_t xi01, real_t BOverBmin, real_t xi_star, real_t xi0Cutoff,
    real_t &xi1_over_xi01
) {
    // perform the checks on whether the theta is viable
    // the "true" condition is z1 <= xi1 <= z2 for us to
    // take this theta into account. However, roundoff
    // errors when xi1 == z1, z2 can cause invalid values
    // to slip through and crash other parts of the calculation.
    // For the important edge case xi01 == +/- 1,
    // which define the Chiu-Harvey and Rosenbluth-Putvinski
    // approximations, we seek these bounding values, and therefore
    // employ an exacter constraint by inserting values explicitly.
    real_t t1;
    real_t t2;
    AnalyticDelta::ComputeXiIntegrationBounds(t1, t2, xi0_f1, xi0_f2, BOverBmin, xi0Cutoff);

    // handle the xi01 == 1 case (used in the approximate collision terms--RP and CH)
    if (fabs(fabs(xi01) - 1.0) < RELAXED_EPS) {
        xi1_over_xi01 = 1.0;
        // exact solution of z2==1 (since z = xi1 = 1)
        if (xi01 > 0.0 && !(t1 < xi_star && xi_star < t2)) {
            return 0.0;
        }
        // exact solution of z1==-1 (since z = xi1 = -1)
        if (xi01 < 0.0 && !(t1 < -xi_star && -xi_star < t2)) {
            return 0.0;
        }
        return 1.0 / (xi0_f2 - xi0_f1);
    }
    real_t z1;
    real_t z2;
    AnalyticDelta::ComputeXi1Bounds(z1, z2, t1, t2, xi_star);
    real_t xi1 = FVM::MomentumGrid::XiFromXi0(xi01, BOverBmin);
    if (fabs(xi01) > RELAXED_EPS) {
        xi1_over_xi01 = xi1 / xi01;
    } else {
        xi1_over_xi01 = 1.0;
    }
    if (!(z1 + RELAXED_EPS < xi1 && xi1 < z2 - RELAXED_EPS)) {
        return 0.0;
    }
    return AnalyticDelta::EvaluateXiIntegral(t1, t2, xi1, xi_star) / (xi0_f2 - xi0_f1);
}

/**
 * Evaluates the geometric quantities BOverBmin, Jacobian and xi0Cutoff at this poloidal location.
 * Returns true if this poloidal location can be reached both by some xi0 in [xi0_f1, xi0_f2] and by
 * xi01.
 */
bool KnockOnUtilities::CheckIfReachableAndSetGeometricQuantities(
    len_t ir, real_t theta, real_t xi0_f1, real_t xi0_f2, real_t xi01, real_t &BOverBmin,
    real_t &Jacobian, real_t &xi0Cutoff, FVM::FluxSurfaceAverager *fsa
) {
    real_t Bmin = fsa->GetBmin(ir, FVM::FLUXGRIDTYPE_DISTRIBUTION);
    real_t B_;
    real_t ROverR0_;
    real_t NablaR2_;
    fsa->GeometricQuantitiesAtTheta(
        ir, theta, B_, Jacobian, ROverR0_, NablaR2_, FVM::FLUXGRIDTYPE_DISTRIBUTION
    );
    BOverBmin = B_ / Bmin;
    xi0Cutoff = sqrt(1 - 1 / BOverBmin);

    if (xi0Cutoff >= fabs(xi01)) {
        // xi01 orbit cannot reach this poloidal location => return 0.
        return false;
    }
    if (xi0Cutoff >= std::max(fabs(xi0_f1), fabs(xi0_f2))) {
        // no xi0 orbit in [xi0_f1, xi0_f2] can reach this poloidal location => return 0.
        return false;
    }
    return true;
}

/**
 * Orbit-dependent integrand of the analytically gyro- and FVM-averaged
 * kinematic delta function.
 *
 * Returns zero whenever kinematic or trapping constraints are violated
 * at the given poloidal angle theta.
 */
real_t KnockOnUtilities::OrbitDeltaIntegrand(
    real_t theta, len_t ir, real_t xi0_f1, real_t xi0_f2, real_t xi01, real_t xi_star,
    const FVM::RadialGrid *rg
) {
    real_t BOverBmin;
    real_t Jacobian;
    real_t xi0Cutoff;
    bool isReachable = CheckIfReachableAndSetGeometricQuantities(
        ir, theta, xi0_f1, xi0_f2, xi01, BOverBmin, Jacobian, xi0Cutoff,
        rg->GetFluxSurfaceAverager()
    );
    if (!isReachable) {
        return 0.0;
    }

    real_t xi1_over_xi01;
    real_t D = AnalyticDelta::EvaluateLocalContribution(
        xi0_f1, xi0_f2, xi01, BOverBmin, xi_star, xi0Cutoff, xi1_over_xi01
    );

    real_t sqrtG = Jacobian * FVM::MomentumGrid::evaluatePXiMetricOverP2(xi1_over_xi01, BOverBmin);

    return sqrtG * D;
}

namespace {
// Purely numerical helpers for quadrature, containing no knock-on physics.

// check if value if effectively zero
bool nearly_zero(real_t y, real_t tol = RELAXED_EPS) { return fabs(y) < tol; }

real_t trapezoid(real_t x1, real_t x2, real_t f1, real_t f2) { return (x2 - x1) * (f1 + f2) * 0.5; }

// Adaptively integrate a given interval with bisection (if needed) + trapz
template <typename Integrand>
real_t bisect_and_integrate(
    real_t x1, real_t x2, real_t f1, real_t f2, Integrand f, bool at_lower_endpoint,
    bool at_upper_endpoint, real_t atol
) {
    // both endpoints are zero, assume entire interval is zero
    if (nearly_zero(f1) && nearly_zero(f2)) {
        return 0.0;
    }

    // bisect the interval if value changes between zero and non-zero on the
    // interval, or if we're near endpoints
    if (nearly_zero(f1) || nearly_zero(f2) || at_lower_endpoint || at_upper_endpoint) {
        // also only bisect if tolerance is not yet met
        if (!nearly_zero(x2 - x1, atol)) {
            real_t x12 = 0.5 * (x1 + x2);
            real_t f12 = f(x12);

            return bisect_and_integrate(x1, x12, f1, f12, f, at_lower_endpoint, false, atol) +
                   bisect_and_integrate(x12, x2, f12, f2, f, false, at_upper_endpoint, atol);
        }
    }

    // if at endpoint, use midpoint rule to avoid the integrable singularity on
    // the boundaries
    if (at_lower_endpoint || at_upper_endpoint) {
        return f((x1 + x2) / 2) * (x2 - x1);
    }

    return trapezoid(x1, x2, f1, f2);
}

/**
 * Adaptively integrate a function by bisecting intervals where
 * the integrand changes between zero and non-zero.
 * This outperforms general-purpose adaptive quadratures (like GSL)
 * for this application because our integrand has discontinuous support.
 *
 * xmin: lower integration limit
 * xmax: upper integration limit
 * dx: max step at which to sample the integrand
 *     (will be refined where function jumps from/to 0)
 * integrand: the function to be integrated
 * par: parameters passed to the func for its evaluation
 * xtol: minimum interval size for adaptive bisection
 */
template <typename Integrand>
real_t integrate_adaptive(real_t xmin, real_t xmax, real_t dx, Integrand &&f, real_t xtol) {
    real_t x1 = xmin;
    real_t x2 = x1;

    real_t f1 = f(x1);
    real_t integral = 0.0;

    real_t dx_zero = dx / 3;
    while (x2 < xmax) {
        if (nearly_zero(f1)) {
            x2 = x1 + dx_zero;
        } else {
            x2 = x1 + dx;
        }

        if (x2 > xmax) {
            x2 = xmax;
        }

        real_t f2 = f(x2);

        bool at_lower_endpoint = fabs(x1 - xmin) < RELAXED_EPS;
        bool at_upper_endpoint = fabs(x2 - xmax) < RELAXED_EPS;
        integral +=
            bisect_and_integrate(x1, x2, f1, f2, f, at_lower_endpoint, at_upper_endpoint, xtol);

        x1 = x2;
        f1 = f2;
    }

    return integral;
}

#include <type_traits>
// Use a GSL QUADPACK adaptive quadrature to integrate the function.
template <typename Integrand>
real_t integrate_gsl(real_t a, real_t b, Integrand &&f, real_t atol) {
    using F = std::decay_t<Integrand>;

    struct Adapter {
        static double eval(double x, void *p) { return (*static_cast<F *>(p))(x); }
    };

    F f_store = std::forward<Integrand>(f);

    gsl_function G;
    G.function = &Adapter::eval;
    G.params = &f_store;

    // relatively large max number of subdivision to
    // accomodate the discontinuous support of the integrand.
    len_t limit = 1000;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(limit);

    double result, error;
    real_t pts[2] = {a, b};
    gsl_integration_qagp(&G, pts, 2, atol, 0.0, limit, w, &result, &error);

    gsl_integration_workspace_free(w);
    return result;
}

// Integrate the provided function on x in [xmin, xmax] with uniform spacing over n_points.
template <typename Integrand>
real_t integrate_midpoint(real_t xmin, real_t xmax, len_t n_points, Integrand &&f) {
    real_t dx = (xmax - xmin) / (real_t)n_points;
    real_t integral = 0;
    real_t x = xmin + 0.5 * dx;
    for (len_t i = 0; i < n_points; i++) {
        integral += dx * f(x);
        x += dx;
    }
    return integral;
}

}  // namespace

/**
 * Calculation of the fundamental \delta_j function in the notes.
 *
 * Performs the orbit (theta) integral of the analytically FVM- and gyro-averaged
 * delta function. Requires the caller to provide a bounding interval (theta1,
 * theta2) outside which the integrand is 0. This function will set the
 * integrand to 0 wherever the kinematic constraints are not satisfied, so
 * providing a good (theta1, theta2) guess helps with convergence
 * and computation time.
 *
 * The caller is also responsible for providing the BA volume element Vp at
 * xi01. NOTE: Vp ought to be normalized by p1^2, e.g. should have no momentum
 * dependence.
 */
real_t KnockOnUtilities::EvaluateOrbitAveragedDelta(
    len_t ir, real_t xi_star, real_t xi01, real_t xi0_f1, real_t xi0_f2, real_t Vp1, real_t theta1,
    real_t theta2, const FVM::RadialGrid *rg, len_t n_points_integral,
    orbit_integration_method quad_method
) {
    if (fabs(xi01) < RELAXED_EPS) {
        // when xi01 == 0, the function reduced to 0/0, but the correct limit is known:
        real_t BOverBmin = 1;
        real_t xi0Cutoff = 0;
        real_t xi1_over_xi01;
        return AnalyticDelta::EvaluateLocalContribution(
            xi0_f1, xi0_f2, xi01, BOverBmin, xi_star, xi0Cutoff, xi1_over_xi01
        );
    }

    auto integrand = [&](real_t theta) {
        return KnockOnUtilities::OrbitDeltaIntegrand(theta, ir, xi0_f1, xi0_f2, xi01, xi_star, rg);
    };
    real_t integral;
    if (quad_method == ADAPTIVE_TRAPEZOID) {
        real_t dtheta = (theta2 - theta1) / (real_t)n_points_integral;
        integral = integrate_adaptive(theta1, theta2, dtheta, integrand, QUAD_ADAPTIVE_ATOL);
    } else if (quad_method == ADAPTIVE_GSL) {
        integral = integrate_gsl(theta1, theta2, integrand, QUAD_ADAPTIVE_ATOL);
    } else {
        integral = integrate_midpoint(theta1, theta2, n_points_integral, integrand);
    }

    // Note, here only one factor of (2*pi) instead of (2*pi)^2 as in the notes,
    // because here the metric sqrtG (and therefore also Vp) contains an extra
    // factor of 2*pi.
    return 2 * M_PI * integral / Vp1;
}

/**
 * Compiles the true delta with trapping corrections.
 *
 * Here we will use the fortunate symmetry property of the Delta calculation,
 * that under
 *     xi_star -> -xi_star; and xi01 -> -xi0,
 * the value is unchanged.
 * Now, when [xi0_f1, xi0_f2] is in the trapped region, we should add the
 * knock on integral from the mirrored xi -> -xi contributions.
 * This happens to correspond to the case xi_star -> -xi_star.
 * Similarly, when xi01 is in the trapped region, we should add
 * xi1 -> -xi1. Because of the symmetry property, this is identical
 * to the other case.
 * Finally, when both xi01 and xi0 are in the trapped region, we should add the
 * double-flipped xi -> -xi; xi1 -> -xi1, which is the same as the original
 * calculation.
 *
 * Corresponds to \delta_j from the theory notes.
 */
real_t KnockOnUtilities::EvaluateOrbitAveragedDeltaWithTrappingCorrection(
    len_t ir, real_t xi_star, real_t xi01, real_t xi0_f1, real_t xi0_f2, real_t Vp1, real_t theta1,
    real_t theta2, const FVM::RadialGrid *rg, len_t n_points_integral, orbit_integration_method quad
) {
    // no contribution for unphysical parts of the xi01 phase space,
    // i.e. the negative trapped region.
    if (nearly_zero(Vp1)) {
        return 0.0;
    }

    real_t xi0T = rg->GetXi0TrappedBoundary(ir);

    // First, split contributions when cells cross trapped-passing boundary or
    // 0. We want to evaluate with cell faces on special points to do mirroring
    // properly

    // split out negative passing contribution
    if (xi0_f1 < -xi0T - RELAXED_EPS && xi0_f2 > -xi0T + RELAXED_EPS) {
        return (
            EvaluateOrbitAveragedDeltaWithTrappingCorrection(
                ir, xi_star, xi01, xi0_f1, -xi0T, Vp1, theta1, theta2, rg, n_points_integral, quad
            ) +
            EvaluateOrbitAveragedDeltaWithTrappingCorrection(
                ir, xi_star, xi01, -xi0T, xi0_f2, Vp1, theta1, theta2, rg, n_points_integral, quad
            )
        );
    }

    // handle when 0 is inside the cell, discard negative-trapped region
    if (xi0_f1 < -RELAXED_EPS && xi0_f2 > RELAXED_EPS) {
        return EvaluateOrbitAveragedDeltaWithTrappingCorrection(
            ir, xi_star, xi01, 0, xi0_f2, Vp1, theta1, theta2, rg, n_points_integral, quad
        );
    }

    // escape pure negative-trapped calculation
    if (xi0_f1 > -xi0T - RELAXED_EPS && xi0_f2 < RELAXED_EPS) {
        return 0;
    }

    // split out positive passing region
    if (xi0_f1 < xi0T - RELAXED_EPS && xi0_f2 > xi0T + RELAXED_EPS) {
        return (
            EvaluateOrbitAveragedDeltaWithTrappingCorrection(
                ir, xi_star, xi01, xi0_f1, xi0T, Vp1, theta1, theta2, rg, n_points_integral, quad
            ) +
            EvaluateOrbitAveragedDeltaWithTrappingCorrection(
                ir, xi_star, xi01, xi0T, xi0_f2, Vp1, theta1, theta2, rg, n_points_integral, quad
            )
        );
    }

    // Now, we know that [xi0_f1, xi0_f2] is either purely passing or purely in
    // the positive trapped region. Here we can add the unit contributions from
    // a single bounce average integral.

    real_t delta = EvaluateOrbitAveragedDelta(
        ir, xi_star, xi01, xi0_f1, xi0_f2, Vp1, theta1, theta2, rg, n_points_integral, quad
    );
    bool isTrapped = (xi0_f1 > -RELAXED_EPS) && (xi0_f2 < xi0T + RELAXED_EPS);
    bool is1Trapped = xi01 > -RELAXED_EPS && xi01 < xi0T;
    if (isTrapped || is1Trapped) {
        delta += EvaluateOrbitAveragedDelta(
            ir, -xi_star, xi01, xi0_f1, xi0_f2, Vp1, theta1, theta2, rg, n_points_integral, quad
        );
    }
    if (isTrapped && is1Trapped) {
        delta *= 2;
    }
    return delta;
}

/**
 * Estimates a conservative poloidal angle interval [theta1, theta2]
 * outside which the orbit-averaged delta integrand is guaranteed to vanish
 * for the given momentum-space cell indices (j, l).
 */
void KnockOnUtilities::EstimateBoundingTheta(
    len_t ir, len_t j, len_t l, real_t &theta1, real_t &theta2, const FVM::Grid *grid_knockon,
    const FVM::Grid *grid_primary
) {
    if (!grid_knockon->HasTrapped()) {
        // cylindrical case, in which case thetabounce have not been initialized
        theta1 = 0;
        theta2 = 2 * M_PI;
        return;
    }
    real_t xi0_f2 = grid_knockon->GetMomentumGrid(ir)->GetP2_f(j + 1);
    if (xi0_f2 <= RELAXED_EPS) {
        theta1 = grid_knockon->GetThetaBounce1_f2(ir, 0, j);
        theta2 = grid_knockon->GetThetaBounce2_f2(ir, 0, j);
    } else {
        theta1 = grid_knockon->GetThetaBounce1_f2(ir, 0, j + 1);
        theta2 = grid_knockon->GetThetaBounce2_f2(ir, 0, j + 1);
    }
    theta1 = std::max(theta1, grid_primary->GetThetaBounce1(ir, 0, l));
    theta2 = std::min(theta2, grid_primary->GetThetaBounce2(ir, 0, l));
    // we nudge the endpoints into the interval to ensure that we won't
    // be culled due to roundoff
    real_t theta_range = theta2 - theta1;
    constexpr real_t roundoff_safety_factor = 1e-8;
    theta1 = theta1 + roundoff_safety_factor * theta_range;
    theta2 = theta2 - roundoff_safety_factor * theta_range;
}

/**
 * Evaluates the orbit-averaged Delta including trapping corrections on
 * the provided grid.
 *
 * Uses xi01 on grid centers at index l and xi0_{j-1/2}, xi0_{j+1/2} at index j.
 * This function is a grid-aware convenience wrapper around
 * EvaluateOrbitAveragedDeltaWithTrappingCorrection().
 *
 * Corresponds to the matrix element Delta_{j l} in the theory notes.
 */
real_t KnockOnUtilities::EvaluateDeltaMatrixElementOnGrid(
    len_t ir, real_t xi_star, len_t j, len_t l, const FVM::Grid *grid_knockon,
    const FVM::Grid *grid_primary, len_t n_points_integral, orbit_integration_method quad
) {
    real_t xi01 = grid_primary->GetMomentumGrid(ir)->GetXi0(0, l);
    FVM::MomentumGrid *mg = grid_knockon->GetMomentumGrid(ir);
    real_t xi0_f1 = mg->GetXi0_f2(0, j);
    real_t xi0_f2 = mg->GetXi0_f2(0, j + 1);
    real_t Vp1 = grid_primary->GetVpOverP2AtZero(ir)[l];

    FVM::RadialGrid *rGrid = grid_primary->GetRadialGrid();

    // Handle the case in which xi01 is near the trapped-passing boundary.
    // In the case that xi01 is exactly on the boundary, the bounce average diverges.
    // For Vp this is solved by averaging the bounce integral over the entire cell--
    // that seems unnecessarily expensive here. Let's just evaluate on the passing side
    // of the TP boundary. This will cause the delta not to be conservative at this xi01,
    // but this error can be normalized away.
    real_t xi01_f1 = grid_primary->GetMomentumGrid(ir)->GetXi0_f2(0, l);
    real_t xi01_f2 = grid_primary->GetMomentumGrid(ir)->GetXi0_f2(0, l + 1);
    if (rGrid->GetFluxSurfaceAverager()->shouldCellAverageBounceIntegral(
            ir, xi01_f1, xi01_f2, DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION
        )) {
        if (xi01_f2 > rGrid->GetXi0TrappedBoundary(ir)) {
            xi01 = xi01_f2;
        } else {
            xi01 = xi01_f1;
        }
    }

    // Set the outer theta interval as the smallest of the orbits reached by xi_{01} and
    // [xi_{0,j-1/2}, xi_{0, j+1/2}]
    real_t theta1, theta2;
    EstimateBoundingTheta(ir, j, l, theta1, theta2, grid_knockon, grid_primary);
    real_t delta = EvaluateOrbitAveragedDeltaWithTrappingCorrection(
        ir, xi_star, xi01, xi0_f1, xi0_f2, Vp1, theta1, theta2, rGrid, n_points_integral, quad
    );
    return delta;
}

/**
 * Evaluates and normalizes a column of the orbit-averaged knock-on delta
 * function on the given FVM grids.
 *
 * For a fixed primary pitch index l (xi01) and xi_star, this routine
 * computes Delta_{j l} for all secondary pitch cells j using
 * EvaluateDeltaMatrixElementOnGrid(), and enforces the exact analytic
 * normalization property
 *
 *     sum_j Delta_{j l} * dxi0_j = 1,
 *
 * whenever a non-zero contribution is physically expected.
 *
 * If the initial quadrature resolution fails to capture any non-zero
 * contribution (e.g. due to very narrow contributing poloidal-angle
 * intervals), the number of poloidal integration points is increased
 * adaptively until either a non-zero contribution is found or a
 * hard upper limit is reached.
 *
 * This function encodes physics-domain knowledge about the definition
 * and normalization of the orbit-averaged knock-on delta function, and
 * is intended to be used by collision operators when assembling the
 * knock-on kernel.
 *
 * Returns the unnormalized column integral sum_j Delta_{j l} * dxi0_j.
 */
real_t KnockOnUtilities::SetDeltaMatrixColumnOnGrid(
    len_t ir, real_t xi_star, len_t l, const FVM::Grid *grid_knockon, const FVM::Grid *grid_primary,
    real_t *deltaCol, len_t n_points_integral, orbit_integration_method quad
) {
    if (grid_primary->GetVpOverP2AtZero(ir)[l] == 0) {
        return 0;
    }
    FVM::MomentumGrid *mg = grid_knockon->GetMomentumGrid(ir);
    len_t Nxi = mg->GetNp2();

    real_t norm = 0;
    for (len_t j = 0; j < Nxi; j++) {
        real_t dxi0j = mg->GetDp2(j);
        deltaCol[j] = EvaluateDeltaMatrixElementOnGrid(
            ir, xi_star, j, l, grid_knockon, grid_primary, n_points_integral, quad
        );
        norm += dxi0j * deltaCol[j];
    }

    // If no non-zero elements were discovered, but we expected them to,
    // it may mean that only very small poloidal-angle intervals contributed
    // to the orbit integrals, and were not sampled at the current resolution.
    // This should be exceedingly rare.
    // We expect non-zeroes when Vp(xi01) != 0. When == 0, this represents
    // negative-trapped orbits that should not be counted.
    bool nonzero_expected = !nearly_zero(grid_primary->GetVpOverP2AtZero(ir)[l]);
    if (norm < RELAXED_EPS && nonzero_expected) {
        // increase number of poloidal samples and try again, but terminate if 1000 is reached.
        constexpr len_t n_points_max = 1000;
        if (n_points_integral >= n_points_max) {
            throw DREAMException(
                "Knock-on operator expected non-zero elements at ir=%ld, l=%ld, xi_star=%.8g but "
                "none were found",
                ir, l, xi_star
            );
        }
        // try again with increased resolution in theta quadrature
        return SetDeltaMatrixColumnOnGrid(
            ir, xi_star, l, grid_knockon, grid_primary, deltaCol, n_points_integral * 4, quad
        );
    }
    for (len_t j = 0; j < Nxi; j++) {
        deltaCol[j] = deltaCol[j] / norm;
    }
    return norm;
}

/**
 * Evaluates the momentum-space knock-on kernel S_{ik} on a DREAM grid.
 *
 * Assumes a p-xi grid that is identical at all radii.
 *
 * Uses p1 on primary grid centers at index k and p_{i-1/2}, p_{i+1/2}
 * on the knock-on momentum grid at index i.
 *
 * This function is a grid-aware convenience wrapper around
 * Kinematics::EvaluateMollerFlux(), including the primary speed factor v1,
 * with kinematic clamping to the relativistically allowed knock-on domain.
 *
 * Corresponds to the matrix element S_{ik} in the theory notes.
 */
real_t KnockOnUtilities::EvaluateMollerFluxMatrixElementOnGrid(
    len_t i, len_t k, const FVM::Grid *grid_knockon, const FVM::Grid *grid_primary, real_t pCutoff
) {
    if (!(pCutoff > RELAXED_EPS)) {
        throw DREAMException(
            "EvaluateMollerFluxMatrixElementOnGrid(): invalid pCutoff=%.16g (must satisfy pCutoff "
            "> %.3g). "
            "This is required because the Møller flux antiderivative is singular at p=0.",
            pCutoff, RELAXED_EPS
        );
    }

    real_t p_f1 = grid_knockon->GetMomentumGrid(0)->GetP1_f(i);
    real_t p_f2 = grid_knockon->GetMomentumGrid(0)->GetP1_f(i + 1);

    if (p_f2 <= pCutoff) {
        // entire momentum cell below the admissible interval
        return 0.0;
    }

    real_t p1 = grid_primary->GetMomentumGrid(0)->GetP1(k);
    real_t pMax = KnockOnUtilities::Kinematics::MaximumKnockOnMomentum(p1);

    if (p_f1 > pMax) {
        // entire momentum cell above the admissible interval
        return 0.0;
    }

    // Clamp range to [pCutoff, pMax]
    real_t pLower = std::max(p_f1, pCutoff);
    real_t pUpper = std::min(p_f2, pMax);
    if (!(pUpper > pLower + RELAXED_EPS)) {
        return 0.0;
    }

    real_t S2 = KnockOnUtilities::Kinematics::EvaluateMollerFlux(pUpper, p1);
    real_t S1 = KnockOnUtilities::Kinematics::EvaluateMollerFlux(pLower, p1);
    real_t v1 = p1 / sqrt(1 + p1 * p1);
    return v1 * (S2 - S1) / (p_f2 - p_f1);
}

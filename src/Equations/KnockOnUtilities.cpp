#include "DREAM/Equations/KnockOnUtilities.hpp"
#include <cmath>

using namespace DREAM;

auto KnockOnUtilities::evaluateF(real_t y, real_t xi1, real_t xi_star) -> real_t {
    return asin((y - xi1 * xi_star) / sqrt((1 - xi1 * xi1) * (1 - xi_star * xi_star))) / M_PI;
}

auto KnockOnUtilities::evaluateD(real_t t1, real_t t2, real_t xi1, real_t xi_star) -> real_t {
    real_t a = xi_star * xi1;
    real_t b = sqrt(1 - xi_star * xi_star) * sqrt(1 - xi1 * xi1);
    real_t F1;
    real_t F2;
    if (t2 >= a + b) {
        F2 = 0.5;
    } else {
        F2 = evaluateF(t2, xi1, xi_star);
    }
    if (t1 <= a - b) {
        F1 = -0.5;
    } else {
        F1 = evaluateF(t1, xi1, xi_star);
    }

    return F2 - F1;
}

/**
 * Evaluate and set the `t1` and `t2` variables defined as in
 * docs/notes/theory.tex.
 */
void KnockOnUtilities::evaluateT1T2(
    real_t &t1, real_t &t2, real_t xi0_f1, real_t xi0_f2, real_t BOverBmin,
    real_t xi0Cutoff
) {
    if (xi0_f2 < 0) {
        // negative xi0 intervals:
        t1 = FVM::MomentumGrid::XiFromXi0(xi0_f1, BOverBmin);
        if (xi0_f2 < -xi0Cutoff) {
            t2 = FVM::MomentumGrid::XiFromXi0(xi0_f2, BOverBmin);
        } else {
            t2 = 0;
        }
    } else {
        // positive xi0 intervals:
        if (xi0_f1 > xi0Cutoff) {
            t1 = FVM::MomentumGrid::XiFromXi0(xi0_f1, BOverBmin);
        } else {
            t1 = 0;
        }
        t2 = FVM::MomentumGrid::XiFromXi0(xi0_f2, BOverBmin);
    }
}

/**
 * Evaluate and set the `z1` and `z2` variables defined as in
 * docs/notes/theory.tex. Note, that z essentially represents xi1, and [z1, z2]
 * is the range of values xi1 takes when xi0 is swept over [xi_{0,j-1/2},
 * xi_{0,j+1/2}].
 */
void KnockOnUtilities::evaluateZ1Z2(real_t &z1, real_t &z2, real_t t1, real_t t2, real_t xi_star) {
#define Z_upper(x) (x * xi_star + sqrt(1 - (x) * (x)) * sqrt(1 - xi_star * xi_star))
#define Z_lower(x) (x * xi_star - sqrt(1 - (x) * (x)) * sqrt(1 - xi_star * xi_star))
    if (t2 < -xi_star) {
        z1 = Z_lower(t2);
    } else if (t1 <= -xi_star && -xi_star <= t2) {
        z1 = -1;
    } else {
        z1 = Z_lower(t1);
    }

    if (t2 < xi_star) {
        z2 = Z_upper(t2);
    } else if (t1 <= xi_star && xi_star <= t2) {
        z2 = 1;
    } else {
        z2 = Z_upper(t1);
    }
}

real_t KnockOnUtilities::deltaIntegrand(
    real_t theta, len_t ir, real_t xi0_f1, real_t xi0_f2, real_t xi01, real_t xi_star,
    const FVM::Grid *grid
) {
    real_t eps = 100 * std::numeric_limits<real_t>::epsilon();
    FVM::RadialGrid *rg = grid->GetRadialGrid();
    FVM::FluxSurfaceAverager *fsa = rg->GetFluxSurfaceAverager();

    real_t Bmin = rg->GetBmin(ir);

    real_t B;
    real_t Jacobian;
    real_t xi0Cutoff;
    real_t ROverR0_;
    real_t NablaR2_;
    fsa->GeometricQuantitiesAtTheta(
        ir, theta, B, Jacobian, ROverR0_, NablaR2_, FVM::FLUXGRIDTYPE_DISTRIBUTION
    );
    real_t BOverBmin = B / Bmin;
    xi0Cutoff = sqrt(1 - 1 / BOverBmin);

    // Ensure that xi01 represent an orbit that
    // can reach this poloidal location. Otherwise, the function
    // has been sampled at a theta that shouldn't contribute => return 0.
    if (xi0Cutoff >= fabs(xi01)) {
        return 0.0;
    }

    real_t t1;
    real_t t2;
    evaluateT1T2(t1, t2, xi0_f1, xi0_f2, BOverBmin, xi0Cutoff);

    // perform the checks on whether the theta is viable
    // the "true" condition is z1 <= xi1 <= z2 for us to
    // take this theta into account. However, roundoff
    // errors when xi1 == z1, z2 can cause invalid values
    // to slip through and crash other parts of the calculation.
    // For the important edge case xi01 == +/- 1,
    // which define the Chiu-Harvey and Rosenbluth-Putvinski
    // approximations, we seek these bounding values, and therefore
    // employ an exacter constraint by inserting values explicitly.
    real_t Xi1OverXi01;
    real_t D;
    if (fabs(fabs(xi01) - 1) < eps) {
        // exact solution of z2==1 (since z = xi1 = 1)
        if (xi01 > 0 && !(t1 < xi_star && xi_star < t2)) {
            return 0;
        }
        // exact solution of z1==-1 (since z = xi1 = -1)
        if (xi01 < 0 && !(t1 < -xi_star && -xi_star < t2)) {
            return 0;
        }
        Xi1OverXi01 = 1;
        D = 1;
    } else {
        real_t z1;
        real_t z2;
        evaluateZ1Z2(z1, z2, t1, t2, xi_star);
        real_t xi1 = FVM::MomentumGrid::XiFromXi0(xi01, BOverBmin);
        Xi1OverXi01 = xi1 / xi01;
        if (!(z1 + eps < xi1 && xi1 < z2 - eps)) {
            return 0.0;
        }

        D = evaluateD(t1, t2, xi1, xi_star);
    }
    // from the integration limits on xi1 before we (analytically) exchanged
    // integration order with theta (which cannot be inverted analytically).
    real_t sqrtG = Jacobian * FVM::MomentumGrid::evaluatePXiMetricOverP2(Xi1OverXi01, BOverBmin);

    return sqrtG * D;
}

// Collect functions for an adaptive quadrature of the orbit averaged delta
// function
namespace {

// check if value if effectively zero
bool nearly_zero(real_t y, real_t tol = 1e-8) { return fabs(y) < tol; }

real_t trapezoid(real_t x1, real_t x2, real_t f1, real_t f2) { return (x2 - x1) * (f1 + f2) * 0.5; }

// Adaptively integrate a given interval with bisection (if needed) + trapz
real_t bisect_and_integrate(
    real_t x1, real_t x2, real_t f1, real_t f2, real_t (*integrand)(real_t, void *), void *par,
    bool at_lower_endpoint, bool at_upper_endpoint, real_t atol
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
            real_t f12 = integrand(x12, par);

            return bisect_and_integrate(
                       x1, x12, f1, f12, integrand, par, at_lower_endpoint, false, atol
                   ) +
                   bisect_and_integrate(
                       x12, x2, f12, f2, integrand, par, false, at_upper_endpoint, atol
                   );
        }
    }

    // if at endpoint, use midpoint rule to avoid the integrable singularity on
    // the boundaries
    if (at_lower_endpoint || at_upper_endpoint) {
        return integrand((x1 + x2) / 2, par) * (x2 - x1);
    }

    return trapezoid(x1, x2, f1, f2);
}

/**
 * Adaptively integrate a function by bisecting intervals where
 * the integrand changes between zero and non-zero.
 *
 * xmin: lower integration limit
 * xmax: upper integration limit
 * dx: max step at which to sample the integrand
 *     (will be refined where function jumps from/to 0)
 * integrand: the function to be integrated
 * par: parameters passed to the func for its evaluation
 */
real_t integrate_adaptive(
    real_t xmin, real_t xmax, real_t dx, real_t (*integrand)(real_t, void *), void *par,
    real_t xtol = 1e-3
) {
    real_t x1 = xmin;
    real_t x2 = x1;

    real_t f1 = integrand(x1, par);
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

        real_t f2 = integrand(x2, par);

        bool at_lower_endpoint = fabs(x1 - xmin) < 1e-8;
        bool at_upper_endpoint = fabs(x2 - xmax) < 1e-8;
        integral += bisect_and_integrate(
            x1, x2, f1, f2, integrand, par, at_lower_endpoint, at_upper_endpoint, xtol
        );

        x1 = x2;
        f1 = f2;
    }

    return integral;
}

real_t integrate_midpoint(
    real_t xmin, real_t xmax, len_t n_points, real_t (*integrand)(real_t, void *), void *par
) {
    real_t dx = (xmax - xmin) / n_points;
    real_t integral = 0;
    real_t x = xmin + 0.5 * dx;
    for (len_t i = 0; i < n_points; i++) {
        integral += dx * integrand(x, par);
        x += dx;
    }
    return integral;
}

struct DeltaIntegrandParameters {
    len_t ir;
    real_t xi0_f1;
    real_t xi0_f2;
    real_t xi01;
    real_t xi_star;
    const FVM::Grid *grid;
};
real_t delta_integrand(real_t theta, void *deltaParams) {
    DeltaIntegrandParameters *params = (struct DeltaIntegrandParameters *)deltaParams;
    return KnockOnUtilities::deltaIntegrand(
        theta, params->ir, params->xi0_f1, params->xi0_f2, params->xi01, params->xi_star,
        params->grid
    );
}
} // namespace

/**
 * Performs the orbit (theta) integral of the analytically FVM- and gyroaveraged
 * delta function. Requires the caller to provide a bounding interval (theta1,
 * theta2) outside which the integrand is 0. This function will set the
 * integrand to 0 wherever the kinematic constraints are not satisfied, so
 * providing a good (theta1, theta2) guess just helps with convergence.
 *
 * The caller is also responsible for providing the BA volume element Vp at
 * xi01. NOTE: Vp ought to be normalized by p1^2, e.g. should have no momentum
 * dependence.
 */
real_t KnockOnUtilities::EvaluateDeltaContribution(
    len_t ir, real_t xi_star, real_t xi01, real_t xi0_f1, real_t xi0_f2, real_t Vp1, real_t theta1,
    real_t theta2, len_t n_points_integral, const FVM::Grid *grid, KnockOnUtilities::orbitIntegrationQuadrature quadrature
) {

    DeltaIntegrandParameters par = {ir, xi0_f1, xi0_f2, xi01, xi_star, grid};
    
    real_t integral;
    if (quadrature == QUADRATURE_ADAPTIVE){
        real_t dtheta = (theta2 - theta1) / n_points_integral;
        integral = integrate_adaptive(theta1, theta2, dtheta,&(delta_integrand), &par, 0.005);
    } else {
        integral = integrate_midpoint(theta1, theta2, n_points_integral, &(delta_integrand), &par);
    }
    // Note, here only one factor of (2*pi) instead of (2*pi)^2 as in the notes,
    // because here the metric sqrtG (and therefore also Vp) contains an extra
    // factor of 2*pi.
    real_t dxi = xi0_f2 - xi0_f1;

    return 2 * M_PI * integral / (Vp1 * dxi);
}

/**
 * Compiles the true delta with trapping corrections.
 *
 * Here we will use
 * the fortunate symmetry property of the Delta calculation, that under
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
 */
real_t KnockOnUtilities::EvaluateDelta(
    len_t ir, real_t xi_star, real_t xi01, real_t xi0_f1, real_t xi0_f2, real_t Vp1, real_t theta1,
    real_t theta2, len_t n_points_integral, const FVM::Grid *grid
) {
    real_t xi0T = grid->GetRadialGrid()->GetXi0TrappedBoundary(ir);
    real_t eps = 100 * std::numeric_limits<real_t>::epsilon();

    // First, split contributions when cells cross trapped-passing boundary or
    // 0. We want to evaluate with cell faces on special points to do mirroring
    // properly

    // split out negative passing contribution
    if (xi0_f1 < -xi0T - eps && xi0_f2 > -xi0T + eps) {
        return (
            EvaluateDelta(
                ir, xi_star, xi01, xi0_f1, -xi0T, Vp1, theta1, theta2, n_points_integral, grid
            ) +
            EvaluateDelta(
                ir, xi_star, xi01, -xi0T, xi0_f2, Vp1, theta1, theta2, n_points_integral, grid
            )
        );
    }

    // handle when 0 is inside the cell
    if (xi0_f1 < -eps && xi0_f2 > eps) {
        return EvaluateDelta(
            ir, xi_star, xi01, 0, xi0_f2, Vp1, theta1, theta2, n_points_integral, grid
        );
    }

    // escape pure negative-trapped calculation
    if (xi0_f1 > -xi0T - eps && xi0_f2 < eps) {
        return 0;
    }

    // split out positive passing region
    if (xi0_f1 < xi0T - eps && xi0_f2 > xi0T + eps) {
        return (
            EvaluateDelta(
                ir, xi_star, xi01, xi0_f1, xi0T, Vp1, theta1, theta2, n_points_integral, grid
            ) +
            EvaluateDelta(
                ir, xi_star, xi01, xi0T, xi0_f2, Vp1, theta1, theta2, n_points_integral, grid
            )
        );
    }

    // Now, we know that [xi0_f1, xi0_f2] is either purely passing or purely in
    // the positive trapped region

    real_t delta = KnockOnUtilities::EvaluateDeltaContribution(
        ir, xi_star, xi01, xi0_f1, xi0_f2, Vp1, theta1, theta2, n_points_integral, grid
    );
    bool isTrapped = (xi0_f1 > -eps) && (xi0_f2 < xi0T + eps);
    bool is1Trapped = xi01 > -eps && xi01 < xi0T;
    if (isTrapped || is1Trapped) {
        delta += KnockOnUtilities::EvaluateDeltaContribution(
            ir, -xi_star, xi01, xi0_f1, xi0_f2, Vp1, theta1, theta2, n_points_integral, grid
        );
    }
    if (isTrapped && is1Trapped) {
        delta *= 2;
    }
    return delta;
}
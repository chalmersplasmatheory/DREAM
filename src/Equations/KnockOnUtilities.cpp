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

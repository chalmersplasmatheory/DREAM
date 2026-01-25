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

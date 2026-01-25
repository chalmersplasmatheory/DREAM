/**
 * Implements unit tests of the KnockOnUtilities module.
 *
 * Where numerical accuracy is tested, we intentionally allow relatively loose
 * tolerances. This is partly to limit computational cost, but mainly because
 * this level of accuracy is sufficient for production use.
 *
 * In typical simulations, relatively low numerical precision is acceptable
 * since exact particle-number conservation is enforced by construction in
 * the knock-on operator. We primarily need the shapes to be sufficiently resolved.
 */
#include "KnockOn.hpp"

#include "DREAM/Equations/KnockOnUtilities.hpp"
#include "FVM/Grid/Grid.hpp"

using namespace DREAMTESTS::_DREAM;

/**
 * Run this test.
 */
bool KnockOn::Run(bool) {
    bool success = true;

    if (CheckDeltaMirrorProperties())
        this->PrintOK("The delta function has appropriate symmetries in xi and xi1.");
    else {
        success = false;
        this->PrintError(
            "Knock-on test failed: delta function does not have the expected symmetries."
        );
    }

    if (CheckDeltaQuadratureConvergence())
        this->PrintOK("The default delta-function quadrature is sufficiently accurate.");
    else {
        success = false;
        this->PrintError(
            "Knock-on test failed: delta-function quadrature is not sufficiently accurate."
        );
    }

    if (CheckLocalContributionConservationProperty())
        this->PrintOK("The local contribution integrated over xi0 correctly sums to 1.");
    else {
        success = false;
        this->PrintError(
            "Knock-on test failed: The local contribution integrated over xi0 does not sum to 1."
        );
    }

    if (CheckDeltaConservationProperty())
        this->PrintOK(
            "The delta function integrated over xi0 is sufficiently close to 1 when primaries and "
            "secondaries live on the same grid."
        );
    else {
        success = false;
        this->PrintError(
            "Knock-on test failed: the delta-function does not integrate to 1 when primaries and "
            "secondaries live on the same grid."
        );
    }

    if (CheckDeltaConservationPropertyDifferentGrids())
        this->PrintOK(
            "The delta function integrated over xi0 is sufficiently close to 1 when primaries and "
            "secondaries live on different grids."
        );
    else {
        success = false;
        this->PrintError(
            "Knock-on test failed: the delta-function does not integrate to 1 when primaries and "
            "secondaries live on different grids."
        );
    }

    if (CheckAgreementWithOldRPTerm())
        this->PrintOK(
            "The new delta-function implementation agrees with the previous one for xi01=1"
        );
    else {
        success = false;
        this->PrintError(
            "Knock-on test failed: the new delta-function implementation does not "
            "agree with the old one."
        );
    }

    if (CheckCylindricalDeltaCalculation())
        this->PrintOK("The delta matrix elements agree exactly with the analytic result.");
    else {
        success = false;
        this->PrintError(
            "Knock-on test failed: the delta matrix elements do not "
            "agree exactly with theory in the cylindrical limit."
        );
    }

    if (CheckMollerFluxIntegration())
        this->PrintOK("The Moller cross section correctly integrates to the Moller flux.");
    else {
        success = false;
        this->PrintError("The Moller cross section fails to integrate to the Moller flux.");
    }

    if (CheckMollerDifferentialConvergesToInfiniteLimit())
        this->PrintOK("The Moller differential cross-section converges to the infinite limit.");
    else {
        success = false;
        this->PrintError(
            "The Moller differential cross-section does not converge to the infinite limit."
        );
    }

    if (CheckMollerFluxConvergesToInfiniteLimit())
        this->PrintOK("The Moller integrated flux converges to the infinite limit.");
    else {
        success = false;
        this->PrintError("The Moller integrated flux does not converge to the infinite limit.");
    }

    if (CheckMollerFluxConvergesToInfiniteLimit())
        this->PrintOK("The scattering angle xiStar converges to the infinite limit.");
    else {
        success = false;
        this->PrintError("The scattering angle xiStar does not converge to the infinite limit.");
    }

    return success;
}

namespace {}

// The orbit-averaged delta function analytically satisfies the
// exact relation delta(xi_star, xi01, ...) = delta(-xi_star, -xi01),
// representing the mirror symmetry of the collisions as well as the
// particle orbits in the magnetic field. This test verifies that this
// physical invariance is preserved numerically.
bool KnockOn::CheckDeltaMirrorProperties() {
    real_t successRelErrorThreshold = 1e-5;

    len_t nr = 5;
    len_t np = 4;
    len_t nxi = 30;

    len_t ntheta_interp = 100;
    len_t nrProfiles = 8;

    real_t pMin = 0;
    real_t pMax = 2;

    DREAM::FVM::Grid *grid =
        InitializeGridGeneralRPXi(nr, np, nxi, ntheta_interp, nrProfiles, pMin, pMax);
    DREAM::FVM::RadialGrid *rg = grid->GetRadialGrid();

    bool success = true;
    real_t xi_stars[4] = {0.1, 0.5, 0.9, 1.0};
    len_t xi0_indices[6] = {0, 5, 16, 20, 25, 29};
    len_t xi01_indices[6] = {0, 5, 16, 20, 25, 29};
    for (len_t ir = 0; ir < nr; ir++) {
        DREAM::FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir);
        for (len_t n = 0; n < 4; n++) {
            real_t xi_star = xi_stars[n];
            for (len_t idx_j = 0; idx_j < 6; idx_j++) {
                len_t j = xi0_indices[idx_j];
                if (grid->IsNegativePitchTrappedIgnorableCell(ir, j)) {
                    // everything should be 0, not interesting
                    continue;
                }
                real_t xi0_f1 = mg->GetP2_f(j);
                real_t xi0_f2 = mg->GetP2_f(j + 1);
                for (len_t idx_l = 0; idx_l < 6; idx_l++) {
                    len_t l = xi01_indices[idx_l];
                    real_t xi01 = mg->GetP2(l);
                    if (grid->IsNegativePitchTrappedIgnorableCell(ir, l)) {
                        // everything should be 0, not interesting
                        continue;
                    }
                    real_t Vp1 = grid->GetVpOverP2AtZero(ir)[l];
                    real_t theta1, theta2;
                    DREAM::KnockOnUtilities::EstimateBoundingTheta(
                        ir, j, l, theta1, theta2, grid, grid
                    );

                    real_t delta1 = DREAM::KnockOnUtilities::EvaluateOrbitAveragedDelta(
                        ir, xi_star, xi01, xi0_f1, xi0_f2, Vp1, theta1, theta2, rg, 30,
                        DREAM::KnockOnUtilities::MIDPOINT_RULE
                    );
                    real_t delta2 = DREAM::KnockOnUtilities::EvaluateOrbitAveragedDelta(
                        ir, -xi_star, -xi01, xi0_f1, xi0_f2, Vp1, theta1, theta2, rg, 30,
                        DREAM::KnockOnUtilities::MIDPOINT_RULE
                    );
                    if (fabs(delta1 - delta2) > fabs(delta1) * successRelErrorThreshold) {
                        this->PrintError(
                            "Mirror symmetry failed at ir=%ld, j=%ld, l=%ld, xi_star=%.3g\n", ir, j,
                            l, xi_star
                        );
                        success = false;
                    }
                    real_t delta3 = DREAM::KnockOnUtilities::EvaluateOrbitAveragedDelta(
                        ir, -xi_star, xi01, -xi0_f2, -xi0_f1, Vp1, theta1, theta2, rg, 30,
                        DREAM::KnockOnUtilities::MIDPOINT_RULE
                    );
                    if (fabs(delta2 - delta3) > fabs(delta2) * successRelErrorThreshold) {
                        this->PrintError(
                            "Mirror symmetry failed at ir=%ld, j=%ld, l=%ld, xi_star=%.3g\n", ir, j,
                            l, xi_star
                        );
                        success = false;
                    }
                }
            }
        }
    }
    delete grid;
    return success;
}

// Tests that the orbit-averaged delta calculation converges with increasing
// number of quadrature points. This is a numerical consistency test rather
// than a test of a physical invariant.
bool KnockOn::CheckDeltaQuadratureConvergence() {
    real_t successRelErrorThreshold = 0.1;

    len_t nr = 5;
    len_t np = 4;
    len_t nxi = 31;

    len_t ntheta_interp = 100;
    len_t nrProfiles = 8;

    real_t pMin = 0;
    real_t pMax = 2;

    DREAM::FVM::Grid *grid =
        InitializeGridGeneralRPXi(nr, np, nxi, ntheta_interp, nrProfiles, pMin, pMax);

    bool success = true;
    real_t xi_stars[4] = {0.1, 0.5, 0.9, 1.0};
    // len_t xi0_indices[6] = {0, 5, 16, 20, 25, 29};
    len_t xi01_indices[6] = {0, 5, 16, 20, 25, 29};
    // real_t xi01s[6] = {-1, -0.9, 0.02, 0.2, 0.8, 1.0};
    for (len_t ir = 0; ir < nr; ir++) {
        for (len_t n = 0; n < 4; n++) {
            real_t xi_star = xi_stars[n];
            for (len_t j = 0; j < grid->GetNp2(ir); j++) {
                if (grid->IsNegativePitchTrappedIgnorableCell(ir, j)) {
                    continue;
                }
                // len_t j = xi0_indices[idx_j];
                for (len_t idx_l = 0; idx_l < 6; idx_l++) {
                    len_t l = xi01_indices[idx_l];
                    real_t delta_default =
                        DREAM::KnockOnUtilities::EvaluateDeltaMatrixElementOnGrid(
                            ir, xi_star, j, l, grid, grid,
                            DREAM::KnockOnUtilities::N_POINTS_INTEGRAL_DEFAULT,
                            DREAM::KnockOnUtilities::MIDPOINT_RULE
                        );
                    real_t delta_hires = DREAM::KnockOnUtilities::EvaluateDeltaMatrixElementOnGrid(
                        ir, xi_star, j, l, grid, grid, 300, DREAM::KnockOnUtilities::MIDPOINT_RULE
                    );
                    if (fabs(delta_default - delta_hires) >
                        successRelErrorThreshold * (1 + fabs(delta_default))) {
                        this->PrintError("failed (j=%ld, l=%ld)\n", j, l);
                        this->PrintError("delta_default = %.8g\n", delta_default);
                        this->PrintError("delta_hires = %.8g\n", delta_hires);

                        success = false;
                    }
                }
            }
        }
    }
    delete grid;
    return success;
}

// check the conservation property that sum_j dxi_j * D_j == 1 for all combinations of pitches on
// the provided grids, and for a range of poloidal angles theta and values of xi_star.
bool KnockOn::_checkLocalContributionConservationProperty(
    const DREAM::FVM::Grid *grid_knockon, const DREAM::FVM::Grid *grid_primary,
    real_t successRelErrorThreshold
) {
    bool success = true;

    DREAM::FVM::RadialGrid *rg = grid_knockon->GetRadialGrid();
    real_t xi_stars[4] = {0.1, 0.5, 0.9, 1.0};
    len_t n_theta = 10;
    constexpr real_t theta_roundoff_safety_factor = 1e-8;

    for (len_t ir = 0; ir < rg->GetNr(); ir++) {
        const DREAM::FVM::MomentumGrid *mg = grid_primary->GetMomentumGrid(ir);

        len_t Nxi = grid_knockon->GetNp2(ir);
        len_t Nxi1 = grid_primary->GetNp2(ir);
        for (len_t n = 0; n < 4; n++) {
            real_t xi_star = xi_stars[n];
            for (len_t l = 0; l < Nxi1; l++) {
                real_t xi01 = mg->GetXi0(0, l);
                real_t theta1 =
                    grid_primary->GetThetaBounce1(ir, 0, l) + theta_roundoff_safety_factor;
                real_t theta2 =
                    grid_primary->GetThetaBounce2(ir, 0, l) - theta_roundoff_safety_factor;
                for (len_t m = 0; m < n_theta; m++) {
                    real_t theta = theta1 + (theta2 - theta1) * (real_t)m / (n_theta - 1);
                    real_t local_sum = 0;
                    for (len_t j = 0; j < Nxi; j++) {
                        real_t xi0_f1 = grid_knockon->GetMomentumGrid(ir)->GetXi0_f2(0, j);
                        real_t xi0_f2 = grid_knockon->GetMomentumGrid(ir)->GetXi0_f2(0, j + 1);
                        real_t BOverBmin;
                        real_t Jacobian;
                        real_t xi0Cutoff;
                        real_t xi1_over_xi01;
                        bool isReachable =
                            DREAM::KnockOnUtilities::CheckIfReachableAndSetGeometricQuantities(
                                ir, theta, xi0_f1, xi0_f2, xi01, BOverBmin, Jacobian, xi0Cutoff,
                                rg->GetFluxSurfaceAverager()
                            );
                        if (!isReachable) {
                            continue;
                        }
                        real_t dxi = (xi0_f2 - xi0_f1);
                        local_sum +=
                            dxi *
                            DREAM::KnockOnUtilities::AnalyticDelta::EvaluateLocalContribution(
                                xi0_f1, xi0_f2, xi01, BOverBmin, xi_star, xi0Cutoff, xi1_over_xi01
                            );
                    }
                    if (fabs(local_sum - 1.0) > successRelErrorThreshold) {
                        this->PrintError("failed at ir=%ld, n=%ld, l=%ld, m=%ld:\n", ir, n, l, m);
                        this->PrintError("  sum: %.8g", local_sum);
                        success = false;
                    }
                }
            }
        }
    }

    return success;
}

bool KnockOn::CheckLocalContributionConservationProperty() {
    real_t successRelErrorThreshold = 10 * std::numeric_limits<real_t>::epsilon();

    len_t nr = 5;
    len_t np = 4;
    len_t nxi = 60;

    len_t ntheta_interp = 50;
    len_t nrProfiles = 8;

    real_t pMin = 0;
    real_t pMax = 2;

    DREAM::FVM::Grid *grid =
        InitializeGridGeneralRPXi(nr, np, nxi, ntheta_interp, nrProfiles, pMin, pMax);

    bool success =
        _checkLocalContributionConservationProperty(grid, grid, successRelErrorThreshold);
    delete grid;
    return success;
}

// Return true if the delta matrix satisfies the exact analytic normalization property:
//     \int dxi0 delta(xi0, xi01) = 1
// within the provided tolerance, on the provided grids.
bool KnockOn::_checkDeltaConservationProperty(
    const DREAM::FVM::Grid *grid_knockon, const DREAM::FVM::Grid *grid_primary,
    real_t successRelErrorThreshold
) {
    DREAM::FVM::RadialGrid *rg = grid_knockon->GetRadialGrid();
    real_t xi_stars[4] = {0.1, 0.5, 0.9, 1.0};

    bool success = true;
    for (len_t ir = 0; ir < rg->GetNr(); ir++) {
        real_t xi0T = rg->GetXi0TrappedBoundary(ir);
        const DREAM::FVM::MomentumGrid *mg = grid_primary->GetMomentumGrid(ir);

        len_t Nxi = grid_knockon->GetNp2(ir);
        len_t Nxi1 = grid_primary->GetNp2(ir);
        real_t *_deltaCol = new real_t[Nxi];
        for (len_t n = 0; n < 4; n++) {
            real_t xi_star = xi_stars[n];
            for (len_t l = 0; l < Nxi1; l++) {
                real_t xi01_f1 = mg->GetP2_f(l);
                real_t xi01_f2 = mg->GetP2_f(l + 1);

                // Skip cells where the phase-space Jacobian Vp(xi01) is evaluated using
                // cell-averaged bounce integrals; in this case the delta normalization
                // is not expected to be exact. This is fine in the final collision operators,
                // where we correct the matrix elements to account for this deviation.
                if (rg->GetFluxSurfaceAverager()->shouldCellAverageBounceIntegral(
                        ir, xi01_f1, xi01_f2, DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION
                    )) {
                    continue;
                }
                // Also skip cells that contain 0 - these are also handled in a special way
                // in the FluxSurfaceAverager, since otherwise Vp==0 in this point.
                if (xi01_f1 < 0 && xi01_f2 > 0) {
                    continue;
                }
                // Continue if primary pitch is in a void FVM cell -- this represents the negative
                // trapped region, and will not contribute to the total knock on collision terms.
                if (grid_primary->GetVpOverP2AtZero(ir)[l] == 0) {
                    continue;
                }
                real_t deltaXiIntegral = DREAM::KnockOnUtilities::SetDeltaMatrixColumnOnGrid(
                    ir, xi_star, l, grid_knockon, grid_primary, _deltaCol,
                    DREAM::KnockOnUtilities::N_POINTS_INTEGRAL_DEFAULT,
                    DREAM::KnockOnUtilities::MIDPOINT_RULE
                );
                if (fabs(deltaXiIntegral - 1) > successRelErrorThreshold) {
                    real_t xi01 = mg->GetP2(l);
                    this->PrintError("failed at ir=%ld, n=%ld, l=%ld:\n", ir, n, l);
                    this->PrintError("  xi_star: %.8g\n", xi_star);
                    this->PrintError("  trapped boundary: %.8g\n", xi0T);
                    this->PrintError("  xi01: %.8g (%.8g, %.8g)\n", xi01, xi01_f1, xi01_f2);
                    this->PrintError("  deltaXiIntegral: %.8g\n", deltaXiIntegral);
                    success = false;
                }
            }
        }
        delete[] _deltaCol;
    }
    return success;
}

// The Delta function in our formulation of the knock-on operator satisfies
// an exact analytic normalization property:
//     \int dxi0 delta(xi0, xi01) = 1.
// This test verifies that this property is reproduced numerically within
// the accuracy expected from the chosen quadrature resolution, for all
// cells that it is expected to hold (e.g. not in the negative-trapped region).
bool KnockOn::CheckDeltaConservationProperty() {
    // Note: with higher nxi, relative error decreases
    real_t successRelErrorThreshold = 0.25;

    len_t nr = 5;
    len_t np = 4;
    len_t nxi = 60;

    len_t ntheta_interp = 50;
    len_t nrProfiles = 8;

    real_t pMin = 0;
    real_t pMax = 2;

    DREAM::FVM::Grid *grid =
        InitializeGridGeneralRPXi(nr, np, nxi, ntheta_interp, nrProfiles, pMin, pMax);

    bool success = _checkDeltaConservationProperty(grid, grid, successRelErrorThreshold);
    delete grid;
    return success;
}

// Like CheckDeltaConservationProperty, but when secondary and primary
// electrons live on different grids.
bool KnockOn::CheckDeltaConservationPropertyDifferentGrids() {
    // Note: with higher nxi, relative error decreases
    real_t successRelErrorThreshold = 0.25;

    len_t nr = 5;
    len_t np = 4;
    len_t nxi = 60;
    len_t nxi1 = 60;

    len_t ntheta_interp = 50;
    len_t nrProfiles = 8;

    real_t pMin = 0;
    real_t pMax = 2;

    DREAM::FVM::Grid *grid_knockon =
        InitializeGridGeneralRPXi(nr, np, nxi, ntheta_interp, nrProfiles, pMin, pMax);
    DREAM::FVM::Grid *grid_primary =
        InitializeGridGeneralRPXi(nr, np, nxi1, ntheta_interp, nrProfiles, pMin, pMax);

    bool success =
        _checkDeltaConservationProperty(grid_knockon, grid_primary, successRelErrorThreshold);
    delete grid_knockon;
    delete grid_primary;
    return success;
}

// In cylindrical geometry (constant B), the orbit-averaged delta function
// reduces to the analytically integrated gyro-averaged delta divided by
// the FVM xi-cell width. This test verifies that the general implementation
// reproduces this limit.
bool KnockOn::CheckCylindricalDeltaCalculation() {
    real_t successRelErrorThreshold = 1e-10;

    len_t nr = 5;
    len_t np = 4;
    len_t nxi = 60;

    real_t B0 = 1;

    real_t pMin = 0;
    real_t pMax = 2;

    DREAM::FVM::Grid *grid = InitializeGridRCylPXi(nr, np, nxi, B0, pMin, pMax);

    constexpr real_t eps = 5 * std::numeric_limits<real_t>::epsilon();
    bool success = true;
    real_t xi_stars[4] = {0.1, 0.5, 0.9, 1.0};

    for (len_t ir = 0; ir < nr; ir++) {
        DREAM::FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir);
        len_t Nxi = mg->GetNp2();
        for (len_t n = 0; n < 4; n++) {
            real_t xi_star = xi_stars[n];
            for (len_t l = 0; l < Nxi; l++) {
                real_t xi01 = mg->GetP2(l);
                for (len_t j = 0; j < Nxi; j++) {
                    real_t xi0_f1 = mg->GetP2_f(j);
                    real_t xi0_f2 = mg->GetP2_f(j + 1);
                    real_t z1;
                    real_t z2;
                    DREAM::KnockOnUtilities::AnalyticDelta::ComputeXi1Bounds(
                        z1, z2, xi0_f1, xi0_f2, xi_star
                    );

                    real_t dxi = mg->GetDp2(j);
                    real_t delta_expected = 0;
                    if (z1 < xi01 + eps && xi01 < z2 - eps) {
                        delta_expected = DREAM::KnockOnUtilities::AnalyticDelta::EvaluateXiIntegral(
                                             xi0_f1, xi0_f2, xi01, xi_star
                                         ) /
                                         dxi;
                    }
                    real_t delta = DREAM::KnockOnUtilities::EvaluateDeltaMatrixElementOnGrid(
                        ir, xi_star, j, l, grid, grid,
                        DREAM::KnockOnUtilities::N_POINTS_INTEGRAL_DEFAULT,
                        DREAM::KnockOnUtilities::MIDPOINT_RULE
                    );
                    if (fabs(delta - delta_expected) > successRelErrorThreshold) {
                        this->PrintError("failed at ir=%ld, n=%ld, j=%ld, l=%ld:\n", ir, n, j, l);
                        this->PrintError("  delta: %.8g\n", delta);
                        this->PrintError("  delta_expected: %.8g\n", delta_expected);
                        this->PrintError("  xi_star: %.8g\n", xi_star);
                        this->PrintError("  xi0 in [%.8g, %.8g]\n", xi0_f1, xi0_f2);
                        this->PrintError("  xi01: %.8g\n", xi01);
                        success = false;
                    }
                }
            }
        }
    }
    delete grid;
    return success;
}

// Compares the new implementation of the collision kernel with the previous
// implementation tailored to the Rosenbluth-Putvinski limit. This regression
// test ensures sufficient agreement between the two approaches in the regime
// where both are applicable.
// Uses the adaptive quadrature to demonstrate agreement to very high precision.
bool KnockOn::CheckAgreementWithOldRPTerm() {
    real_t successRelErrorThreshold = 1e-5;
    bool success = true;

    len_t nr = 5;
    len_t np = 4;
    len_t nxi = 30;

    len_t ntheta_interp = 100;
    len_t nrProfiles = 8;

    real_t pMin = 0;
    real_t pMax = 2;

    DREAM::FVM::Grid *grid =
        InitializeGridGeneralRPXi(nr, np, nxi, ntheta_interp, nrProfiles, pMin, pMax);
    DREAM::FVM::RadialGrid *rg = grid->GetRadialGrid();
    DREAM::FVM::FluxSurfaceAverager *fsa = rg->GetFluxSurfaceAverager();

    len_t idx_p = 1;
    real_t xi01 = 1;
    for (len_t ir = 0; ir < nr; ir++) {
        len_t Nxi = grid->GetNp2(ir);
        DREAM::FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir);
        real_t p = mg->GetP1(idx_p);
        real_t p1 = std::numeric_limits<real_t>::infinity();
        real_t xi_star = DREAM::KnockOnUtilities::Kinematics::EvaluateXiStar(p, p1);
        real_t Vp1 = 2 * M_PI * grid->GetVpVol(ir) * rg->GetFSA_B(ir);  // normalized to p1^2
        for (len_t j = 0; j < Nxi; j++) {
            real_t xi0_f1 = mg->GetP2_f(j);
            real_t xi0_f2 = mg->GetP2_f(j + 1);
            real_t theta1, theta2;
            DREAM::KnockOnUtilities::EstimateBoundingTheta(
                ir, j, Nxi - 1, theta1, theta2, grid, grid
            );
            real_t Vp = grid->GetVp(ir, idx_p, j);
            real_t VpVol = grid->GetVpVol(ir);
            real_t FSA_B = rg->GetFSA_B(ir);
            real_t old_delta =
                FSA_B * Vp / Vp1 * fsa->EvaluateAvalancheDeltaHat(ir, p, xi0_f1, xi0_f2, Vp, VpVol);
            real_t new_delta = DREAM::KnockOnUtilities::EvaluateOrbitAveragedDelta(
                ir, xi_star, xi01, xi0_f1, xi0_f2, Vp1, theta1, theta2, rg, 1000,
                DREAM::KnockOnUtilities::ADAPTIVE_TRAPEZOID
            );
            if (fabs(new_delta - old_delta) > old_delta * successRelErrorThreshold) {
                this->PrintError("FSA_B * Vp / Vp1: %.8g\n", FSA_B * Vp / Vp1);
                this->PrintError("Trapped boundary xi0T: %.8g\n", rg->GetXi0TrappedBoundary(ir));
                this->PrintError("non-zero contribution at xi in [%.3g, %.3g]\n", xi0_f1, xi0_f2);
                this->PrintError("old delta: %.8g\n", old_delta);
                this->PrintError("new delta: %.8g\n", new_delta);
                if (new_delta != 0) {
                    this->PrintError("ratio: %.8g\n", old_delta / new_delta);
                }
                success = false;
            }
        }
    }
    delete grid;
    return success;
}

// verify that integrating the differential cross section yields the "flux"
bool KnockOn::CheckMollerFluxIntegration() {
    real_t rtol = 1e-4;

    real_t p1 = 7;

    real_t p0_lower = 0.2;
    real_t p0_upper = 0.9;

    real_t mollerFluxAnalytic =
        (DREAM::KnockOnUtilities::Kinematics::EvaluateMollerFlux(p0_upper, p1) -
         DREAM::KnockOnUtilities::Kinematics::EvaluateMollerFlux(p0_lower, p1));

    // Manually trapz the differential cross section
    len_t N = 1000;
    real_t dp = (p0_upper - p0_lower) / (N - 1);
    real_t v = p0_lower / sqrt(1 + p0_lower * p0_lower);
    real_t f1 =
        v *
        DREAM::KnockOnUtilities::Kinematics::EvaluateMollerDifferentialCrossSection(p0_lower, p1);
    real_t mollerFluxNumeric = 0;
    for (len_t n = 1; n < N; n++) {
        real_t p0 = p0_lower + dp * n;
        v = p0 / sqrt(1 + p0 * p0);
        real_t f2 =
            v * DREAM::KnockOnUtilities::Kinematics::EvaluateMollerDifferentialCrossSection(p0, p1);
        mollerFluxNumeric += dp * (f2 + f1) / 2;
        f1 = f2;
    }

    return fabs(mollerFluxAnalytic - mollerFluxNumeric) <
           rtol * (fabs(mollerFluxAnalytic) + fabs(mollerFluxNumeric));
}

bool KnockOn::CheckMollerDifferentialConvergesToInfiniteLimit() {
    real_t rtol = 1e-6;
    real_t p = 0.5;
    real_t p1 = 1e8;
    real_t pinf = std::numeric_limits<real_t>::infinity();

    real_t sigma_large =
        DREAM::KnockOnUtilities::Kinematics::EvaluateMollerDifferentialCrossSection(p, p1);
    real_t sigma_inf =
        DREAM::KnockOnUtilities::Kinematics::EvaluateMollerDifferentialCrossSection(p, pinf);

    return fabs(sigma_large - sigma_inf) < rtol * (fabs(sigma_large) + fabs(sigma_inf));
}

bool KnockOn::CheckMollerFluxConvergesToInfiniteLimit() {
    real_t rtol = 1e-6;
    real_t p = 0.5;
    real_t p1 = 1e8;
    real_t pinf = std::numeric_limits<real_t>::infinity();

    real_t S_large = DREAM::KnockOnUtilities::Kinematics::EvaluateMollerFlux(p, p1);
    real_t S_inf = DREAM::KnockOnUtilities::Kinematics::EvaluateMollerFlux(p, pinf);

    return fabs(S_large - S_inf) < rtol * (fabs(S_large) + fabs(S_inf));
}

bool KnockOn::CheckXiStarConvergesToInfiniteLimit() {
    real_t rtol = 1e-6;
    real_t p = 0.5;
    real_t p1 = 1e8;
    real_t pinf = std::numeric_limits<real_t>::infinity();

    real_t xiStar_large = DREAM::KnockOnUtilities::Kinematics::EvaluateXiStar(p, p1);
    real_t xiStar_inf = DREAM::KnockOnUtilities::Kinematics::EvaluateXiStar(p, pinf);

    return fabs(xiStar_inf - xiStar_large) < rtol * (fabs(xiStar_large) + fabs(xiStar_inf));
}

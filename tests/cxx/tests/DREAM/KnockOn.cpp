/**
 * Implements unit tests of the KnockOnUtilities module.
 *
 * Where numerical accuracy is tested, we intentionally allow relatively loose
 * tolerances. This is partly to limit computational cost, but mainly because
 * this level of accuracy is sufficient for production use.
 *
 * In typical simulations, relatively low numerical precision is acceptable
 * since exact particle-number conservation is enforced by construction in
 * the knock-on operator.
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

    if (CheckDeltaConservationProperty())
        this->PrintOK("The delta function integrated over xi0 is sufficiently close to 1.");
    else {
        success = false;
        this->PrintError("Knock-on test failed: the delta-function does not integrate to 1.");
    }

    if (CheckAgreementWithOldRPTerm())
        this->PrintOK(
            "The new delta-function implementation agrees with the previous one for xi01=1"
        );
    else {
        success = false;
        this->PrintError("Knock-on test failed: the new delta-function implementation does not "
                         "agree with the old one.");
    }

    if (CheckCylindricalDeltaCalculation())
        this->PrintOK("The delta matrix elements agree exactly with the analytic result.");
    else {
        success = false;
        this->PrintError("Knock-on test failed: the delta matrix elements do not "
                         "agree exactly with theory in the cylindrical limit.");
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
                if (grid->IsNegativePitchTrappedIgnorableCell(ir, j))
                    continue;
                real_t xi0_f1 = mg->GetP2_f(j);
                real_t xi0_f2 = mg->GetP2_f(j + 1);
                for (len_t idx_l = 0; idx_l < 6; idx_l++) {
                    len_t l = xi01_indices[idx_l];
                    real_t xi01 = mg->GetP2(l);
                    if (grid->IsNegativePitchTrappedIgnorableCell(ir, l))
                        continue;
                    real_t Vp1 = grid->GetVpOverP2AtZero(ir)[l];
                    real_t theta1, theta2;
                    DREAM::KnockOnUtilities::EstimateBoundingTheta(ir, j, l, theta1, theta2, grid);

                    real_t delta1 = DREAM::KnockOnUtilities::EvaluateOrbitAveragedDelta(
                        ir, xi_star, xi01, xi0_f1, xi0_f2, Vp1, theta1, theta2, 30, rg
                    );
                    real_t delta2 = DREAM::KnockOnUtilities::EvaluateOrbitAveragedDelta(
                        ir, -xi_star, -xi01, xi0_f1, xi0_f2, Vp1, theta1, theta2, 30, rg
                    );
                    if (fabs(delta1 - delta2) > fabs(delta1) * successRelErrorThreshold) {
                        this->PrintError(
                            "Mirror symmetry failed at ir=%ld, j=%ld, l=%ld, xi_star=%.3g\n", ir, j,
                            l, xi_star
                        );
                        success = false;
                    }
                    real_t delta3 = DREAM::KnockOnUtilities::EvaluateOrbitAveragedDelta(
                        ir, -xi_star, xi01, -xi0_f2, -xi0_f1, Vp1, theta1, theta2, 30, rg
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
                if (grid->IsNegativePitchTrappedIgnorableCell(ir, j))
                    continue;
                // len_t j = xi0_indices[idx_j];
                for (len_t idx_l = 0; idx_l < 6; idx_l++) {
                    len_t l = xi01_indices[idx_l];
                    real_t delta_default =
                        DREAM::KnockOnUtilities::EvaluateDeltaMatrixElementOnGrid(
                            ir, xi_star, j, l, grid
                        );
                    real_t delta_hires = DREAM::KnockOnUtilities::EvaluateDeltaMatrixElementOnGrid(
                        ir, xi_star, j, l, grid, 300
                    );
                    if (fabs(delta_default - delta_hires) >
                        successRelErrorThreshold * (1 + fabs(delta_default))) {
                        this->PrintError("failed (j=%ld, l=%ld)\n", j, l);
                        this->PrintError("delta_default = %.4g\n", delta_default);
                        this->PrintError("delta_hires = %.4g\n", delta_hires);

                        success = false;
                    }
                }
            }
        }
    }
    delete grid;
    return success;
}

// The Delta function in our formulation of the knock-on operator satisfies
// an exact analytic normalization property:
//     \int dxi0 delta(xi0, xi01) = 1.
// This test verifies that this property is reproduced numerically within
// the accuracy expected from the chosen quadrature resolution.
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
    DREAM::FVM::RadialGrid *rg = grid->GetRadialGrid();

    bool success = true;
    real_t xi_stars[4] = {0.1, 0.5, 0.9, 1.0};

    for (len_t ir = 0; ir < nr; ir++) {
        real_t xi0T = rg->GetXi0TrappedBoundary(ir);
        DREAM::FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir);
        len_t Nxi = mg->GetNp2();
        for (len_t n = 0; n < 4; n++) {
            real_t xi_star = xi_stars[n];
            for (len_t l = 0; l < Nxi; l++) {
                real_t deltaXiIntegral = 0;
                for (len_t j = 0; j < Nxi; j++) {
                    real_t delta = DREAM::KnockOnUtilities::EvaluateDeltaMatrixElementOnGrid(
                        ir, xi_star, j, l, grid
                    );
                    real_t dxi = mg->GetDp2(j);
                    deltaXiIntegral += dxi * delta;
                }
                real_t xi01_f1 = mg->GetP2_f(l);
                real_t xi01_f2 = mg->GetP2_f(l + 1);
                real_t TPDist = 0.05;
                bool awayFromTPBoundaries =
                    (fabs(fabs(xi01_f1) - xi0T) > TPDist && fabs(fabs(xi01_f2) - xi0T) > TPDist &&
                     !(xi01_f1 < -xi0T && xi01_f2 > -xi0T) && !(xi01_f1 < xi0T && xi01_f2 > xi0T));

                if (!grid->IsNegativePitchTrappedIgnorableCell(ir, l) &&
                    (fabs(deltaXiIntegral - 1) > successRelErrorThreshold) &&
                    awayFromTPBoundaries) {
                    real_t xi01 = mg->GetP2(l);
                    this->PrintError("failed:\n");
                    this->PrintError("  trapped boundary: %.4g\n", xi0T);
                    this->PrintError("  xi01: %.4g (%.4g, %.4g)\n", xi01, xi01_f1, xi01_f2);
                    this->PrintError("  deltaXiIntegral: %.4g\n", deltaXiIntegral);
                    success = false;
                }
            }
        }
    }
    delete grid;
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
                    DREAM::KnockOnUtilities::ComputeXi1Bounds(z1, z2, xi0_f1, xi0_f2, xi_star);

                    real_t dxi = mg->GetDp2(j);
                    real_t delta_expected = 0;
                    if (z1 < xi01 + eps && xi01 < z2 - eps) {
                        delta_expected = DREAM::KnockOnUtilities::EvaluateDeltaIntervalContribution(
                                             xi0_f1, xi0_f2, xi01, xi_star
                                         ) /
                                         dxi;
                    }
                    real_t delta = DREAM::KnockOnUtilities::EvaluateDeltaMatrixElementOnGrid(
                        ir, xi_star, j, l, grid
                    );
                    if (fabs(delta - delta_expected) > successRelErrorThreshold) {
                        this->PrintError("failed at ir=%ld, j=%ld, l=%ld:\n", ir, j, l);
                        this->PrintError("  delta: %.4g\n", delta);
                        this->PrintError("  delta_expected: %.4g\n", delta_expected);
                        this->PrintError("  xi_star: %.4g\n", xi_star);
                        this->PrintError("  xi0 in [%.4g, %.4g]\n", xi0_f1, xi0_f2);
                        this->PrintError("  xi01: %.4g\n", xi01);
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
bool KnockOn::CheckAgreementWithOldRPTerm() {
    real_t successRelErrorThreshold = 0.1;
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
        real_t xi_star = DREAM::KnockOnUtilities::EvaluateXiStar(p, p1);
        real_t Vp1 = 2 * M_PI * grid->GetVpVol(ir) * rg->GetFSA_B(ir); // normalized to p1^2
        for (len_t j = 0; j < Nxi; j++) {
            real_t xi0_f1 = mg->GetP2_f(j);
            real_t xi0_f2 = mg->GetP2_f(j + 1);
            real_t theta1, theta2;
            DREAM::KnockOnUtilities::EstimateBoundingTheta(ir, j, Nxi - 1, theta1, theta2, grid);
            real_t Vp = grid->GetVp(ir, idx_p, j);
            real_t VpVol = grid->GetVpVol(ir);
            real_t FSA_B = rg->GetFSA_B(ir);
            real_t old_delta =
                FSA_B * Vp / Vp1 * fsa->EvaluateAvalancheDeltaHat(ir, p, xi0_f1, xi0_f2, Vp, VpVol);
            real_t new_delta = DREAM::KnockOnUtilities::EvaluateOrbitAveragedDelta(
                ir, xi_star, xi01, xi0_f1, xi0_f2, Vp1, theta1, theta2, 10000, rg
            );
            if (fabs(new_delta - old_delta) > old_delta * successRelErrorThreshold) {
                this->PrintError("FSA_B * Vp / Vp1: %.4g\n", FSA_B * Vp / Vp1);
                this->PrintError("Trapped boundary xi0T: %.4g\n", rg->GetXi0TrappedBoundary(ir));
                this->PrintError("non-zero contribution at xi in [%.3g, %.3g]\n", xi0_f1, xi0_f2);
                this->PrintError("old delta: %.4g\n", old_delta);
                this->PrintError("new delta: %.4g\n", new_delta);
                if (new_delta != 0){
                    this->PrintError("ratio: %.4g\n", old_delta / new_delta);
                }
                success = false;
            }
        }
    }
    delete grid;
    return success;
}

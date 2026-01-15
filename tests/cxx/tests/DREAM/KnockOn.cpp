#include "KnockOn.hpp"
#include "DREAM/Equations/KnockOnOperator.hpp"
#include "DREAM/Equations/KnockOnUtilities.hpp"
#include "FVM/Grid/Grid.hpp"
#include <chrono>
#include <cstdio>

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

    return success;
}

namespace {}

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
                    DREAM::KnockOnUtilities::estimateBoundingTheta(ir, j, l, theta1, theta2, grid);

                    real_t delta1 = DREAM::KnockOnUtilities::EvaluateDeltaContribution(
                        ir, xi_star, xi01, xi0_f1, xi0_f2, Vp1, theta1, theta2, 30, grid
                    );
                    real_t delta2 = DREAM::KnockOnUtilities::EvaluateDeltaContribution(
                        ir, -xi_star, -xi01, xi0_f1, xi0_f2, Vp1, theta1, theta2, 30, grid
                    );
                    if (fabs(delta1 - delta2) > fabs(delta1) * successRelErrorThreshold)
                        success = false;

                    real_t delta3 = DREAM::KnockOnUtilities::EvaluateDeltaContribution(
                        ir, -xi_star, xi01, -xi0_f2, -xi0_f1, Vp1, theta1, theta2, 30, grid
                    );
                    if (fabs(delta2 - delta3) > fabs(delta2) * successRelErrorThreshold)
                        success = false;
                }
            }
        }
    }
    // delete [] grid;
    return success;
}

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

    DREAM::KnockOnOperator *boltz_default = new DREAM::KnockOnOperator(grid, 0.1);
    DREAM::KnockOnOperator *boltz_hires = new DREAM::KnockOnOperator(grid, 0.1, 50, 300);

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
                    real_t delta_default = boltz_default->EvaluateDelta(ir, xi_star, j, l);
                    real_t delta_hires = boltz_hires->EvaluateDelta(ir, xi_star, j, l);
                    if (fabs(delta_default - delta_hires) >
                        successRelErrorThreshold * (1 + fabs(delta_default))) {
                        printf("failed (j=%ld, l=%ld)\n", j, l);
                        printf("delta_default = %.4g\n", delta_default);
                        printf("delta_hires = %.4g\n", delta_hires);

                        success = false;
                    }
                }
            }
        }
    }

    // delete [] boltz_default;
    // delete [] boltz_hires;
    // delete [] grid;

    return success;
}

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

    DREAM::KnockOnOperator *boltz = new DREAM::KnockOnOperator(grid, 0.1);

    bool success = true;
    real_t xi_stars[4] = {0.1, 0.5, 0.9, 1.0};
    auto start = std::chrono::high_resolution_clock::now();

    for (len_t ir = 0; ir < nr; ir++) {
        real_t xi0T = grid->GetRadialGrid()->GetXi0TrappedBoundary(ir);
        DREAM::FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir);
        len_t Nxi = mg->GetNp2();
        for (len_t n = 0; n < 4; n++) {
            real_t xi_star = xi_stars[n];
            for (len_t l = 0; l < Nxi; l++) {
                real_t deltaXiIntegral = 0;
                for (len_t j = 0; j < Nxi; j++) {
                    real_t delta = boltz->EvaluateDelta(ir, xi_star, j, l);
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
                    printf("failed:\n");
                    printf("  trapped boundary: %.4g\n", xi0T);
                    printf("  xi01: %.4g (%.4g, %.4g)\n", xi01, xi01_f1, xi01_f2);
                    printf("  deltaXiIntegral: %.4g\n", deltaXiIntegral);
                    success = false;
                }
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<real_t> elapsed = end - start;

    printf("Time spent: %.6f seconds\n", elapsed.count());
    return success;
}

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

    len_t idx_p = 1;
    real_t xi01 = 1;
    for (len_t ir = 0; ir < nr; ir++) {
        len_t Nxi = grid->GetNp2(ir);
        DREAM::FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir);
        real_t p = mg->GetP1(idx_p);
        real_t p1 = std::numeric_limits<real_t>::infinity();
        real_t xi_star = DREAM::KnockOnUtilities::evaluateXiStar(p, p1);
        real_t Vp1 = 2 * M_PI * grid->GetVpVol(ir) *
                     grid->GetRadialGrid()->GetFSA_B(ir); // normalized to p1^2
        for (len_t j = 0; j < Nxi; j++) {
            real_t xi0_f1 = mg->GetP2_f(j);
            real_t xi0_f2 = mg->GetP2_f(j + 1);
            real_t theta1, theta2;
            DREAM::KnockOnUtilities::estimateBoundingTheta(ir, j, Nxi - 1, theta1, theta2, grid);
            real_t Vp = grid->GetVp(ir, idx_p, j);
            real_t VpVol = grid->GetVpVol(ir);
            real_t FSA_B = grid->GetRadialGrid()->GetFSA_B(ir);
            real_t old_delta =
                FSA_B * Vp / Vp1 *
                grid->GetRadialGrid()->GetFluxSurfaceAverager()->EvaluateAvalancheDeltaHat(
                    ir, p, xi0_f1, xi0_f2, Vp, VpVol
                );
            real_t new_delta = DREAM::KnockOnUtilities::EvaluateDeltaContribution(
                ir, xi_star, xi01, xi0_f1, xi0_f2, Vp1, theta1, theta2, 10000, grid
            );
            if (fabs(new_delta - old_delta) > old_delta * successRelErrorThreshold) {
                printf("FSA_B * Vp / Vp1: %.4g\n", FSA_B * Vp / Vp1);
                printf(
                    "Trapped boundary xi0T: %.4g\n",
                    grid->GetRadialGrid()->GetXi0TrappedBoundary(ir)
                );
                printf("non-zero contribution at xi in [%.3g, %.3g]\n", xi0_f1, xi0_f2);
                printf("old delta: %.4g\n", old_delta);
                printf("new delta: %.4g\n", new_delta);
                if (new_delta != 0)
                    printf("ratio: %.4g\n", old_delta / new_delta);
                success = false;
            }
        }
    }
    return success;
}

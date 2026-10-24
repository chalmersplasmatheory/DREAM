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

#include <algorithm>
#include <memory>
#include <utility>
#include <vector>

#include "DREAM/Equations/KnockOnUtilities.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/AnalyticBRadialGridGenerator.hpp"
#include "FVM/Grid/PXiGrid/PXiMomentumGrid.hpp"
#include "FVM/Grid/PXiGrid/PUniformGridGenerator.hpp"
#include "FVM/Grid/PXiGrid/XiCustomGridGenerator.hpp"

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

    if (CheckPoloidalNumerators())
        this->PrintOK("Raw poloidal numerators agree with an independent dense orbit integral.");
    else {
        success = false;
        this->PrintError("Knock-on unnormalized poloidal numerator regression failed.");
    }

    if (CheckPitchCellAverages())
        this->PrintOK("Primary cell averages preserve raw volumes and boundary-row weights.");
    else {
        success=false;
        this->PrintError("Primary pitch cell-average regression failed.");
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

    if (CheckDeltaConservationPropertyConstantB())
        this->PrintOK(
            "The delta function integrated over xi0 is sufficiently close to 1 when primaries and "
            "secondaries live on the same grid in a constant magnetic field."
        );
    else {
        success = false;
        this->PrintError(
            "Knock-on test failed: the delta-function does not integrate to 1 when primaries and "
            "secondaries live on the same grid in a constant magnetic field."
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
        this->PrintOK("The Møller cross section correctly integrates to the Møller flux.");
    else {
        success = false;
        this->PrintError("The Møller cross section fails to integrate to the Møller flux.");
    }

    if (CheckMollerDifferentialConvergesToInfiniteLimit())
        this->PrintOK("The Møller differential cross-section converges to the infinite limit.");
    else {
        success = false;
        this->PrintError(
            "The Møller differential cross-section does not converge to the infinite limit."
        );
    }

    if (CheckMollerFluxConvergesToInfiniteLimit())
        this->PrintOK("The Møller integrated flux converges to the infinite limit.");
    else {
        success = false;
        this->PrintError("The Møller integrated flux does not converge to the infinite limit.");
    }

    if (CheckXiStarConvergesToInfiniteLimit())
        this->PrintOK("The scattering angle xiStar converges to the infinite limit.");
    else {
        success = false;
        this->PrintError("The scattering angle xiStar does not converge to the infinite limit.");
    }

    if (CheckMollerSConservationProperty())
        this->PrintOK(
            "The Møller momentum S matrix integrated over p reduces to the correct total cross "
            "section"
        );
    else {
        success = false;
        this->PrintError(
            "The Møller momentum S matrix integrated over p does not produce the correct total "
            "cross section"
        );
    }

    return success;
}

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
    real_t xi_stars[3] = {0.0, 0.5, 1.0};
    len_t xi0_indices[6] = {0, 5, 16, 20, 25, 29};
    len_t xi01_indices[6] = {0, 5, 16, 20, 25, 29};
    for (len_t ir = 0; ir < nr; ir++) {
        DREAM::FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir);
        for (len_t n = 0; n < 3; n++) {
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

                    real_t delta1 = DREAM::KnockOnUtilities::EvaluateOrbitAveragedDelta(
                        ir, xi_star, xi01, xi0_f1, xi0_f2, Vp1, rg, 8
                    );
                    real_t delta2 = DREAM::KnockOnUtilities::EvaluateOrbitAveragedDelta(
                        ir, -xi_star, -xi01, xi0_f1, xi0_f2, Vp1, rg, 8
                    );
                    if (fabs(delta1 - delta2) > fabs(delta1) * successRelErrorThreshold) {
                        this->PrintError(
                            "Mirror symmetry failed at ir=%ld, j=%ld, l=%ld, xi_star=%.3g", ir, j,
                            l, xi_star
                        );
                        success = false;
                    }
                    real_t delta3 = DREAM::KnockOnUtilities::EvaluateOrbitAveragedDelta(
                        ir, -xi_star, xi01, -xi0_f2, -xi0_f1, Vp1, rg, 8
                    );
                    if (fabs(delta2 - delta3) > fabs(delta2) * successRelErrorThreshold) {
                        this->PrintError(
                            "Mirror symmetry failed at ir=%ld, j=%ld, l=%ld, xi_star=%.3g", ir, j,
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

    const real_t poloidalOrder = 8;

    bool success = true;
    real_t xi_stars[3] = {0.0, 0.5, 1.0};
    // len_t xi0_indices[6] = {0, 5, 16, 20, 25, 29};
    len_t xi01_indices[6] = {0, 5, 16, 20, 25, 29};
    // real_t xi01s[6] = {-1, -0.9, 0.02, 0.2, 0.8, 1.0};
    for (len_t ir = 0; ir < nr; ir++) {
        for (len_t n = 0; n < 3; n++) {
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
                            ir, xi_star, j, l, grid, grid, poloidalOrder
                        );
                    real_t delta_hires = DREAM::KnockOnUtilities::EvaluateDeltaMatrixElementOnGrid(
                        ir, xi_star, j, l, grid, grid, 32
                    );
                    if (fabs(delta_default - delta_hires) >
                        successRelErrorThreshold * (1 + fabs(delta_default))) {
                        this->PrintError("failed (j=%ld, l=%ld)", j, l);
                        this->PrintError("delta_default = %.8g", delta_default);
                        this->PrintError("delta_hires = %.8g", delta_hires);

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
    real_t xi_stars[3] = {0.0, 0.5, 1.0};
    len_t n_theta = 10;
    constexpr real_t theta_roundoff_safety_factor = 1e-8;

    for (len_t ir = 0; ir < rg->GetNr(); ir++) {
        const DREAM::FVM::MomentumGrid *mg = grid_primary->GetMomentumGrid(ir);

        len_t Nxi = grid_knockon->GetNp2(ir);
        len_t Nxi1 = grid_primary->GetNp2(ir);
        for (len_t n = 0; n < 3; n++) {
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
                        this->PrintError("failed at ir=%ld, n=%ld, l=%ld, m=%ld:", ir, n, l, m);
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
    const auto *rg = grid_knockon->GetRadialGrid();
    for (len_t ir=0; ir<rg->GetNr(); ++ir) {
        const auto *mg=grid_knockon->GetMomentumGrid(ir);
        std::vector<real_t> column(mg->GetNp2());
        for (real_t chi : {0.0,0.5,1.0})
            for (len_t l=0; l<grid_primary->GetNp2(ir); ++l) {
                const bool isVoid=grid_primary->GetVpOverP2AtZero(ir)[l]==0;
                const real_t rawNorm=DREAM::KnockOnUtilities::SetDeltaMatrixColumnOnGrid(
                    ir,chi,l,grid_knockon,grid_primary,column.data(),8,4
                );
                real_t mass=0;
                for (len_t j=0; j<column.size(); ++j) {
                    if (!std::isfinite(column[j]) || column[j]<0) return false;
                    if (grid_knockon->GetVpOverP2AtZero(ir)[j]==0 && column[j]!=0) return false;
                    mass+=mg->GetDp2(j)*column[j];
                }
                // Hquad/Hgrid is NOT generally one: the grid stores its own
                // bounce-volume approximation. Check the normalized invariant
                // here, including boundary/origin cells and native voids.
                if (isVoid ? (mass!=0 || rawNorm!=0)
                    : (!(rawNorm>0) || !std::isfinite(rawNorm)
                        || std::abs(mass-1)>successRelErrorThreshold)) return false;
            }
    }
    return true;
}

// The Delta function in our formulation of the knock-on operator satisfies
// an exact analytic normalization property:
//     \int dxi0 delta(xi0, xi01) = 1.
// This test verifies that this property is reproduced numerically within
// the accuracy expected from the chosen quadrature resolution, for all
// cells that it is expected to hold (e.g. not in the negative-trapped region).
bool KnockOn::CheckDeltaConservationProperty() {
    // Conservation is enforced on each completed non-void column.
    real_t successRelErrorThreshold = 1e-12;

    len_t nr = 5;
    len_t np = 4;
    len_t nxi = 61;

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

// The Delta function in our formulation of the knock-on operator satisfies
// an exact analytic normalization property:
//     \int dxi0 delta(xi0, xi01) = 1.
// This test verifies that this property is reproduced numerically within
// the accuracy expected from the chosen quadrature resolution, for all
// cells that it is expected to hold (e.g. not in the negative-trapped region).
bool KnockOn::CheckDeltaConservationPropertyConstantB() {
    // Conservation is enforced on each completed non-void column.
    real_t successRelErrorThreshold = 1e-12;

    len_t nr = 5;
    len_t np = 4;
    len_t nxi = 60;

    real_t pMin = 0;
    real_t pMax = 2;

    real_t B0 = 1.0;

    DREAM::FVM::Grid *grid = InitializeGridRCylPXi(nr, np, nxi, B0, pMin, pMax);

    bool success = _checkDeltaConservationProperty(grid, grid, successRelErrorThreshold);
    delete grid;
    return success;
}

// Like CheckDeltaConservationProperty, but when secondary and primary
// electrons live on different grids.
bool KnockOn::CheckDeltaConservationPropertyDifferentGrids() {
    // Conservation is enforced on each completed non-void column.
    real_t successRelErrorThreshold = 1e-12;

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
    namespace K=DREAM::KnockOnUtilities;
    // For chi=1 in a cylinder the exact cell integral is the interval overlap.
    // Unlike, odd/even grids test the primary cell average, not a centre sample.
    for (const auto sizes : {std::pair<len_t,len_t>{13,18},{18,13}}) {
        std::unique_ptr<DREAM::FVM::Grid> target(InitializeGridRCylPXi(1,2,sizes.first,1,0,2));
        std::unique_ptr<DREAM::FVM::Grid> primary(InitializeGridRCylPXi(1,2,sizes.second,1,0,2));
        const auto *s=target->GetMomentumGrid(0), *p=primary->GetMomentumGrid(0);
        for (len_t order : {len_t(1),len_t(4)})
            for (len_t l=0; l<p->GetNp2(); ++l)
                for (len_t j=0; j<s->GetNp2(); ++j) {
                    const real_t overlap=std::max(real_t(0),std::min(s->GetP2_f(j+1),p->GetP2_f(l+1))
                        -std::max(s->GetP2_f(j),p->GetP2_f(l)));
                    const real_t expected=overlap/(p->GetDp2(l)*s->GetDp2(j));
                    const real_t value=K::EvaluateDeltaMatrixElementOnGrid(
                        0,1,j,l,target.get(),primary.get(),8,order
                    );
                    if (!std::isfinite(value) || std::abs(value-expected)>1e-11) return false;
                }
    }
    // Regression: contact rounding in a coarse toroidal table used to pass an
    // asin argument just outside [-1,1]. The saturated probability is finite.
    for (real_t x : {-1.0,1.0}) {
        const real_t value=K::AnalyticDelta::EvaluateXiAntiderivative(
            std::nextafter(x,2*x),0,0
        );
        if (value!=0.5*x) return false;
    }
    return true;
}

// Compares the new implementation of the collision kernel with the previous
// implementation tailored to the Rosenbluth-Putvinski limit. This regression
// test ensures sufficient agreement between the two approaches in the regime
// where both are applicable.
// Uses a higher mapped Gauss order to compare with the independent RP implementation.
bool KnockOn::CheckAgreementWithOldRPTerm() {
    real_t successRelErrorThreshold = 1e-3;
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
            real_t Vp = grid->GetVp(ir, idx_p, j);
            real_t VpVol = grid->GetVpVol(ir);
            real_t FSA_B = rg->GetFSA_B(ir);
            real_t old_delta =
                FSA_B * Vp / Vp1 * fsa->EvaluateAvalancheDeltaHat(ir, p, xi0_f1, xi0_f2, Vp, VpVol);
            real_t new_delta = DREAM::KnockOnUtilities::EvaluateOrbitAveragedDelta(
                ir, xi_star, xi01, xi0_f1, xi0_f2, Vp1, rg, 16
            );
            if (fabs(new_delta - old_delta) > old_delta * successRelErrorThreshold) {
                this->PrintError("FSA_B * Vp / Vp1: %.8g", FSA_B * Vp / Vp1);
                this->PrintError("Trapped boundary xi0T: %.8g", rg->GetXi0TrappedBoundary(ir));
                this->PrintError("non-zero contribution at xi in [%.3g, %.3g]", xi0_f1, xi0_f2);
                this->PrintError("old delta: %.8g", old_delta);
                this->PrintError("new delta: %.8g", new_delta);
                if (new_delta != 0) {
                    this->PrintError("ratio: %.8g", old_delta / new_delta);
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

// Test the conservation property that if the Møller flux matrix between the two grids
// is integrated over knock-on momenta, yields the expected total cross section.
bool KnockOn::_checkMollerSConservationProperty(
    const DREAM::FVM::Grid *grid_knockon, const DREAM::FVM::Grid *grid_primary
) {
    real_t rtol = 100 * std::numeric_limits<real_t>::epsilon();
    real_t atol = 100 * std::numeric_limits<real_t>::epsilon();
    bool success = true;
    for (len_t ir = 0; ir < grid_knockon->GetNr(); ir++) {
        auto *mgK = grid_knockon->GetMomentumGrid(ir);
        auto *mgP = grid_primary->GetMomentumGrid(ir);

        real_t pMaxK =
            DREAM::KnockOnUtilities::Kinematics::MaximumKnockOnMomentum(mgK->GetP1_f(mgK->GetNp1())
            );
        // test for various cutoffs that are challenging in different ways
        std::vector<real_t> pCutoffs;
        pCutoffs.push_back(mgK->GetP1_f(2));
        pCutoffs.push_back(mgK->GetP1(2));
        pCutoffs.push_back(pMaxK);                    // on grid boundary
        pCutoffs.push_back(pMaxK - 1e-10);            // just inside the grid
        pCutoffs.push_back(mgK->GetP1_f(2) - 1e-10);  // near cell edges
        pCutoffs.push_back(mgK->GetP1_f(2) + 1e-10);

        for (len_t k = 0; k < mgP->GetNp1(); k++) {
            for (real_t pCutoff : pCutoffs) {
                real_t expected =
                    DREAM::KnockOnUtilities::EvaluateMollerFluxIntegratedOverKnockonGrid(
                        k, grid_knockon, grid_primary, pCutoff
                    );
                real_t sum = 0;
                for (len_t i = 0; i < mgK->GetNp1(); i++) {
                    real_t dp = mgK->GetDp1(i);
                    sum += dp * DREAM::KnockOnUtilities::EvaluateMollerFluxMatrixElementOnGrid(
                                    i, k, grid_knockon, grid_primary, pCutoff
                                );
                }
                if (fabs(sum - expected) > atol + rtol * (fabs(sum) + fabs(expected))) {
                    this->PrintError("Failed at ir=%ld, k=%ld, pCutoff=%.8g", ir, k, pCutoff);
                    this->PrintError("    sum: %.8g", sum);
                    this->PrintError("    expected: %.8g", expected);
                    success = false;
                }
            }
        }
    }
    return success;
}

// Test conservation property of the "Møller S matrix"
// when knock-on and runaway grids are the same.
bool KnockOn::CheckMollerSConservationProperty() {
    len_t nr = 2;
    len_t np = 50;
    len_t nxi = 2;

    len_t ntheta_interp = 50;
    len_t nrProfiles = 8;

    real_t pMin = 0;
    real_t pMax = 4;

    DREAM::FVM::Grid *grid =
        InitializeGridGeneralRPXi(nr, np, nxi, ntheta_interp, nrProfiles, pMin, pMax);

    bool success = _checkMollerSConservationProperty(grid, grid);
    delete grid;
    return success;
}

/**
 * Reference uses a dense midpoint rule on one globally mapped orbit and the
 * elementary gyro-cone probability. It uses neither contact cuts nor the
 * production local-delta/metric helpers. Check raw values before normalization.
 */
bool KnockOn::CheckPoloidalNumerators() {
    namespace K = DREAM::KnockOnUtilities;
    namespace Q = K::Quadrature;
    // Pole sources on faces, including +/-1, must belong to exactly one
    // outgoing cell in cylindrical geometry.
    for (real_t u : {-1.0,1.0})
        for (real_t chi : {-1.0,-0.5,0.0,0.5,1.0}) {
            real_t sum=0;
            for (len_t j=0; j<4; ++j) {
                real_t ratio;
                sum+=0.5*K::AnalyticDelta::EvaluateLocalContribution(
                    -1+0.5*j,-0.5+0.5*j,u,1,chi,0,ratio
                );
            }
            if (sum!=1) return false;
        }
    std::unique_ptr<DREAM::FVM::Grid> grid(InitializeGridGeneralRPXi(2,2,12,100,8));
    auto *rg=grid->GetRadialGrid();
    auto *fsa=rg->GetFluxSurfaceAverager();
    const Q::GaussLegendreRule rule(8);
    for (len_t ir=0; ir<grid->GetNr(); ++ir) {
        Q::OrbitGeometry geometry(fsa,ir);
        const real_t T=std::sqrt(geometry.GetBetaMax());
        for (real_t u : {0.6*T,T*(1-1e-4),T*(1+1e-4),(1+T)/2,1.0}) {
            const auto bounds=geometry.GetOrbitBounds(u);
            auto reference=[&](real_t lo, real_t hi, real_t chi, len_t n) {
                real_t result=0;
                const real_t Bmin=fsa->GetBmin(ir,DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION);
                for (len_t k=0; k<n; ++k) {
                    const real_t t=(k+0.5)/n;
                    const real_t theta=bounds.lower+(bounds.upper-bounds.lower)*std::pow(std::sin(M_PI*t/2),2);
                    const real_t weight=(bounds.upper-bounds.lower)*M_PI/2*std::sin(M_PI*t)/n;
                    real_t B,J,R,nabla;
                    fsa->GeometricQuantitiesAtTheta(ir,theta,B,J,R,nabla,DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION);
                    const real_t b=std::max(1.0,B/Bmin), beta=1-1/b;
                    if (beta>=u*u) continue;
                    const real_t eta=std::sqrt((u*u-beta)/(1-beta));
                    auto local=[&](real_t v) { return std::copysign(std::sqrt(std::max(0.0,(v*v-beta)/(1-beta))),v); };
                    const real_t a=chi*eta, d=std::sqrt(std::max(0.0,(1-eta*eta)*(1-chi*chi)));
                    real_t probability;
                    if (d==0) probability=local(lo)<a && a<=local(hi) ? 1 : 0;
                    else {
                        auto angle=[](real_t x) { return std::acos(std::max(-1.0,std::min(1.0,x))); };
                        probability=(angle((local(lo)-a)/d)-angle((local(hi)-a)/d))/M_PI;
                    }
                    result+=weight*4*M_PI*M_PI*J*b*u/eta*probability/(hi-lo);
                }
                return result;
            };
            const real_t volume=2*reference(-1,1,0.55,32768);
            for (real_t chi : {0.1,0.55,0.9999})
                for (const auto cell : {Q::Interval{-1,1},Q::Interval{0.15,0.85}}) {
                    const real_t ref=reference(cell.lower,cell.upper,chi,32768);
                    const real_t check=reference(cell.lower,cell.upper,chi,16384);
                    const real_t value=K::EvaluateOrbitNumerator(ir,chi,u,cell.lower,cell.upper,rg,geometry,rule);
                    const real_t scale=(cell.upper-cell.lower)/volume;
                    if (!std::isfinite(value) || value<0 || std::abs(ref-check)*scale>2e-4
                        || std::abs(value-ref)*scale>1e-3) {
                        PrintError("Raw orbit: ir=%ld u=%.16g chi=%.8g value=%.8g reference=%.8g",(long)ir,u,chi,value,ref);
                        return false;
                    }
                }
        }
    }
    return true;
}

namespace {
// Grid deliberately does not own its radial grid. These fixtures own both.
struct PitchGridDeleter {
    void operator()(DREAM::FVM::Grid *grid) const {
        auto *rg=grid->GetRadialGrid();
        delete grid;
        delete rg;
    }
};
using OwnedPitchGrid=std::unique_ptr<DREAM::FVM::Grid,PitchGridDeleter>;

// A circular surface at r/R0=epsilon with independently chosen pitch faces.
// The generator owns the profile descriptor; its constant arrays live statically.
DREAM::FVM::Grid *PitchTestGrid(real_t epsilon, const std::vector<real_t>& faces) {
    using namespace DREAM::FVM;
    static const real_t r[]={0,1}, one[]={1,1}, zero[]={0,0};
    auto *profiles=new AnalyticBRadialGridGenerator::shape_profiles{
        2,2,2,2,2, one,r,zero,r,one,r,zero,r,zero,r
    };
    auto *generator=new AnalyticBRadialGridGenerator(1,epsilon-0.01,epsilon+0.01,1,128,profiles);
    auto *rg=new RadialGrid(generator);
    auto *mg=new PXiGrid::PXiMomentumGrid(new PXiGrid::MomentumGridGenerator(
        new PXiGrid::PUniformGridGenerator(2,0,2),
        new PXiGrid::XiCustomGridGenerator(faces.data(),faces.size()-1)
    ),0,rg);
    auto *grid=new Grid(rg,mg);
    grid->Rebuild(0);
    return grid;
}
}

bool KnockOn::CheckPitchCellAverages() {
    namespace K=DREAM::KnockOnUtilities;
    namespace Q=K::Quadrature;
    // No sign symmetry or shared faces: the negative boundary row and its
    // receiving positive row have different widths; one cell crosses zero.
    const std::vector<real_t> sf={-1,-.82,-.55,-.2,.12,.43,.72,.86,1};
    const std::vector<real_t> pf={-1,-.79,-.51,-.11,.17,.57,.83,1};
    for (real_t epsilon : {0.05,0.3}) {
        OwnedPitchGrid target(PitchTestGrid(epsilon,sf));
        OwnedPitchGrid primary(PitchTestGrid(epsilon,pf));
        auto *rg=primary->GetRadialGrid();
        auto *fsa=rg->GetFluxSurfaceAverager();
        Q::OrbitGeometry geometry(fsa,0);
        const real_t T=std::sqrt(geometry.GetBetaMax());
        const real_t Bmin=fsa->GetBmin(0,DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION);
        // Independent continuum cell volume: integrate in theta AFTER changing
        // u to local pitch eta. The singular u/eta metric cancels du exactly.
        // No production contacts, orbit integrals, or stored Vp enter this reference.
        auto volume=[&](real_t lo,real_t hi,len_t n=32768) {
            if (hi<=lo) return real_t(0);
            real_t result=0;
            for (len_t k=0; k<n; ++k) {
                const real_t theta=2*M_PI*(k+0.5)/n;
                real_t B,J,R,nabla;
                fsa->GeometricQuantitiesAtTheta(0,theta,B,J,R,nabla,DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION);
                const real_t beta=std::max(real_t(0),1-Bmin/B);
                auto eta=[&](real_t u) {
                    return std::copysign(std::sqrt(std::max(real_t(0),(u*u-beta)/(1-beta))),u);
                };
                real_t length=0;
                // Negative passing, positive trapped (both orbit signs), positive passing.
                for (const auto part : {Q::Interval{-1,-T},{0,T},{T,1}}) {
                    const real_t a=std::max(lo,part.lower), b=std::min(hi,part.upper);
                    if (b>a) length+=(eta(b)-eta(a))*(part.lower==0 ? 2 : 1);
                }
                result+=4*M_PI*M_PI*J*length*(2*M_PI/n);
            }
            return result;
        };
        std::vector<real_t> column(sf.size()-1);
        for (len_t l=0; l+1<pf.size(); ++l) {
            if (primary->GetVpOverP2AtZero(0)[l]==0) continue;
            const real_t H=volume(pf[l],pf[l+1]);
            if (!(H>0) || std::abs(volume(pf[l],pf[l+1],16384)/H-1)>1e-4) return false;
            const real_t Hgrid=(pf[l+1]-pf[l])*primary->GetVpOverP2AtZero(0)[l];
            for (real_t chi : {0.0,0.55,1.0}) {
                const real_t norm=K::SetDeltaMatrixColumnOnGrid(
                    0,chi,l,target.get(),primary.get(),column.data(),8,4
                );
                if (!std::isfinite(norm) || std::abs(norm*Hgrid/H-1)>0.005) {
                    PrintError("Raw pitch volume: eps=%g l=%ld chi=%g Hquad/H=%g",epsilon,(long)l,chi,norm*Hgrid/H);
                    return false;
                }
                if (chi!=1) continue;
                for (len_t j=0; j<column.size(); ++j) {
                    real_t expected=0;
                    if (target->GetVpOverP2AtZero(0)[j]!=0) {
                        expected=volume(std::max(pf[l],sf[j]),std::min(pf[l+1],sf[j+1]));
                        if (sf[j]<=T && T<sf[j+1])
                            for (len_t k=0; k+1<sf.size(); ++k)
                                if (target->GetVpOverP2AtZero(0)[k]==0 && sf[k]<-T && sf[k+1]>-T)
                                    expected+=volume(std::max(pf[l],sf[k]),std::min(pf[l+1],-T));
                    }
                    const real_t actual=column[j]*(sf[j+1]-sf[j])*norm*Hgrid;
                    if (!std::isfinite(actual) || std::abs(actual-expected)/H>0.005) {
                        PrintError("Pitch overlap: eps=%g j=%ld l=%ld error/H=%g",epsilon,(long)j,(long)l,(actual-expected)/H);
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

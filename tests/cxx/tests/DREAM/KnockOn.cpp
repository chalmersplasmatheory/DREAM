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
#include "KnockOnTestGrid.hpp"
#include "KnockOnTestReference.hpp"
#include <map>

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
    namespace K = DREAM::KnockOnUtilities;
    namespace Q = K::Quadrature;
    auto grid = KnockOnTest::ToroidalGrid(0.3,{-1,1});
    auto *radial = grid->GetRadialGrid();
    Q::OrbitGeometry geometry(radial->GetFluxSurfaceAverager(),0);
    Q::GaussLegendreRule rule(8);
    const real_t T = std::sqrt(geometry.GetBetaMax());
    struct Case { const char *name; real_t u, lo, hi, chi; };
    const Case cases[] = {
        {"passing",0.9,0.1,0.8,0.55}, {"trapped",0.6*T,0.1,0.9,0.55},
        {"near separatrix",T*(1-1e-4),0.15,0.85,0.1},
        {"near perpendicular",0.9,-0.1,0.1,0}, {"near forward",0.9,0.8,1,0.9999},
        {"pole",1,0.1,0.8,0.55}
    };
    for (const auto& c : cases) {
        const real_t a = K::EvaluateOrbitNumerator(0,c.chi,c.u,c.lo,c.hi,radial,geometry,rule);
        const real_t b = K::EvaluateOrbitNumerator(0,-c.chi,-c.u,c.lo,c.hi,radial,geometry,rule);
        const real_t d = K::EvaluateOrbitNumerator(0,-c.chi,c.u,-c.hi,-c.lo,radial,geometry,rule);
        const real_t scale = std::max({real_t(1e-100),std::abs(a),std::abs(b),std::abs(d)});
        if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(d)
            || std::abs(a-b)>1e-5*scale || std::abs(b-d)>1e-5*scale) {
            PrintError("Mirror symmetry: %s",c.name);
            return false;
        }
    }
    return true;
}

// Tests that the orbit-averaged delta calculation converges with increasing
// number of quadrature points. This is a numerical consistency test rather
// than a test of a physical invariant.
bool KnockOn::CheckDeltaQuadratureConvergence() {
    namespace K = DREAM::KnockOnUtilities;
    namespace Q = K::Quadrature;
    for (real_t epsilon : {0.05,0.3}) {
        auto grid = KnockOnTest::ToroidalGrid(epsilon,{-1,1});
        auto *radial = grid->GetRadialGrid();
        Q::OrbitGeometry geometry(radial->GetFluxSurfaceAverager(),0);
        const real_t T = std::sqrt(geometry.GetBetaMax());
        const Q::GaussLegendreRule theta(8), pitch(4), fineTheta(16), finePitch(8);
        const KnockOnTest::SurfaceReference volume(radial,0,8192);
        struct Case { const char *name; Q::Interval primary, secondary; real_t chi; };
        const Case cases[] = {
            {"negative passing",{-0.95,-0.85},{-1,-0.8},0.55},
            {"both trapped",{0.1*T,0.4*T},{0.1*T,0.8*T},0.55},
            {"crossing separatrix",{0.9*T,1.1*T},{0.8*T,1.2*T},0.55},
            {"crossing zero",{-0.1*T,0.1*T},{0,0.5*T},0},
            {"near perpendicular",{0.8,1},{-0.08,0.07},0},
            {"near forward",{0.8,1},{0.8,1},0.9999}
        };
        for (const auto& c : cases) {
            const real_t a = K::EvaluatePitchCellNumerator(0,c.chi,c.primary,c.secondary,radial,geometry,theta,pitch);
            const real_t b = K::EvaluatePitchCellNumerator(0,c.chi,c.primary,c.secondary,radial,geometry,fineTheta,finePitch);
            const real_t H = volume.Volume(c.primary.lower,c.primary.upper);
            const real_t error = std::abs(a-b)*(c.secondary.upper-c.secondary.lower)/H;
            // This is resolution consistency, not the independent accuracy
            // oracle: raw poloidal and cell-volume tests below supply that.
            if (!(H>0) || !std::isfinite(error) || a<0 || b<0 || error>0.005) {
                PrintError("Cell quadrature: %s epsilon=%g error=%g",c.name,epsilon,error);
                return false;
            }
        }
    }
    return true;
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
    auto grid = KnockOnTest::ToroidalGrid(0.3,KnockOnTest::UniformPitchFaces(8));
    return _checkLocalContributionConservationProperty(grid.get(),grid.get(),
        10*std::numeric_limits<real_t>::epsilon());
}

// Return true if the delta matrix satisfies the exact analytic normalization property:
//     \int dxi0 delta(xi0, xi01) = 1
// within the provided tolerance, on the provided grids.
bool KnockOn::_checkDeltaConservationProperty(
    const DREAM::FVM::Grid *grid_knockon, const DREAM::FVM::Grid *grid_primary,
    real_t successRelErrorThreshold
) {
    const auto *rg = grid_knockon->GetRadialGrid();
    namespace Q = DREAM::KnockOnUtilities::Quadrature;
    const Q::GaussLegendreRule theta(8), pitch(4);
    for (len_t ir=0; ir<rg->GetNr(); ++ir) {
        Q::OrbitGeometry geometry(rg->GetFluxSurfaceAverager(),ir);
        const auto *mg=grid_knockon->GetMomentumGrid(ir);
        std::vector<real_t> column(mg->GetNp2());
        for (real_t chi : {0.0,0.5,1.0})
            for (len_t l=0; l<grid_primary->GetNp2(ir); ++l) {
                const bool isVoid=grid_primary->GetVpOverP2AtZero(ir)[l]==0;
                const real_t rawNorm=DREAM::KnockOnUtilities::SetDeltaMatrixColumnOnGrid(
                    ir,chi,l,grid_knockon,grid_primary,column.data(),geometry,theta,pitch
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
    auto grid = KnockOnTest::ToroidalGrid(0.3,KnockOnTest::UniformPitchFaces(9));
    return _checkDeltaConservationProperty(grid.get(),grid.get(),1e-12);
}

// The Delta function in our formulation of the knock-on operator satisfies
// an exact analytic normalization property:
//     \int dxi0 delta(xi0, xi01) = 1.
// This test verifies that this property is reproduced numerically within
// the accuracy expected from the chosen quadrature resolution, for all
// cells that it is expected to hold (e.g. not in the negative-trapped region).
bool KnockOn::CheckDeltaConservationPropertyConstantB() {
    KnockOnTest::OwnedGrid grid(InitializeGridRCylPXi(1,1,7,1,0,2));
    return _checkDeltaConservationProperty(grid.get(),grid.get(),1e-12);
}

// Like CheckDeltaConservationProperty, but when secondary and primary
// electrons live on different grids.
bool KnockOn::CheckDeltaConservationPropertyDifferentGrids() {
    auto target = KnockOnTest::ToroidalGrid(0.05,{-1,-0.8,-0.5,-0.2,0.1,0.4,0.75,1});
    auto primary = KnockOnTest::ToroidalGrid(0.05,KnockOnTest::UniformPitchFaces(9));
    return _checkDeltaConservationProperty(target.get(),primary.get(),1e-12);
}

// In cylindrical geometry (constant B), the orbit-averaged delta function
// reduces to the analytically integrated gyro-averaged delta divided by
// the FVM xi-cell width. This test verifies that the general implementation
// reproduces this limit.
bool KnockOn::CheckCylindricalDeltaCalculation() {
    namespace K=DREAM::KnockOnUtilities;
    // For chi=1 in a cylinder the exact cell integral is the interval overlap.
    // Unlike, odd/even grids test the primary cell average, not a centre sample.
    for (const auto &sizes : {std::pair<len_t,len_t>{5,8},{8,5}}) {
        KnockOnTest::OwnedGrid target(InitializeGridRCylPXi(1,1,sizes.first,1,0,2));
        KnockOnTest::OwnedGrid primary(InitializeGridRCylPXi(1,1,sizes.second,1,0,2));
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

    auto owned = KnockOnTest::ToroidalGrid(0.3,KnockOnTest::UniformPitchFaces(12),1,2,0,2);
    auto *grid = owned.get();
    const len_t nr = grid->GetNr();
    auto *rg = grid->GetRadialGrid();
    auto *fsa = rg->GetFluxSurfaceAverager();

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
    KnockOnTest::OwnedGrid grid(InitializeGridRCylPXi(1,12,1,1,0,4));
    return _checkMollerSConservationProperty(grid.get(),grid.get());
}

/**
 * Reference uses a dense midpoint rule on one globally mapped orbit and the
 * elementary gyro-cone probability. It uses neither contact cuts nor the
 * production local-delta/metric helpers. Check raw values before normalization.
 */
bool KnockOn::CheckPoloidalNumerators(bool extended) {
    namespace K = DREAM::KnockOnUtilities;
    namespace Q = K::Quadrature;
    for (real_t u : {-1.0,1.0})
        for (real_t chi : {-1.0,-0.5,0.0,0.5,1.0}) {
            real_t sum=0;
            for (len_t j=0; j<4; ++j) {
                real_t ratio;
                sum += 0.5*K::AnalyticDelta::EvaluateLocalContribution(
                    -1+0.5*j,-0.5+0.5*j,u,1,chi,0,ratio);
            }
            if (sum!=1) return false;
        }
    for (real_t epsilon : {0.05,0.3}) {
        auto grid = KnockOnTest::ToroidalGrid(epsilon,{-1,1});
        auto *radial = grid->GetRadialGrid();
        Q::OrbitGeometry geometry(radial->GetFluxSurfaceAverager(),0);
        const real_t T = std::sqrt(geometry.GetBetaMax());
        const Q::GaussLegendreRule rule(8);
        struct Case { const char *name; real_t u, lo, hi, chi; };
        std::vector<Case> cases;
        if (extended) {
            for (real_t u : {0.6*T,T*(1-1e-4),T*(1+1e-4),(1+T)/2,1.0})
                for (real_t chi : {0.1,0.55,0.9999})
                    for (const auto cell : {Q::Interval{-1,1},Q::Interval{0.15,0.85}})
                        cases.push_back({"extended",u,cell.lower,cell.upper,chi});
        } else {
            cases = {
                {"trapped",0.6*T,0.15,0.85,0.55},
                {"below separatrix",T*(1-1e-4),0.15,0.85,0.1},
                {"above separatrix",T*(1+1e-4),0.15,0.85,0.55},
                {"passing full measure",(1+T)/2,-1,1,0.9999},
                {"narrow outgoing support",0.9,0.855,0.86,0.55},
                {"pole",1,0.15,0.85,0.55},
                {"near perpendicular",0.9,-0.1,0.1,0}
            };
        }
        const len_t n = extended ? 32768 : 8192;
        using Reference = KnockOnTest::OrbitReference;
        std::map<real_t,std::unique_ptr<Reference>> fine, coarse;
        for (const auto& c : cases) {
            if (!fine.count(c.u)) {
                fine[c.u] = std::make_unique<Reference>(radial,0,c.u,n);
                coarse[c.u] = std::make_unique<Reference>(radial,0,c.u,n/2);
            }
            const real_t ref = fine[c.u]->Cell(c.lo,c.hi,c.chi);
            const real_t check = coarse[c.u]->Cell(c.lo,c.hi,c.chi);
            const real_t H = fine[c.u]->Volume();
            const real_t value = K::EvaluateOrbitNumerator(0,c.chi,c.u,c.lo,c.hi,radial,geometry,rule);
            const real_t scale = (c.hi-c.lo)/H;
            if (!(H>0) || !std::isfinite(value) || value<0
                || std::abs(ref-check)*scale>2e-4 || std::abs(value-ref)*scale>1e-3) {
                PrintError("Raw orbit: %s eps=%g u=%.16g chi=%g error=%g reference convergence=%g",
                    c.name,epsilon,c.u,c.chi,(value-ref)*scale,(ref-check)*scale);
                return false;
            }
        }
    }
    return true;
}

bool KnockOn::CheckPitchCellAverages(bool extended) {
    namespace K=DREAM::KnockOnUtilities;
    namespace Q=K::Quadrature;
    // The default faces have no sign symmetry: negative/positive boundary
    // rows have different widths, and one primary cell crosses zero.
    const std::vector<real_t> sf = extended ? KnockOnTest::UniformPitchFaces(40)
        : std::vector<real_t>{-1,-.82,-.55,-.2,.12,.43,.72,.86,1};
    const std::vector<real_t> pf = extended ? KnockOnTest::UniformPitchFaces(31)
        : std::vector<real_t>{-1,-.79,-.51,-.11,.17,.57,.83,1};
    const Q::GaussLegendreRule theta(8), pitch(4);
    for (real_t epsilon : {0.05,0.3}) {
        auto target = KnockOnTest::ToroidalGrid(epsilon,sf);
        auto primary = KnockOnTest::ToroidalGrid(epsilon,pf);
        auto *rg=primary->GetRadialGrid();
        auto *fsa=rg->GetFluxSurfaceAverager();
        Q::OrbitGeometry geometry(fsa,0);
        const real_t T=std::sqrt(geometry.GetBetaMax());
        // Geometry is independent of the primary/secondary intervals. Cache
        // each reference node grid once, including a separate convergence grid.
        const len_t n = extended ? 65536 : 16384;
        const KnockOnTest::SurfaceReference reference(rg,0,n), coarse(rg,0,n/2);
        auto volume = [&](real_t lo, real_t hi) { return reference.Volume(lo,hi); };
        std::vector<real_t> column(sf.size()-1);
        for (len_t l=0; l+1<pf.size(); ++l) {
            if (primary->GetVpOverP2AtZero(0)[l]==0) continue;
            const real_t H=volume(pf[l],pf[l+1]);
            if (!(H>0) || std::abs(coarse.Volume(pf[l],pf[l+1])/H-1)>1e-4) return false;
            const real_t Hgrid=(pf[l+1]-pf[l])*primary->GetVpOverP2AtZero(0)[l];
            for (real_t chi : {0.0,0.55,1.0}) {
                const real_t norm=K::SetDeltaMatrixColumnOnGrid(
                    0,chi,l,target.get(),primary.get(),column.data(),geometry,theta,pitch
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


// Called only by the explicitly built extended executable, never by Run().
bool KnockOn::RunExtended() {
    bool success = CheckPoloidalNumerators(true);
    success = CheckPitchCellAverages(true) && success;
    auto target = KnockOnTest::ToroidalGrid(0.3,KnockOnTest::UniformPitchFaces(40));
    auto primary = KnockOnTest::ToroidalGrid(0.3,KnockOnTest::UniformPitchFaces(31));
    success = _checkDeltaConservationProperty(target.get(),primary.get(),1e-12) && success;
    if (success) PrintOK("Extended raw quadrature, cell-volume/overlap and unlike-grid checks passed.");
    return success;
}

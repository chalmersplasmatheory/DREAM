#include "MollerBoltzmannOperator.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include "DREAM/Equations/Kinetic/MollerBoltzmannOperator.hpp"
#include "DREAM/Equations/Kinetic/MollerDeltaAngleKernel.hpp"
#include "DREAM/Equations/Kinetic/MollerEnergyKernel.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "KnockOnTestHelpers.hpp"

using namespace DREAMTESTS::_DREAM;

namespace {

void ApplyTerm(
    DREAM::MollerBoltzmannOperator &op, const DREAM::FVM::Grid *grid_knockon, const real_t ntot,
    std::vector<real_t> &out
) {
    const len_t Nr = grid_knockon->GetNr();
    std::vector<real_t> ntot_arr(Nr, ntot);

    const len_t NK = grid_knockon->GetNCells();
    out.assign(NK, 0.0);
    op.SetVectorElements(out.data(), ntot_arr.data());
}

}  // anonymous namespace

bool MollerBoltzmannOperator::Run(bool) {
    bool success = true;

    if (CheckBO_LinearityInFPrimary())
        this->PrintOK("Assembled contribution is linear in f_primary.");
    else {
        success = false;
        this->PrintError("Test failed: assembled contribution is not linear in f_primary.");
    }

    if (CheckBO_JacobianFiniteDifferenceNt())
        this->PrintOK("Jacobian w.r.t. n_tot matches finite-difference derivative.");
    else {
        success = false;
        this->PrintError("Test failed: Jacobian w.r.t. n_tot does not match FD.");
    }

    if (CheckBO_GlobalProductionIdentity())
        this->PrintOK("Global production identity matches Moller-S prediction.");
    else {
        success = false;
        this->PrintError("Test failed: global production identity does not match prediction.");
    }

    if (CheckBO_HotRunawayGlobalProductionIdentity())
        this->PrintOK("hot<-runaway (cross-grid) global identity passed.");
    else {
        success = false;
        this->PrintError("Test failed: hot<-runaway (cross-grid) global identity.");
    }

    if (CheckBO_TimeCachingRegression())
        this->PrintOK("Time caching regression test passed.");
    else {
        success = false;
        this->PrintError("Test failed: time caching regression.");
    }

    if (CheckBO_NonNegativity())
        this->PrintOK("Non-negativity sanity check passed.");
    else {
        success = false;
        this->PrintError("Test failed: non-negativity sanity check.");
    }

    if (CheckBO_RadiusLocality())
        this->PrintOK("Radius locality test passed.");
    else {
        success = false;
        this->PrintError("Test failed: radius locality.");
    }

    return success;
}

bool MollerBoltzmannOperator::CheckBO_LinearityInFPrimary() {
    const real_t tol = 1e-10;

    len_t nr = 2;
    len_t np = 6;
    len_t nxi = 20;

    len_t ntheta_interp = 50;
    len_t nrProfiles = 8;
    real_t pMin = 0;
    real_t pMax = 3;

    real_t B0 = 1.0;
    auto *gridF = InitializeFluidGrid(nr, B0);
    auto *gridK = InitializeGridGeneralRPXi(nr, np, nxi, ntheta_interp, nrProfiles, pMin, pMax);
    auto *gridP = InitializeGridGeneralRPXi(nr, np, nxi, ntheta_interp, nrProfiles, pMin, pMax);

    len_t id_ntot, id_E, id_f;
    auto *uqh = BuildUQH_Minimal(gridF, gridP, id_ntot, id_E, id_f);

    const real_t pCutoff = gridK->GetMomentumGrid(0)->GetP1_f(2);

    const real_t scaleFactor = 1.0;
    const len_t nXiStars = 60;
    const len_t nIntPts = 80;
    DREAM::KnockOnUtilities::orbit_integration_method integrationMethod = DREAM::KnockOnUtilities::MIDPOINT_RULE;

    auto *energy = new DREAM::MollerEnergyKernel(gridK, gridP, pCutoff);
    auto *angle = new DREAM::MollerDeltaAngleKernel(gridK, gridP, pCutoff, nXiStars, nIntPts, integrationMethod);
    DREAM::MollerBoltzmannOperator op(gridK, gridP, uqh, id_f, energy, angle, scaleFactor);

    const len_t NP = gridP->GetNCells();
    real_t *fA = new real_t[NP];
    real_t *fB = new real_t[NP];
    real_t *fAB = new real_t[NP];

    for (len_t i = 0; i < NP; i++) {
        fA[i] = 1e20 + 1e18 * (real_t)i;
        fB[i] = 1e20 + 2e18 * (real_t)i;
        fAB[i] = fA[i] + fB[i];
    }

    const real_t ntot = 1.0;
    const len_t NK = gridK->GetNCells();
    std::vector<real_t> RA, RB, RAB;

    SetPreviousUnknownData(uqh, id_f, fA, 0.0);
    op.Rebuild(0.0, 1.0, uqh);
    ApplyTerm(op, gridK, ntot, RA);

    SetPreviousUnknownData(uqh, id_f, fB, 1.0);
    op.Rebuild(1.0, 1.0, uqh);
    ApplyTerm(op, gridK, ntot, RB);

    SetPreviousUnknownData(uqh, id_f, fAB, 2.0);
    op.Rebuild(2.0, 1.0, uqh);
    ApplyTerm(op, gridK, ntot, RAB);

    bool success = true;
    for (len_t q = 0; q < NK; q++) {
        const real_t lhs = RAB[q];
        const real_t rhs = RA[q] + RB[q];
        const real_t diff = fabs(lhs - rhs);
        if (diff > tol * (1 + fabs(lhs))) {
            this->PrintError(
                "Linearity failed at q=%ld: lhs=%.8g rhs=%.8g diff=%.8g", (long)q, lhs, rhs, diff
            );
            success = false;
            break;
        }
    }

    delete[] fA;
    delete[] fB;
    delete[] fAB;

    delete uqh;
    delete gridF;
    delete gridK;
    delete gridP;

    delete energy;
    delete angle;

    return success;
}

bool MollerBoltzmannOperator::CheckBO_JacobianFiniteDifferenceNt() {
    const real_t tol = 1e-6;
    const real_t eps = 1e-7;

    len_t nr = 2;
    len_t np = 6;
    len_t nxi = 10;

    len_t ntheta_interp = 50;
    len_t nrProfiles = 8;
    real_t pMin = 0;
    real_t pMax = 3;

    real_t B0 = 1.0;
    auto *gridF = InitializeFluidGrid(nr, B0);
    auto *gridK = InitializeGridGeneralRPXi(nr, np, nxi, ntheta_interp, nrProfiles, pMin, pMax);
    auto *gridP = InitializeGridGeneralRPXi(nr, np, nxi, ntheta_interp, nrProfiles, pMin, pMax);

    len_t id_ntot, id_E, id_f;
    auto *uqh = BuildUQH_Minimal(gridF, gridP, id_ntot, id_E, id_f);

    const real_t pCutoff = gridK->GetMomentumGrid(0)->GetP1_f(2);
    constexpr real_t scaleFactor = 1.0;
    constexpr len_t nXiStars = 80;
    constexpr len_t nIntPts = 80;
    DREAM::KnockOnUtilities::orbit_integration_method integrationMethod = DREAM::KnockOnUtilities::MIDPOINT_RULE;

    auto *energy = new DREAM::MollerEnergyKernel(gridK, gridP, pCutoff);
    auto *angle = new DREAM::MollerDeltaAngleKernel(gridK, gridP, pCutoff, nXiStars, nIntPts, integrationMethod);
    DREAM::MollerBoltzmannOperator op(gridK, gridP, uqh, id_f, energy, angle, scaleFactor);

    const len_t NP = gridP->GetNCells();
    real_t *f0 = new real_t[NP];
    for (len_t i = 0; i < NP; i++) {
        f0[i] = 1e20;
    }
    SetPreviousUnknownData(uqh, id_f, f0, 0.0);
    op.Rebuild(0.0, 1.0, uqh);

    const len_t NK = gridK->GetNCells();
    real_t *R0 = new real_t[NK];
    real_t *R1 = new real_t[NK];
    for (len_t q = 0; q < NK; q++) {
        R0[q] = 0.0;
        R1[q] = 0.0;
    }

    real_t *nt = new real_t[nr];
    for (len_t ir = 0; ir < nr; ir++) nt[ir] = 1.0;

    op.SetVectorElements(R0, nt);

    const len_t ir0 = 0;
    nt[ir0] += eps;
    op.SetVectorElements(R1, nt);

    std::vector<real_t> dR(NK, 0.0);
    for (len_t q = 0; q < NK; q++) {
        dR[q] = (R1[q] - R0[q]) / eps;
    }
    DREAM::FVM::Matrix jac(NK, NK, 1);
    jac.Zero();
    jac.PartialAssemble();
    op.SetJacobianBlock(id_ntot, id_ntot, &jac, nullptr);
    jac.Assemble();

    bool success = true;
    for (len_t q = 0; q < NK; q++) {
        const real_t J = jac.GetElement(q, ir0);
        const real_t diff = fabs(J - dR[q]);
        if (diff > tol * (1 + fabs(J) + fabs(dR[q]))) {
            success = false;
            break;
        }
    }

    delete[] f0;
    delete[] R0;
    delete[] R1;
    delete[] nt;

    delete uqh;
    delete gridF;
    delete gridK;
    delete gridP;

    delete energy;
    delete angle;

    return success;
}

bool MollerBoltzmannOperator::CheckBO_GlobalProductionIdentity() {
    const real_t rtol = 1e-14;
    const real_t atol = 1e-14;

    len_t nr = 2;
    len_t np = 15;
    len_t nxi = 30;

    len_t ntheta_interp = 50;
    len_t nrProfiles = 8;
    real_t pMin = 0;
    real_t pMax = 4;

    real_t B0 = 1.0;
    auto *gridF = InitializeFluidGrid(nr, B0);
    auto *gridK = InitializeGridGeneralRPXi(nr, np, nxi, ntheta_interp, nrProfiles, pMin, pMax);
    auto *gridP = InitializeGridGeneralRPXi(nr, np, nxi, ntheta_interp, nrProfiles, pMin, pMax);

    len_t id_ntot, id_E, id_f;
    auto *uqh = BuildUQH_Minimal(gridF, gridP, id_ntot, id_E, id_f);

    const real_t pCutoff = gridK->GetMomentumGrid(0)->GetP1_f(2);

    const real_t scaleFactor = 1.0;
    const len_t nXiStars = 80;
    const len_t nIntPts = 80;
    DREAM::KnockOnUtilities::orbit_integration_method integrationMethod = DREAM::KnockOnUtilities::MIDPOINT_RULE;

    auto *energy = new DREAM::MollerEnergyKernel(gridK, gridP, pCutoff);
    auto *angle = new DREAM::MollerDeltaAngleKernel(gridK, gridP, pCutoff, nXiStars, nIntPts, integrationMethod);
    DREAM::MollerBoltzmannOperator op(gridK, gridP, uqh, id_f, energy, angle, scaleFactor);

    const len_t NP = gridP->GetNCells();
    real_t *f = new real_t[NP];
    for (len_t i = 0; i < NP; i++) {
        f[i] = 1e20;
    }
    SetPreviousUnknownData(uqh, id_f, f, 0.0);
    op.Rebuild(0.0, 1.0, uqh);

    std::vector<real_t> vec;
    ApplyTerm(op, gridK, /*ntot*/ 1.0, vec);

    const real_t LHS = IntegrateTotalProductionOverKnockonGrid(gridK, vec.data());

    const real_t RHS0 = PredictTotalProductionFromMollerS(gridK, gridP, f, pCutoff);
    const real_t RHS = scaleFactor * RHS0;

    const real_t diff = fabs(LHS - RHS);
    const bool success = (diff <= atol + rtol * (fabs(LHS) + fabs(RHS)));

    if (!success) {
        this->PrintError("Global identity failed:");
        this->PrintError("  LHS: %.8g", LHS);
        this->PrintError("  RHS: %.8g", RHS);
        this->PrintError("  diff: %.8g", diff);
    }

    delete[] f;

    delete uqh;
    delete gridF;
    delete gridK;
    delete gridP;

    delete energy;
    delete angle;

    return success;
}

bool MollerBoltzmannOperator::CheckBO_HotRunawayGlobalProductionIdentity() {
    const real_t rtol = 1e-14;
    const real_t atol = 1e-14;
    const real_t absTolNonNeg = 1e-13;

    const len_t nr = 2;

    const len_t npHot = 18;
    const len_t nxiHot = 12;

    const len_t npRe = 22;
    const len_t nxiRe = 14;

    const len_t ntheta_interp = 50;
    const len_t nrProfiles = 8;

    const real_t pMin = 0.0;
    const real_t pMaxHot = 4.0;
    const real_t pMaxRe = 20.0;

    const real_t B0 = 1.0;
    auto *gridF = InitializeFluidGrid(nr, B0);

    auto *gridHot =
        InitializeGridGeneralRPXi(nr, npHot, nxiHot, ntheta_interp, nrProfiles, pMin, pMaxHot);
    auto *gridRe =
        InitializeGridGeneralRPXi(nr, npRe, nxiRe, ntheta_interp, nrProfiles, pMin, pMaxRe);

    len_t id_ntot, id_E, id_f_re;
    auto *uqh = BuildUQH_MinimalFRe(gridF, gridRe, id_ntot, id_E, id_f_re);

    const real_t pCutoff = gridHot->GetMomentumGrid(0)->GetP1_f(2);

    const real_t ntot = 1.0;
    const real_t scaleFactor = 1.0;
    const len_t nXiStars = 80;
    const len_t nIntPts = 80;
    DREAM::KnockOnUtilities::orbit_integration_method integrationMethod = DREAM::KnockOnUtilities::MIDPOINT_RULE;

    auto *energy = new DREAM::MollerEnergyKernel(gridHot, gridRe, pCutoff);
    auto *angle = new DREAM::MollerDeltaAngleKernel(gridHot, gridRe, pCutoff, nXiStars, nIntPts, integrationMethod);
    DREAM::MollerBoltzmannOperator op(gridHot, gridRe, uqh, id_f_re, energy, angle, scaleFactor);

    const len_t NRe = gridRe->GetNCells();
    real_t *f = new real_t[NRe];
    for (len_t idx = 0; idx < NRe; idx++) {
        f[idx] = 1e20 * (1.0 + 0.1 * (real_t)(idx % 7));
    }
    const real_t dt = 1.0;
    const real_t t = 1.0;
    SetPreviousUnknownData(uqh, id_f_re, f, t - dt);

    op.Rebuild(t, dt, uqh);

    std::vector<real_t> vec;
    ApplyTerm(op, gridHot, ntot, vec);

    bool success = true;

    {
        real_t minVal = 1e300;
        len_t minIdx = 0;
        for (len_t q = 0; q < (len_t)vec.size(); q++) {
            if (vec[q] < minVal) {
                minVal = vec[q];
                minIdx = q;
            }
        }
        if (minVal < -absTolNonNeg) {
            this->PrintError(
                "hot<-re non-negativity failed: min(vec)=%.8g at q=%ld", minVal, minIdx
            );
            success = false;
        }
    }

    const real_t LHS = IntegrateTotalProductionOverKnockonGrid(gridHot, vec.data());

    const real_t RHS0 = PredictTotalProductionFromMollerS(gridHot, gridRe, f, pCutoff);
    const real_t RHS = scaleFactor * ntot * RHS0;

    const real_t diff = fabs(LHS - RHS);
    if (diff > atol + rtol * (fabs(LHS) + fabs(RHS))) {
        this->PrintError("hot<-re global identity failed:");
        this->PrintError("  LHS: %.8g", LHS);
        this->PrintError("  RHS: %.8g", RHS);
        this->PrintError("  diff: %.8g", diff);
        success = false;
    }

    delete[] f;

    delete uqh;
    delete gridF;
    delete gridHot;
    delete gridRe;

    delete energy;
    delete angle;

    return success;
}

bool MollerBoltzmannOperator::CheckBO_TimeCachingRegression() {
    const real_t rtol = 1e-12;

    len_t nr = 2;
    len_t np = 6;
    len_t nxi = 12;

    len_t ntheta_interp = 50;
    len_t nrProfiles = 8;
    real_t pMin = 0;
    real_t pMax = 3;

    real_t B0 = 1.0;
    auto *gridF = InitializeFluidGrid(nr, B0);
    auto *gridK = InitializeGridGeneralRPXi(nr, np, nxi, ntheta_interp, nrProfiles, pMin, pMax);
    auto *gridP = InitializeGridGeneralRPXi(nr, np, nxi, ntheta_interp, nrProfiles, pMin, pMax);

    len_t id_ntot, id_E, id_f;
    auto *uqh = BuildUQH_Minimal(gridF, gridP, id_ntot, id_E, id_f);

    const real_t pCutoff = gridK->GetMomentumGrid(0)->GetP1_f(2);
    const real_t scaleFactor = 1.0;
    const len_t nXiStars = 60;
    const len_t nIntPts = 60;
    DREAM::KnockOnUtilities::orbit_integration_method integrationMethod = DREAM::KnockOnUtilities::MIDPOINT_RULE;

    auto *energy = new DREAM::MollerEnergyKernel(gridK, gridP, pCutoff);
    auto *angle = new DREAM::MollerDeltaAngleKernel(gridK, gridP, pCutoff, nXiStars, nIntPts, integrationMethod);
    DREAM::MollerBoltzmannOperator op(gridK, gridP, uqh, id_f, energy, angle, scaleFactor);

    const len_t NP = gridP->GetNCells();
    real_t *fA = new real_t[NP];
    real_t *fB = new real_t[NP];
    for (len_t idx = 0; idx < NP; idx++) {
        fA[idx] = 1e20;
        fB[idx] = 2e20;
    }

    std::vector<real_t> vecA, vecSameT, vecNewT;
    const real_t ntot = 1.0;

    const real_t dt = 1.0;
    const real_t t1 = 1.0;
    const real_t t2 = 2.0;

    SetPreviousUnknownData(uqh, id_f, fA, t1 - dt);
    op.Rebuild(t1, dt, uqh);
    ApplyTerm(op, gridK, ntot, vecA);

    SetPreviousUnknownData(uqh, id_f, fB, t1);
    op.Rebuild(t1, dt, uqh);
    ApplyTerm(op, gridK, ntot, vecSameT);

    op.Rebuild(t2, dt, uqh);
    ApplyTerm(op, gridK, ntot, vecNewT);

    real_t maxRef = 0.0, maxDiffSame = 0.0, maxDiffNew = 0.0;
    for (len_t q = 0; q < (len_t)vecA.size(); q++) {
        maxRef = std::max(maxRef, fabs(vecA[q]));
        maxDiffSame = std::max(maxDiffSame, fabs(vecSameT[q] - vecA[q]));
        maxDiffNew = std::max(maxDiffNew, fabs(vecNewT[q] - 2.0 * vecA[q]));
    }

    bool success = true;
    if (maxDiffSame > rtol * (1.0 + maxRef)) {
        this->PrintError("Time-caching failed: vec changed despite same t.");
        this->PrintError("  maxDiffSame=%.8g maxRef=%.8g", maxDiffSame, maxRef);
        success = false;
    }
    if (maxDiffNew > rtol * (1.0 + 2.0 * maxRef)) {
        this->PrintError("Time-caching failed: vec after t-change not consistent with f scaling.");
        this->PrintError("  maxDiffNew=%.8g maxRef=%.8g", maxDiffNew, maxRef);
        success = false;
    }

    delete[] fA;
    delete[] fB;

    delete uqh;
    delete gridF;
    delete gridK;
    delete gridP;

    delete energy;
    delete angle;

    return success;
}

bool MollerBoltzmannOperator::CheckBO_NonNegativity() {
    const real_t absTol = 1e-13;

    len_t nr = 2;
    len_t np = 10;
    len_t nxi = 16;

    len_t ntheta_interp = 50;
    len_t nrProfiles = 8;
    real_t pMin = 0;
    real_t pMax = 3;

    real_t B0 = 1.0;
    auto *gridF = InitializeFluidGrid(nr, B0);
    auto *gridK = InitializeGridGeneralRPXi(nr, np, nxi, ntheta_interp, nrProfiles, pMin, pMax);
    auto *gridP = InitializeGridGeneralRPXi(nr, np, nxi, ntheta_interp, nrProfiles, pMin, pMax);

    len_t id_ntot, id_E, id_f;
    auto *uqh = BuildUQH_Minimal(gridF, gridP, id_ntot, id_E, id_f);

    const real_t pCutoff = gridK->GetMomentumGrid(0)->GetP1_f(2);
    const real_t scaleFactor = 1.0;
    const len_t nXiStars = 60;
    const len_t nIntPts = 60;
    DREAM::KnockOnUtilities::orbit_integration_method integrationMethod = DREAM::KnockOnUtilities::MIDPOINT_RULE;

    auto *energy = new DREAM::MollerEnergyKernel(gridK, gridP, pCutoff);
    auto *angle = new DREAM::MollerDeltaAngleKernel(gridK, gridP, pCutoff, nXiStars, nIntPts, integrationMethod);
    DREAM::MollerBoltzmannOperator op(gridK, gridP, uqh, id_f, energy, angle, scaleFactor);

    const len_t NP = gridP->GetNCells();
    real_t *f = new real_t[NP];
    for (len_t idx = 0; idx < NP; idx++) {
        f[idx] = 1e20;
    }
    SetPreviousUnknownData(uqh, id_f, f, 0.0);
    op.Rebuild(0.0, 1.0, uqh);

    std::vector<real_t> vec;
    ApplyTerm(op, gridK, 1.0, vec);

    bool success = true;

    real_t minVal = 1e300;
    len_t minIdx = 0;
    for (len_t q = 0; q < (len_t)vec.size(); q++) {
        if (vec[q] < minVal) {
            minVal = vec[q];
            minIdx = q;
        }
    }

    if (minVal < -absTol) {
        this->PrintError("Non-negativity failed: min(vec)=%.8g at q=%ld", minVal, (long)minIdx);
        success = false;
    }

    // Check source==0 in void region (Vp==0)
    len_t offset = 0;
    for (len_t ir = 0; ir < gridK->GetNr(); ir++) {
        for (len_t i = 0; i < gridK->GetNp1(ir); i++) {
            for (len_t j = 0; j < gridK->GetNp2(ir); j++) {
                const len_t idx = offset + gridK->GetNp1(ir) * j + i;
                if (gridK->GetVp(ir, i, j) == 0 && vec[idx] != 0) {
                    this->PrintError(
                        "Non-negativity failed: non-zero in void region (Vp==0): "
                        "vec=%.16g at ir=%ld, i=%ld, j=%ld",
                        vec[idx], ir, i, j
                    );
                    success = false;
                }
            }
        }
        offset += gridK->GetMomentumGrid(ir)->GetNCells();
    }

    delete[] f;

    delete uqh;
    delete gridF;
    delete gridK;
    delete gridP;

    delete energy;
    delete angle;

    return success;
}

bool MollerBoltzmannOperator::CheckBO_RadiusLocality() {
    const real_t absTol = 1e-14;
    const real_t relTol = 1e-12;

    len_t nr = 3;
    len_t np = 8;
    len_t nxi = 12;

    len_t ntheta_interp = 50;
    len_t nrProfiles = 8;
    real_t pMin = 0;
    real_t pMax = 3;

    real_t B0 = 1.0;
    auto *gridF = InitializeFluidGrid(nr, B0);
    auto *gridK = InitializeGridGeneralRPXi(nr, np, nxi, ntheta_interp, nrProfiles, pMin, pMax);
    auto *gridP = InitializeGridGeneralRPXi(nr, np, nxi, ntheta_interp, nrProfiles, pMin, pMax);

    len_t id_ntot, id_E, id_f;
    auto *uqh = BuildUQH_Minimal(gridF, gridP, id_ntot, id_E, id_f);

    const real_t pCutoff = gridK->GetMomentumGrid(0)->GetP1_f(2);
    const real_t scaleFactor = 1.0;
    const len_t nXiStars = 60;
    const len_t nIntPts = 60;
    DREAM::KnockOnUtilities::orbit_integration_method integrationMethod = DREAM::KnockOnUtilities::MIDPOINT_RULE;

    auto *energy = new DREAM::MollerEnergyKernel(gridK, gridP, pCutoff);
    auto *angle = new DREAM::MollerDeltaAngleKernel(gridK, gridP, pCutoff, nXiStars, nIntPts, integrationMethod);
    DREAM::MollerBoltzmannOperator op(gridK, gridP, uqh, id_f, energy, angle, scaleFactor);

    const len_t NP = gridP->GetNCells();
    real_t *f = new real_t[NP];
    for (len_t idx = 0; idx < NP; idx++) {
        f[idx] = 1e20;
    }
    SetPreviousUnknownData(uqh, id_f, f, 0.0);
    op.Rebuild(0.0, 1.0, uqh);

    std::vector<real_t> vecAll;
    ApplyTerm(op, gridK, 1.0, vecAll);

    const len_t NK = gridK->GetNCells();
    std::vector<real_t> vecLocal(NK, 0.0);
    std::vector<real_t> ntot_arr(nr, 0.0);

    const len_t ir0 = 1;
    ntot_arr[ir0] = 1.0;
    op.SetVectorElements(vecLocal.data(), ntot_arr.data());

    bool success = true;
    len_t offset = 0;
    for (len_t ir = 0; ir < nr; ir++) {
        const auto *mg = gridK->GetMomentumGrid(ir);
        const len_t nCells = mg->GetNCells();
        const len_t start = offset;
        const len_t end = offset + nCells;

        if (ir != ir0) {
            real_t maxAbs = 0.0;
            for (len_t q = start; q < end; q++) {
                maxAbs = std::max(maxAbs, fabs(vecLocal[q]));
            }
            if (maxAbs > absTol) {
                this->PrintError(
                    "Radius locality failed: ir=%ld leakage maxAbs=%.8g", (long)ir, maxAbs
                );
                success = false;
                break;
            }
        } else {
            real_t maxDiff = 0.0;
            real_t maxRef = 0.0;
            for (len_t q = start; q < end; q++) {
                maxDiff = std::max(maxDiff, fabs(vecLocal[q] - vecAll[q]));
                maxRef = std::max(maxRef, fabs(vecAll[q]));
            }
            const real_t thresh = absTol + relTol * (1.0 + maxRef);
            if (maxDiff > thresh) {
                this->PrintError("Radius locality failed: ir0 block mismatch.");
                this->PrintError(
                    "  maxDiff=%.8g thresh=%.8g (maxRef=%.8g)", maxDiff, thresh, maxRef
                );
                success = false;
                break;
            }
        }

        offset += nCells;
    }

    delete[] f;

    delete uqh;
    delete gridF;
    delete gridK;
    delete gridP;

    delete energy;
    delete angle;

    return success;
}

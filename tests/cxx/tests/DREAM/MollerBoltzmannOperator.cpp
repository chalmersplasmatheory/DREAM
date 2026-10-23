#include "MollerBoltzmannOperator.hpp"

#include <algorithm>
#include <cmath>
#include <vector>
#include <memory>

#include "DREAM/Equations/Kinetic/MollerBoltzmannOperator.hpp"
#include "DREAM/Equations/Fluid/MollerProductionRateTerm.hpp"
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

    if (CheckBO_MatrixAssemblyAndFluidIdentity())
        this->PrintOK("Kinetic/fluid matrix assembly and production identity passed.");
    else {
        success = false;
        this->PrintError("Test failed: kinetic/fluid matrix assembly or production identity.");
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

    if (CheckBO_IsotropicZeroCurrentAndReflection())
        this->PrintOK("Assembled source preserves zero isotropic current and reflection.");
    else {
        success = false;
        this->PrintError("Test failed: assembled isotropic current/reflection.");
    }

    return success;
}

bool MollerBoltzmannOperator::CheckBO_IsotropicZeroCurrentAndReflection() {
    std::unique_ptr<DREAM::FVM::Grid> fluid(InitializeFluidGrid(1, 1.0));
    std::unique_ptr<DREAM::FVM::Grid> target(InitializeGridRCylPXi(1, 20, 45, 1.0, 0, 2));
    std::unique_ptr<DREAM::FVM::Grid> primary(InitializeGridRCylPXi(1, 30, 30, 1.0, 0, 3));
    len_t idN, idE, idF;
    std::unique_ptr<DREAM::FVM::UnknownQuantityHandler> uqh(
        BuildUQH_Minimal(fluid.get(), primary.get(), idN, idE, idF));
    DREAM::MollerEnergyKernel energy(target.get(), primary.get(), 0.1);
    DREAM::MollerDeltaAngleKernel angle(target.get(), primary.get(), 0.1,
        80, 8);
    DREAM::MollerBoltzmannOperator term(target.get(), primary.get(), uqh.get(),
        idF, &energy, &angle, 1.0);
    auto *mgP = primary->GetMomentumGrid(0), *mgK = target->GetMomentumGrid(0);
    std::vector<real_t> f(primary->GetNCells()), output[3];
    const real_t theta = 20000.0/510998.95;
    const real_t amplitudes[] = {0, 0.1, -0.1};
    for (len_t s=0; s<3; ++s) {
        for (len_t l=0; l<mgP->GetNp2(); ++l) for (len_t k=0; k<mgP->GetNp1(); ++k) {
            const real_t p = mgP->GetP1(k);
            f[l*mgP->GetNp1()+k] = 1e20*exp(-p*p/((sqrt(1+p*p)+1)*theta))
                *(1+amplitudes[s]*mgP->GetXi0(0,l));
        }
        SetPreviousUnknownData(uqh.get(), idF, f.data(), (real_t)s);
        term.Rebuild(s+1, 1, uqh.get());
        ApplyTerm(term, target.get(), 1, output[s]);
    }
    real_t gain = 0, drift = 0, reflectionError = 0;
    for (len_t j=0; j<mgK->GetNp2(); ++j) for (len_t i=0; i<mgK->GetNp1(); ++i) {
        const len_t q = j*mgK->GetNp1()+i;
        const len_t reflected = (mgK->GetNp2()-1-j)*mgK->GetNp1()+i;
        const real_t measure = mgK->GetDp1(i)*mgK->GetDp2(j)*target->GetVp(0,i,j);
        const real_t p = mgK->GetP1(i);
        gain += measure*output[0][q];
        drift += measure*output[0][q]*p/sqrt(1+p*p)*mgK->GetXi0(0,j);
        reflectionError += measure*fabs(output[1][q]-output[2][reflected]);
    }
    return std::isfinite(gain) && gain>0 && std::isfinite(drift)
        && std::isfinite(reflectionError) && fabs(drift)/gain<1e-12
        && reflectionError/gain<1e-12;
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
    const len_t poloidalOrder = 8;

    auto *energy = new DREAM::MollerEnergyKernel(gridK, gridP, pCutoff);
    auto *angle = new DREAM::MollerDeltaAngleKernel(gridK, gridP, pCutoff, nXiStars, poloidalOrder);
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
    constexpr len_t poloidalOrder = 8;

    auto *energy = new DREAM::MollerEnergyKernel(gridK, gridP, pCutoff);
    auto *angle = new DREAM::MollerDeltaAngleKernel(gridK, gridP, pCutoff, nXiStars, poloidalOrder);
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
    const len_t poloidalOrder = 8;

    auto *energy = new DREAM::MollerEnergyKernel(gridK, gridP, pCutoff);
    auto *angle = new DREAM::MollerDeltaAngleKernel(gridK, gridP, pCutoff, nXiStars, poloidalOrder);
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
    const len_t poloidalOrder = 8;

    auto *energy = new DREAM::MollerEnergyKernel(gridHot, gridRe, pCutoff);
    auto *angle = new DREAM::MollerDeltaAngleKernel(gridHot, gridRe, pCutoff, nXiStars, poloidalOrder);
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
    const len_t poloidalOrder = 8;

    auto *energy = new DREAM::MollerEnergyKernel(gridK, gridP, pCutoff);
    auto *angle = new DREAM::MollerDeltaAngleKernel(gridK, gridP, pCutoff, nXiStars, poloidalOrder);
    DREAM::MollerBoltzmannOperator op(gridK, gridP, uqh, id_f, energy, angle, scaleFactor);

    const len_t NP = gridP->GetNCells();
    real_t *fA = new real_t[NP];
    real_t *fB = new real_t[NP];
    for (len_t idx = 0; idx < NP; idx++) {
        fA[idx] = 1e20;
        fB[idx] = 2e20;
    }

    std::vector<real_t> vecA, vecSameT, vecNewT, vecRollback;
    const real_t ntot = 1.0;

    const real_t dt = 0.5;
    const real_t t1 = 0.5;
    const real_t t2 = 1.0;

    SetPreviousUnknownData(uqh, id_f, fA, t1 - dt);
    op.Rebuild(t1, dt, uqh);
    ApplyTerm(op, gridK, ntot, vecA);

    // A repeated Newton rebuild must use the same previous distribution.
    op.Rebuild(t1, dt, uqh);
    ApplyTerm(op, gridK, ntot, vecSameT);

    SetPreviousUnknownData(uqh, id_f, fB, t1);
    op.Rebuild(t2, dt, uqh);
    ApplyTerm(op, gridK, ntot, vecNewT);

    // Adaptive full-step trial: same endpoint as the second half-step,
    // but the previous solution has been restored to t=0.
    SetPreviousUnknownData(uqh, id_f, fA, 0.0);
    op.Rebuild(t2, 2*dt, uqh);
    ApplyTerm(op, gridK, ntot, vecRollback);

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
    for (len_t q = 0; q < vecRollback.size(); ++q) {
        if (!std::isfinite(vecRollback[q]) ||
            fabs(vecRollback[q] - vecA[q]) > rtol * (1.0 + maxRef)) {
            this->PrintError("Time-caching failed after adaptive rollback at the same endpoint.");
            success = false;
            break;
        }
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

bool MollerBoltzmannOperator::CheckBO_MatrixAssemblyAndFluidIdentity() {
    using DREAM::FVM::Grid;
    const len_t nr = 2;
    std::unique_ptr<Grid> fluid(InitializeFluidGrid(nr, 1.0));
    std::unique_ptr<Grid> target(InitializeGridRCylPXi(nr, 12, 20, 1.0, 0.0, 3.0));
    std::unique_ptr<Grid> primary(InitializeGridRCylPXi(nr, 12, 20, 1.0, 0.0, 3.0));
    len_t id_ntot, id_E, id_f;
    std::unique_ptr<DREAM::FVM::UnknownQuantityHandler> uqh(
        BuildUQH_Minimal(fluid.get(), primary.get(), id_ntot, id_E, id_f)
    );
    DREAM::MollerEnergyKernel energy(target.get(), primary.get(), 0.5);
    DREAM::MollerDeltaAngleKernel angle(target.get(), primary.get(), 0.5, 80, 8);
    DREAM::MollerBoltzmannOperator kinetic(target.get(), primary.get(), uqh.get(),
        id_f, &energy, &angle, -1.0);
    DREAM::MollerProductionRateTerm fluidRate(fluid.get(), target.get(), primary.get(),
        uqh.get(), id_f, &energy, -1.0);
    std::vector<real_t> f(primary->GetNCells(), 1e20);
    SetPreviousUnknownData(uqh.get(), id_f, f.data(), 0.0);
    kinetic.Rebuild(0.5, 0.5, uqh.get());
    fluidRate.Rebuild(0.5, 0.5, uqh.get());

    const len_t nk = target->GetNCells();
    real_t ntot[] = {2.0, 3.0};
    std::vector<real_t> kvec(nk, 0), fvec(nr, 0), krhs(nk, 0), frhs(nr, 0);
    kinetic.SetVectorElements(kvec.data(), ntot);
    fluidRate.SetVectorElements(fvec.data(), ntot);
    // Embed the nk-by-nr density block in a square test matrix, as in other
    // coupling tests. The core constructor currently seeds one diagonal
    // entry per row, which is invalid for a tall rectangular matrix.
    // Reserve room for both that diagonal and the density coupling.
    DREAM::FVM::Matrix kmat(nk, nk, 2), fmat(nr, nr, 1);
    kmat.Zero(); fmat.Zero();
    kmat.PartialAssemble(); fmat.PartialAssemble();
    kinetic.SetMatrixElements(&kmat, krhs.data());
    fluidRate.SetMatrixElements(&fmat, frhs.data());
    kmat.Assemble(); fmat.Assemble();

    // Nonzero padding also exercises the requirement that only the first nr
    // columns couple to n_tot; unintended entries must affect this product.
    std::vector<real_t> matrixInput(nk, -7.0);
    std::copy(ntot, ntot+nr, matrixInput.begin());
    std::unique_ptr<real_t[]> matrixProduct(kmat.Multiply(nk, matrixInput.data()));

    const auto agrees = [](real_t a, real_t b) {
        return std::isfinite(a) && std::isfinite(b) &&
            fabs(a-b) <= 1e-12 * std::max({real_t(1e-100), fabs(a), fabs(b)});
    };
    bool success = true;
    len_t offset = 0;
    for (len_t ir=0; ir<nr; ++ir) {
        auto *mg = target->GetMomentumGrid(ir);
        real_t integral = 0;
        for (len_t j=0; j<mg->GetNp2(); ++j) {
            for (len_t i=0; i<mg->GetNp1(); ++i) {
                const len_t q = offset + j*mg->GetNp1() + i;
                success &= agrees(kmat.GetElement(q, ir)*ntot[ir] + krhs[q], kvec[q]);
                success &= agrees(matrixProduct[q] + krhs[q], kvec[q]);
                integral += kvec[q] * target->GetVp(ir, i, j) *
                    mg->GetDp1(i)*mg->GetDp2(j) / target->GetVpVol(ir);
            }
        }
        success &= agrees(fmat.GetElement(ir, ir)*ntot[ir] + frhs[ir], fvec[ir]);
        success &= agrees(integral, fvec[ir]);
        offset += mg->GetNCells();
    }

    // The fluid reaction must also follow an adaptive rollback at fixed t.
    for (auto &x : f) x *= 2;
    SetPreviousUnknownData(uqh.get(), id_f, f.data(), 0.5);
    fluidRate.Rebuild(1.0, 0.5, uqh.get());
    for (auto &x : f) x /= 2;
    SetPreviousUnknownData(uqh.get(), id_f, f.data(), 0.0);
    fluidRate.Rebuild(1.0, 1.0, uqh.get());
    std::vector<real_t> rolledBack(nr, 0);
    fluidRate.SetVectorElements(rolledBack.data(), ntot);
    for (len_t ir=0; ir<nr; ++ir)
        success &= agrees(rolledBack[ir], fvec[ir]);
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
    const len_t poloidalOrder = 8;

    auto *energy = new DREAM::MollerEnergyKernel(gridK, gridP, pCutoff);
    auto *angle = new DREAM::MollerDeltaAngleKernel(gridK, gridP, pCutoff, nXiStars, poloidalOrder);
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
    const len_t poloidalOrder = 8;

    auto *energy = new DREAM::MollerEnergyKernel(gridK, gridP, pCutoff);
    auto *angle = new DREAM::MollerDeltaAngleKernel(gridK, gridP, pCutoff, nXiStars, poloidalOrder);
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

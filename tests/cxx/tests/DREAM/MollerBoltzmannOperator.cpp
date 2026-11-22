#include "MollerBoltzmannOperator.hpp"

#include <algorithm>
#include <cmath>
#include <vector>
#include <memory>

#include "DREAM/Equations/Kinetic/MollerBoltzmannOperator.hpp"
#include "DREAM/Constants.hpp"
#include "DREAM/MollerKernelHandler.hpp"
#include "DREAM/Equations/Fluid/MollerProductionRateTerm.hpp"
#include "DREAM/Equations/Kinetic/MollerDeltaAngleKernel.hpp"
#include "DREAM/Equations/Kinetic/MollerEnergyKernel.hpp"
#include "DREAM/Equations/RosenbluthPutvinskiBoltzmann.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "KnockOnTestHelpers.hpp"
#include "KnockOnTestGrid.hpp"

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

struct MollerBoltzmannOperator::Fixture {
    KnockOnTest::OwnedGrid fluid, target, primary;
    DREAM::MollerEnergyKernel energy;
    DREAM::MollerDeltaAngleKernel angle;
    static constexpr real_t cutoff = 0.2;

    Fixture(DREAM::FVM::Grid *f, DREAM::FVM::Grid *k, DREAM::FVM::Grid *p)
        : fluid(f), target(k), primary(p), energy(k,p,cutoff),
          angle(k,p,cutoff,3,8,4) {}
};

bool MollerBoltzmannOperator::Run(bool) {
    // This table is shared only as immutable data; every check creates fresh
    // unknowns and equation terms. Unlike grids and two radii expose indexing.
    Fixture fixture(InitializeFluidGrid(2,1),
        InitializeGridRCylPXi(2,3,7,1,0,3),
        InitializeGridRCylPXi(2,4,9,1,0,4));
    bool success = true;

    if (CheckBO_LinearityInFPrimary(fixture))
        this->PrintOK("Assembled contribution is linear in f_primary.");
    else {
        success = false;
        this->PrintError("Test failed: assembled contribution is not linear in f_primary.");
    }

    if (CheckBO_JacobianFiniteDifferenceNt(fixture))
        this->PrintOK("Jacobian w.r.t. n_tot matches finite-difference derivative.");
    else {
        success = false;
        this->PrintError("Test failed: Jacobian w.r.t. n_tot does not match FD.");
    }

    if (CheckBO_GlobalProductionIdentity(fixture))
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

    if (CheckBO_TimeCachingRegression(fixture))
        this->PrintOK("Time caching regression test passed.");
    else {
        success = false;
        this->PrintError("Test failed: time caching regression.");
    }

    if (CheckBO_MatrixAssemblyAndFluidIdentity(fixture))
        this->PrintOK("Kinetic/fluid matrix assembly and production identity passed.");
    else {
        success = false;
        this->PrintError("Test failed: kinetic/fluid matrix assembly or production identity.");
    }

    if (CheckBO_NonNegativity(fixture))
        this->PrintOK("Non-negativity sanity check passed.");
    else {
        success = false;
        this->PrintError("Test failed: non-negativity sanity check.");
    }

    if (CheckBO_RadiusLocality(fixture))
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

    if (CheckBO_ChiuHarveyGain())
        this->PrintOK("Chiu–Harvey gain splits the two primary poles and preserves production.");
    else {
        success = false;
        this->PrintError("Test failed: Chiu–Harvey pole weights or production.");
    }

    if (CheckBO_RosenbluthPutvinskiGain())
        this->PrintOK("Boltzmann RP uses partial primary cells, signed poles, and conservative particle bookkeeping.");
    else {
        success = false;
        this->PrintError("Test failed: Boltzmann RP primary density or source.");
    }

    if (CheckBO_SharedExternalGain())
        this->PrintOK("All Boltzmann models share the external gain and close the grid/fluid birth balance.");
    else {
        success = false;
        this->PrintError("Test failed: shared external gain or particle balance.");
    }

    return success;
}

bool MollerBoltzmannOperator::CheckBO_ChiuHarveyGain() {
    constexpr real_t pCutoff = 0.1;
    constexpr real_t tolerance = 1e-12;
    const std::vector<real_t> faces = {-1, -0.9, -0.6, -0.1, 0, 0.3, 0.6, 0.8, 0.9, 1};
    for (bool toroidal : {false, true}) {
        auto primary = toroidal
            ? KnockOnTest::ToroidalGrid(0.3, faces, 2, 8, 0, 6)
            : KnockOnTest::OwnedGrid(InitializeGridRCylPXi(2, 8, 9, 1, 0, 6));
        auto target = toroidal
            ? KnockOnTest::ToroidalGrid(0.3, KnockOnTest::UniformPitchFaces(20), 2, 7, 0, 3)
            : KnockOnTest::OwnedGrid(InitializeGridRCylPXi(2, 7, 20, 1, 0, 3));
        len_t idN, idE, idF;
        std::unique_ptr<DREAM::FVM::UnknownQuantityHandler> unknowns(
            BuildUQH_Minimal(target.get(), primary.get(), idN, idE, idF));
        DREAM::MollerEnergyKernel energy(target.get(), primary.get(), pCutoff);
        DREAM::MollerDeltaAngleKernel angle(
            target.get(), primary.get(), pCutoff, 12, 8, 4,
            DREAM::OptionConstants::BOLTZMANN_MODEL_CHIU_HARVEY
        );
        DREAM::MollerBoltzmannOperator source(
            target.get(), primary.get(), unknowns.get(), idF, &energy, &angle
        );

        const auto *mgP = primary->GetMomentumGrid(0);
        const auto *mgK = target->GetMomentumGrid(0);
        const len_t np = mgP->GetNp1();
        const len_t nxi = mgP->GetNp2();
        std::vector<real_t> f(primary->GetNCells(), 0);
        std::vector<real_t> positive, negative, trapped, mixed;
        const auto evaluate = [&](len_t pitch, real_t time, std::vector<real_t> &gain) {
            std::fill(f.begin(), f.end(), 0);
            for (len_t ir = 0; ir < primary->GetNr(); ++ir) {
                const len_t offset = ir * mgP->GetNCells() + pitch*np;
                for (len_t k = 0; k < np; ++k)
                    f[offset + k] = 1 / (mgP->GetDp2(pitch) * primary->GetVp(ir, k, pitch));
            }
            SetPreviousUnknownData(unknowns.get(), idF, f.data(), time);
            source.Rebuild(time+1, 1, unknowns.get());
            ApplyTerm(source, target.get(), 1, gain);
        };
        evaluate(nxi-1, 0, positive);
        angle.GridRebuilt();
        source.GridRebuilt();
        evaluate(0, 1, negative);
        // A fully trapped toroidal cell or the cylinder's zero-straddling cell.
        evaluate(toroidal ? 5 : nxi/2, 2, trapped);
        if (toroidal)
            evaluate(6, 3, mixed);

        for (len_t ir = 0; ir < target->GetNr(); ++ir) {
            const len_t base = ir * mgK->GetNCells();
            real_t rate = 0;
            for (len_t k = 0; k < np; ++k)
                rate += mgP->GetDp1(k) * energy.TotalCS(k) / target->GetVpVol(ir);
            if (!(rate > 0))
                return false;
            for (const auto *gain : {&positive, &negative, &trapped, &mixed}) {
                if (gain->empty())
                    continue;
                const real_t integral = target->IntegralMomentumAtRadius(ir, gain->data()+base);
                if (!std::isfinite(integral) || std::abs(integral/rate - 1) > tolerance)
                    return false;
            }

            const auto *rg = primary->GetRadialGrid();
            const real_t xiT = rg->GetXi0TrappedBoundary(ir);
            real_t mixedFraction = 0;
            if (toroidal) {
                // The test grid has a cell [.6,.8] straddling +xiT.
                // Get each subinterval's independent native orbit volume.
                if (!(0.6 < xiT && xiT < 0.8))
                    return false;
                auto *fsa = rg->GetFluxSurfaceAverager();
                const auto volume = [&](real_t lo, real_t hi) {
                    return (hi - lo) * fsa->EvaluateCellAveragedBounceIntegralOverP2(
                        ir, lo, hi, DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION,
                        DREAM::FVM::RadialGrid::BA_FUNC_UNITY, nullptr, DREAM::FVM::RadialGrid::BA_PARAM_UNITY
                    );
                };
                const real_t trappedVolume = volume(0.6, xiT);
                const real_t passingVolume = volume(xiT, 0.8);
                mixedFraction = 0.5 + 0.5 * passingVolume/(passingVolume+trappedVolume);
            }
            for (len_t j = 0; j < mgK->GetNp2(); ++j) {
                for (len_t i = 0; i < mgK->GetNp1(); ++i) {
                    const len_t q = base + j * mgK->GetNp1() + i;
                    const real_t a = positive[q], b = negative[q];
                    const real_t scale = std::max(a, b);
                    if (!(std::isfinite(a) && std::isfinite(b) && a >= 0 && b >= 0))
                        return false;
                    if (std::abs(trapped[q] - 0.5*(a+b)) > tolerance*scale)
                        return false;
                    if (mgK->GetP2_f(j+1) <= -xiT && a != 0)
                        return false;
                    if (mgK->GetP2_f(j) >= xiT && b != 0)
                        return false;
                    if (toroidal) {
                        if (std::abs(mixed[q] - (mixedFraction*a + (1-mixedFraction)*b)) > tolerance*scale)
                            return false;
                        if (target->GetVp(ir, i, j) == 0 && (a != 0 || b != 0))
                            return false;
                    } else {
                        const len_t mirror = base + (mgK->GetNp2()-1-j) * mgK->GetNp1() + i;
                        if (std::abs(a-negative[mirror]) > tolerance*std::max(a, negative[mirror]))
                            return false;
                    }
                }
            }
        }
    }
    return true;
}

bool MollerBoltzmannOperator::CheckBO_RosenbluthPutvinskiGain() {
    using RP = DREAM::RosenbluthPutvinskiBoltzmann;
    using RPTerm = DREAM::RosenbluthPutvinskiBoltzmannTerm;
    auto fluid = KnockOnTest::OwnedGrid(InitializeFluidGrid(1, 1));
    auto hot = KnockOnTest::OwnedGrid(InitializeGridRCylPXi(1, 4, 10, 1, 0, 2));
    auto runaway = KnockOnTest::OwnedGrid(InitializeGridRCylPXi(1, 4, 12, 1, 2, 6));
    constexpr real_t pCutoff = 0.3;
    constexpr real_t nTot[] = {2e19};
    const auto close = [](real_t a, real_t b) {
        return std::isfinite(a) && std::isfinite(b)
            && std::abs(a-b) <= 1e-12 * std::max(std::abs(a), std::abs(b));
    };
    // RP spectrum independently integrated in gamma, with finite and
    // infinite upper boundaries.
    const auto integrated = [](real_t lower, real_t upper) {
        const auto invEnergy = [](real_t p) {
            return std::isinf(p) ? 0 : (std::sqrt(1+p*p)+1)/(p*p);
        };
        return 2 * M_PI * DREAM::Constants::r0 * DREAM::Constants::r0 * DREAM::Constants::c
            * (invEnergy(lower) - invEnergy(upper));
    };
    for (bool twoGrids : {false, true}) {
        const real_t pRe = twoGrids ? 2.25 : 0.625;
        const real_t pMax = twoGrids ? 6 : 2;
        const real_t hotDensity = 100, reDensity = 40, externalDensity = 30;
        const real_t nRe[] = {externalDensity + (twoGrids ? reDensity : 0)};
        const real_t electric[] = {-1};
        len_t idN, idE, idF;
        std::unique_ptr<DREAM::FVM::UnknownQuantityHandler> u(
            BuildUQH_Minimal(fluid.get(), hot.get(), idN, idE, idF));
        const len_t idNRe = u->InsertUnknown(DREAM::OptionConstants::UQTY_N_RE, "0", fluid.get());
        u->SetInitialValue(idNRe, nRe);
        u->SetInitialValue(idE, electric);
        const auto beam = [](const DREAM::FVM::Grid *grid, len_t j, len_t k, real_t density) {
            const auto *mg = grid->GetMomentumGrid(0);
            const len_t index = j*mg->GetNp1() + k;
            std::vector<real_t> f(grid->GetNCells(), 0);
            f[index] = density * grid->GetVpVol(0)
                / (grid->GetVp(0)[index] * mg->GetDp1(k) * mg->GetDp2(j));
            return f;
        };
        auto fHot = beam(hot.get(), 9, 1, hotDensity);
        // Also seed a below-threshold hot cell to ensure it never contributes.
        fHot[0] = 1e8;
        u->SetInitialValue(idF, fHot.data());
        if (twoGrids) {
            const len_t idFRe = u->InsertUnknown(DREAM::OptionConstants::UQTY_F_RE, "0", runaway.get());
            auto fRe = beam(runaway.get(), 0, 0, reDensity);
            u->SetInitialValue(idFRe, fRe.data());
        }
        RP source(hot.get(), twoGrids ? runaway.get() : nullptr, u.get(), pCutoff, pRe, 4, true);
        RPTerm hotTerm(hot.get(), u.get(), &source, RP::HOTTAIL, 1);
        RPTerm reDensityTerm(fluid.get(), u.get(), &source, RP::RUNAWAY_DENSITY, 1);
        RPTerm totalTerm(fluid.get(), u.get(), &source, RP::ALL, 1);
        hotTerm.Rebuild(1, 1, u.get());

        const real_t nResolved = 0.75 * (twoGrids ? reDensity : hotDensity);
        if (!close(source.GetResolvedDensity(!twoGrids)[0], nResolved)
            || source.GetResolvedDensity(twoGrids)[0] != 0
            || !close(source.GetExternalDensity()[0], externalDensity))
            return false;

        std::vector<real_t> fromResolved(hot->GetNCells()), fromExternal(hot->GetNCells());
        hotTerm.SetResolvedSource(fromResolved.data(), nTot);
        hotTerm.SetExternalSource(fromExternal.data(), nTot);
        const real_t onGrid = integrated(pCutoff, 2);
        if (!close(hot->IntegralMomentumAtRadius(0, fromResolved.data()), nTot[0]*nResolved*onGrid)
            || !close(hot->IntegralMomentumAtRadius(0, fromExternal.data()), nTot[0]*externalDensity*onGrid))
            return false;
        const auto *mg = hot->GetMomentumGrid(0);
        for (len_t j = 0; j < mg->GetNp2(); ++j) {
            for (len_t i = 0; i < mg->GetNp1(); ++i) {
                const len_t index = j*mg->GetNp1()+i;
                // External primaries stay at the positive pole even for E<0.
                if (j < mg->GetNp2()/2 && fromExternal[index] != 0)
                    return false;
                if ((j >= mg->GetNp2()/2) == twoGrids && fromResolved[index] != 0)
                    return false;
            }
        }

        real_t nPrimary = 0, reBirths = 0, allBirths = 0;
        hotTerm.GetPrimaryDensity(&nPrimary);
        reDensityTerm.SetVectorElements(&reBirths, nTot);
        totalTerm.SetVectorElements(&allBirths, nTot);
        const real_t coefficient = nTot[0] * (nResolved + externalDensity);
        const real_t infinity = std::numeric_limits<real_t>::infinity();
        if (!close(nPrimary, nResolved + externalDensity)
            || !close(reBirths, coefficient * integrated(2, infinity))
            || !close(allBirths, coefficient * integrated(pCutoff, infinity)))
            return false;
        if (twoGrids) {
            RPTerm reTerm(runaway.get(), u.get(), &source, RP::RUNAWAY, 1);
            std::vector<real_t> gain(runaway->GetNCells());
            reTerm.SetVectorElements(gain.data(), nTot);
            if (!close(runaway->IntegralMomentumAtRadius(0, gain.data()), coefficient * integrated(2, pMax)))
                return false;
        }
        // A current-step change in n_re must not masquerade as a new
        // external population while the kinetic moment is still lagged.
        const real_t currentNRe[] = {nRe[0] + 10};
        u->Store(idNRe, currentNRe);
        source.GridRebuilt();
        hotTerm.Rebuild(1, 1, u.get());
        hotTerm.GetPrimaryDensity(&nPrimary);
        if (!close(nPrimary, nResolved + externalDensity))
            return false;

        if (twoGrids) {
            const real_t lowNRe[] = {reDensity - 1};
            SetPreviousUnknownData(u.get(), idNRe, lowNRe, 1);
            SetPreviousUnknownData(u.get(), idF, fHot.data(), 1);
            hotTerm.Rebuild(2, 1, u.get());
            if (source.GetExternalDensity()[0] != 0 || !close(source.GetRawExternalDensity()[0], -1))
                return false;
        }
    }
    return true;
}

/**
 * Check the common infinite-energy component independently of the resolved
 * model. Cover both geometries and all supported kinetic-grid layouts.
 */
bool MollerBoltzmannOperator::CheckBO_SharedExternalGain() {
    using RP = DREAM::RosenbluthPutvinskiBoltzmann;
    using RPTerm = DREAM::RosenbluthPutvinskiBoltzmannTerm;
    using OC = DREAM::OptionConstants;
    constexpr real_t pCutoff = 0.3;
    constexpr real_t nTot[] = {2e19, 3e19};
    constexpr real_t rawExternal[] = {7, -0.25};
    const auto close = [](real_t a, real_t b) {
        return std::isfinite(a) && std::isfinite(b)
            && std::abs(a-b) <= 2e-12 * std::max({real_t(1), std::abs(a), std::abs(b)});
    };

    for (bool toroidal : {false, true}) {
        auto fluid = KnockOnTest::OwnedGrid(InitializeFluidGrid(2, 1));
        const auto faces = KnockOnTest::UniformPitchFaces(10);
        auto hot = toroidal ? KnockOnTest::ToroidalGrid(0.3, faces, 2, 4, 0, 2)
            : KnockOnTest::OwnedGrid(InitializeGridRCylPXi(2, 4, 10, 1, 0, 2));
        auto runaway = toroidal ? KnockOnTest::ToroidalGrid(0.3, faces, 2, 4, 2, 6)
            : KnockOnTest::OwnedGrid(InitializeGridRCylPXi(2, 4, 10, 1, 2, 6));

        // Hot only, runaway only, and two adjoining kinetic grids.
        for (len_t layout = 0; layout < 3; ++layout) {
            auto *h = layout == 1 ? nullptr : hot.get();
            auto *r = layout == 0 ? nullptr : runaway.get();
            DREAM::FVM::UnknownQuantityHandler unknowns;
            const len_t idNRe = unknowns.InsertUnknown(OC::UQTY_N_RE, "0", fluid.get());
            std::vector<real_t> fHot(hot->GetNCells(), 2), fRe(runaway->GetNCells(), 3);
            len_t idFHot = 0, idFRe = 0;
            if (h) {
                idFHot = unknowns.InsertUnknown(OC::UQTY_F_HOT, "0", h);
                unknowns.SetInitialValue(idFHot, fHot.data());
            }
            if (r) {
                idFRe = unknowns.InsertUnknown(OC::UQTY_F_RE, "0", r);
                unknowns.SetInitialValue(idFRe, fRe.data());
            }
            real_t nRe[2];
            for (len_t ir = 0; ir < 2; ++ir) {
                nRe[ir] = rawExternal[ir];
                if (r)
                    nRe[ir] += r->IntegralMomentumAtRadius(ir, fRe.data() + ir*r->GetMomentumGrid(0)->GetNCells());
            }
            unknowns.SetInitialValue(idNRe, nRe);
            std::array<std::vector<real_t>, 2> referenceExternal;

            for (auto model : {OC::BOLTZMANN_MODEL_FULL_MOLLER, OC::BOLTZMANN_MODEL_CHIU_HARVEY,
                               OC::BOLTZMANN_MODEL_ROSENBLUTH_PUTVINSKI}) {
                const bool resolvedRP = model == OC::BOLTZMANN_MODEL_ROSENBLUTH_PUTVINSKI;
                DREAM::MollerKernelHandler kernels(
                    h, r, &unknowns, pCutoff, 3, 4, 2, model, resolvedRP ? 0.625 : 0
                );
                auto *source = kernels.GetRP();
                RPTerm all(fluid.get(), &unknowns, source, RP::ALL, 1);
                RPTerm density(fluid.get(), &unknowns, source, RP::RUNAWAY_DENSITY, 1);
                all.Rebuild(1, 1, &unknowns);
                real_t raw[2], external[2], primary[2], above[2], allBirths[2] = {}, reBirths[2] = {};
                all.GetRawExternalDensity(raw);
                all.GetExternalDensity(external);
                all.GetPrimaryDensity(primary);
                all.GetAboveGridProduction(above, nTot);
                all.SetVectorElements(allBirths, nTot);
                density.SetVectorElements(reBirths, nTot);
                for (len_t ir = 0; ir < 2; ++ir) {
                    if (!close(raw[ir], rawExternal[ir])
                        || !close(external[ir], std::max(real_t(0), rawExternal[ir]))
                        || (!resolvedRP && !close(primary[ir], external[ir])))
                        return false;
                }

                real_t kineticBirths[2] = {}, runawayBirths[2] = {};
                for (len_t target = 0; target < 2; ++target) {
                    auto *grid = target == 0 ? h : r;
                    if (!grid)
                        continue;
                    RPTerm term(grid, &unknowns, source, target == 0 ? RP::HOTTAIL : RP::RUNAWAY, 1);
                    std::vector<real_t> gain(grid->GetNCells()), resolved(gain.size()), offGrid(gain.size());
                    term.SetVectorElements(gain.data(), nTot);
                    term.SetResolvedSource(resolved.data(), nTot);
                    term.SetExternalSource(offGrid.data(), nTot);
                    if (model == OC::BOLTZMANN_MODEL_FULL_MOLLER)
                        referenceExternal[target] = offGrid;
                    for (len_t cell = 0; cell < gain.size(); ++cell) {
                        if (!close(gain[cell], resolved[cell] + offGrid[cell])
                            || !close(offGrid[cell], referenceExternal[target][cell])
                            || (!resolvedRP && resolved[cell] != 0))
                            return false;
                    }
                    for (len_t ir = 0, offset = 0; ir < 2; ++ir) {
                        const real_t integral = grid->IntegralMomentumAtRadius(ir, gain.data() + offset);
                        kineticBirths[ir] += integral;
                        if (target == 1)
                            runawayBirths[ir] = integral;
                        // The negative raw inventory at radius 1 contributes no gain.
                        if (ir == 1 && grid->IntegralMomentumAtRadius(ir, offGrid.data() + offset) != 0)
                            return false;
                        offset += grid->GetMomentumGrid(ir)->GetNCells();
                    }
                }
                for (len_t ir = 0; ir < 2; ++ir) {
                    if (!close(kineticBirths[ir] + above[ir], allBirths[ir])
                        || !close(runawayBirths[ir] + above[ir], reBirths[ir]))
                        return false;
                }

                // Newton iterates must not mix current n_re/f_re with lagged
                // primaries, even when a grid rebuild invalidates the cache.
                const real_t currentNRe[] = {nRe[0]+100, nRe[1]+100};
                unknowns.Store(idNRe, currentNRe);
                if (r) {
                    std::vector<real_t> currentF(fRe.size(), 30);
                    unknowns.Store(idFRe, currentF.data());
                }
                kernels.GridRebuilt();
                all.Rebuild(1, 1, &unknowns);
                all.GetRawExternalDensity(raw);
                for (len_t ir = 0; ir < 2; ++ir)
                    if (!close(raw[ir], rawExternal[ir]))
                        return false;
                unknowns.Store(idNRe, nRe);
                if (r)
                    unknowns.Store(idFRe, fRe.data());
            }
        }
    }
    return true;
}

bool MollerBoltzmannOperator::CheckBO_IsotropicZeroCurrentAndReflection() {
    KnockOnTest::OwnedGrid fluid(InitializeFluidGrid(1, 1.0));
    KnockOnTest::OwnedGrid target(InitializeGridRCylPXi(1, 5, 15, 1.0, 0, 2));
    KnockOnTest::OwnedGrid primary(InitializeGridRCylPXi(1, 6, 12, 1.0, 0, 3));
    len_t idN, idE, idF;
    std::unique_ptr<DREAM::FVM::UnknownQuantityHandler> uqh(
        BuildUQH_Minimal(fluid.get(), primary.get(), idN, idE, idF));
    DREAM::MollerEnergyKernel energy(target.get(), primary.get(), 0.1);
    DREAM::MollerDeltaAngleKernel angle(target.get(), primary.get(), 0.1,
        3, 8);
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

bool MollerBoltzmannOperator::CheckBO_LinearityInFPrimary(const Fixture& fixture) {
    const real_t tol = 1e-10;

    auto *gridF = fixture.fluid.get(), *gridK = fixture.target.get(), *gridP = fixture.primary.get();
    len_t id_ntot, id_E, id_f;
    std::unique_ptr<DREAM::FVM::UnknownQuantityHandler> unknowns(
        BuildUQH_Minimal(gridF,gridP,id_ntot,id_E,id_f));
    auto *uqh = unknowns.get();
    const real_t scaleFactor = 1;
    const auto *energy = &fixture.energy;
    const auto *angle = &fixture.angle;
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

    // Explicitly test the zero input as well as additivity.
    std::fill(fAB,fAB+NP,0);
    SetPreviousUnknownData(uqh,id_f,fAB,3.0);
    op.Rebuild(4,1,uqh);
    std::vector<real_t> zero;
    ApplyTerm(op,gridK,ntot,zero);
    for (real_t value : zero) success &= std::isfinite(value) && value==0;
    success &= std::any_of(RA.begin(),RA.end(),[](real_t value) { return value>0; });

    delete[] fA;
    delete[] fB;
    delete[] fAB;

    return success;
}

bool MollerBoltzmannOperator::CheckBO_JacobianFiniteDifferenceNt(const Fixture& fixture) {
    const real_t tol = 1e-6;
    const real_t eps = 1e-7;

    const len_t nr = fixture.target->GetNr();
    auto *gridF = fixture.fluid.get(), *gridK = fixture.target.get(), *gridP = fixture.primary.get();
    len_t id_ntot, id_E, id_f;
    std::unique_ptr<DREAM::FVM::UnknownQuantityHandler> unknowns(
        BuildUQH_Minimal(gridF,gridP,id_ntot,id_E,id_f));
    auto *uqh = unknowns.get();
    const real_t scaleFactor = 1;
    const auto *energy = &fixture.energy;
    const auto *angle = &fixture.angle;
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
    DREAM::FVM::Matrix jac(NK, NK, 2);
    jac.Zero();
    jac.PartialAssemble();
    if (!op.SetJacobianBlock(id_ntot,id_ntot,&jac,nullptr)) return false;
    if (op.SetJacobianBlock(id_ntot,id_f,&jac,nullptr)) return false;
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

    return success;
}

bool MollerBoltzmannOperator::CheckBO_GlobalProductionIdentity(const Fixture& fixture) {
    const real_t rtol = 1e-14;
    const real_t atol = 1e-14;

    auto *gridF = fixture.fluid.get(), *gridK = fixture.target.get(), *gridP = fixture.primary.get();
    len_t id_ntot, id_E, id_f;
    std::unique_ptr<DREAM::FVM::UnknownQuantityHandler> unknowns(
        BuildUQH_Minimal(gridF,gridP,id_ntot,id_E,id_f));
    auto *uqh = unknowns.get();
    const real_t pCutoff = Fixture::cutoff, scaleFactor = 1;
    const auto *energy = &fixture.energy;
    const auto *angle = &fixture.angle;
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
    const bool success = RHS>0 && std::isfinite(LHS) && (diff <= atol + rtol * (fabs(LHS) + fabs(RHS)));

    if (!success) {
        this->PrintError("Global identity failed:");
        this->PrintError("  LHS: %.8g", LHS);
        this->PrintError("  RHS: %.8g", RHS);
        this->PrintError("  diff: %.8g", diff);
    }

    delete[] f;

    return success;
}

bool MollerBoltzmannOperator::CheckBO_HotRunawayGlobalProductionIdentity() {
    const real_t rtol = 1e-14;
    const real_t atol = 1e-14;
    const real_t absTolNonNeg = 1e-13;

    const len_t nr = 2;
    KnockOnTest::OwnedGrid fluid(InitializeFluidGrid(nr,1));
    auto hot = KnockOnTest::ToroidalGrid(0.3,KnockOnTest::UniformPitchFaces(7),nr,3,0,3);
    auto runaway = KnockOnTest::ToroidalGrid(0.3,KnockOnTest::UniformPitchFaces(9),nr,4,0,6);
    auto *gridHot = hot.get(), *gridRe = runaway.get();
    len_t id_ntot, id_E, id_f_re;
    std::unique_ptr<DREAM::FVM::UnknownQuantityHandler> unknowns(
        BuildUQH_MinimalFRe(fluid.get(),gridRe,id_ntot,id_E,id_f_re));
    auto *uqh = unknowns.get();
    const real_t pCutoff = 0.2, ntot = 1, scaleFactor = 1;
    DREAM::MollerEnergyKernel energy(gridHot,gridRe,pCutoff);
    DREAM::MollerDeltaAngleKernel angle(gridHot,gridRe,pCutoff,3,8,4);
    DREAM::MollerBoltzmannOperator op(gridHot,gridRe,uqh,id_f_re,&energy,&angle,scaleFactor);

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

    len_t offset = 0;
    for (len_t ir=0; ir<nr; ++ir) {
        const auto *mg = gridHot->GetMomentumGrid(ir);
        for (len_t j=0; j<mg->GetNp2(); ++j)
            for (len_t i=0; i<mg->GetNp1(); ++i)
                if (gridHot->GetVp(ir,i,j)==0 && vec[offset+j*mg->GetNp1()+i]!=0)
                    success = false;
        offset += mg->GetNCells();
    }
    const real_t LHS = IntegrateTotalProductionOverKnockonGrid(gridHot, vec.data());

    const real_t RHS0 = PredictTotalProductionFromMollerS(gridHot, gridRe, f, pCutoff);
    const real_t RHS = scaleFactor * ntot * RHS0;

    success &= RHS>0 && std::isfinite(LHS);
    const real_t diff = fabs(LHS - RHS);
    if (diff > atol + rtol * (fabs(LHS) + fabs(RHS))) {
        this->PrintError("hot<-re global identity failed:");
        this->PrintError("  LHS: %.8g", LHS);
        this->PrintError("  RHS: %.8g", RHS);
        this->PrintError("  diff: %.8g", diff);
        success = false;
    }

    delete[] f;

    return success;
}

bool MollerBoltzmannOperator::CheckBO_TimeCachingRegression(const Fixture& fixture) {
    const real_t rtol = 1e-12;

    auto *gridF = fixture.fluid.get(), *gridK = fixture.target.get(), *gridP = fixture.primary.get();
    len_t id_ntot, id_E, id_f;
    std::unique_ptr<DREAM::FVM::UnknownQuantityHandler> unknowns(
        BuildUQH_Minimal(gridF,gridP,id_ntot,id_E,id_f));
    auto *uqh = unknowns.get();
    const real_t scaleFactor = 1;
    const auto *energy = &fixture.energy;
    const auto *angle = &fixture.angle;
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

    bool success = maxRef>0 && std::isfinite(maxRef);
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

    return success;
}

bool MollerBoltzmannOperator::CheckBO_MatrixAssemblyAndFluidIdentity(const Fixture& fixture) {
    using DREAM::FVM::Grid;
    const len_t nr = fixture.target->GetNr();
    auto *fluid = fixture.fluid.get(), *target = fixture.target.get(), *primary = fixture.primary.get();
    len_t id_ntot, id_E, id_f;
    std::unique_ptr<DREAM::FVM::UnknownQuantityHandler> uqh(
        BuildUQH_Minimal(fluid,primary,id_ntot,id_E,id_f));
    const auto& energy = fixture.energy;
    const auto& angle = fixture.angle;
    DREAM::MollerBoltzmannOperator kinetic(target, primary, uqh.get(),
        id_f, &energy, &angle, -1.0);
    DREAM::MollerProductionRateTerm fluidRate(fluid, target, primary,
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

bool MollerBoltzmannOperator::CheckBO_NonNegativity(const Fixture& fixture) {
    const real_t absTol = 1e-13;

    auto *gridF = fixture.fluid.get(), *gridK = fixture.target.get(), *gridP = fixture.primary.get();
    len_t id_ntot, id_E, id_f;
    std::unique_ptr<DREAM::FVM::UnknownQuantityHandler> unknowns(
        BuildUQH_Minimal(gridF,gridP,id_ntot,id_E,id_f));
    auto *uqh = unknowns.get();
    const real_t scaleFactor = 1;
    const auto *energy = &fixture.energy;
    const auto *angle = &fixture.angle;
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

    return success;
}

bool MollerBoltzmannOperator::CheckBO_RadiusLocality(const Fixture& fixture) {
    const real_t absTol = 1e-14;
    const real_t relTol = 1e-12;

    const len_t nr = fixture.target->GetNr();
    auto *gridF = fixture.fluid.get(), *gridK = fixture.target.get(), *gridP = fixture.primary.get();
    len_t id_ntot, id_E, id_f;
    std::unique_ptr<DREAM::FVM::UnknownQuantityHandler> unknowns(
        BuildUQH_Minimal(gridF,gridP,id_ntot,id_E,id_f));
    auto *uqh = unknowns.get();
    const real_t scaleFactor = 1;
    const auto *energy = &fixture.energy;
    const auto *angle = &fixture.angle;
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

    bool success = std::any_of(vecAll.begin(),vecAll.end(),[](real_t value) { return value>0; });
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

    return success;
}

#include "DREAM/Equations/RosenbluthPutvinskiBoltzmann.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "DREAM/DREAMException.hpp"
#include "DREAM/Equations/KnockOnUtilities.hpp"
#include "DREAM/Settings/OptionConstants.hpp"

using namespace DREAM;

namespace {
namespace K = KnockOnUtilities::Kinematics;
namespace Q = KnockOnUtilities::Quadrature;

constexpr real_t infiniteMomentum = std::numeric_limits<real_t>::infinity();

/** Infinite-energy v*sigma integrated over a momentum interval. */
real_t IntegrateRPSpectrum(real_t lower, real_t upper) {
    if (upper <= lower)
        return 0;
    const real_t atUpper = std::isinf(upper) ? 0 : K::EvaluateMollerFlux(upper, infiniteMomentum);
    return atUpper - K::EvaluateMollerFlux(lower, infiniteMomentum);
}
}

/**
 * Prepare the infinite-energy pole shapes and, in RP mode, primary weights.
 * The highest kinetic grid sets the start of the off-grid population.
 */
RosenbluthPutvinskiBoltzmann::RosenbluthPutvinskiBoltzmann(
    const FVM::Grid *hottail, const FVM::Grid *runaway, FVM::UnknownQuantityHandler *unknowns,
    real_t pCutoff, real_t pRe, len_t poloidalOrder, bool includeResolvedPrimaries
) : hottail(hottail), runaway(runaway), unknowns(unknowns),
    idNRe(unknowns->GetUnknownID(OptionConstants::UQTY_N_RE)),
    pCutoff(pCutoff), pRe(pRe), poloidalOrder(poloidalOrder),
    includeResolvedPrimaries(includeResolvedPrimaries) {
    if ((!hottail && !runaway) || poloidalOrder < 1
        || !std::isfinite(pCutoff) || pCutoff <= 5*std::numeric_limits<real_t>::epsilon())
        throw DREAMException("Boltzmann RP spectrum requires a kinetic grid, poloidal order >= 1 and a finite positive outgoing cutoff.");
    if (includeResolvedPrimaries && (!std::isfinite(pRe) || pRe <= 0))
        throw DREAMException("Resolved Boltzmann RP primaries require a finite positive momentum threshold.");

    if (hottail && runaway && hottail->GetNr() != runaway->GetNr())
        throw DREAMException("Boltzmann RP requires the kinetic grids to share the same radii.");

    if (runaway)
        idFRe = unknowns->GetUnknownID(OptionConstants::UQTY_F_RE);

    // Full Moller and CH already integrate their resolved primaries. Avoid
    // building threshold/pole weights or repeating that moment for them.
    if (includeResolvedPrimaries) {
        if (hottail)
            primaryGrids.push_back({hottail, unknowns->GetUnknownID(OptionConstants::UQTY_F_HOT), {}, {}});
        if (runaway)
            primaryGrids.push_back({runaway, idFRe, {}, {}});
    }
    GridRebuilt();
}

void RosenbluthPutvinskiBoltzmann::GridRebuilt() {
    const len_t nr = (runaway ? runaway : hottail)->GetNr();
    for (auto &pole : resolved) pole.assign(nr, 0);
    external.resize(nr);
    externalRaw.resize(nr);
    densitiesValid = false;

    for (auto &primary : primaryGrids) {
        KnockOnUtilities::BuildPrimaryPoleFractions(primary.grid, primary.positiveFraction);
        primary.momentumWidth.clear();
        for (len_t ir = 0; ir < nr; ++ir) {
            const auto *mg = primary.grid->GetMomentumGrid(ir);
            for (len_t k = 0; k < mg->GetNp1(); ++k) {
                // The stored Vprime and f are constant in this finite-volume
                // moment. A cut cell contributes only its admitted width.
                const real_t lower = std::max(pRe, mg->GetP1_f(k));
                primary.momentumWidth.push_back(std::max(real_t(0), mg->GetP1_f(k+1) - lower));
            }
        }
    }

    if (hottail) BuildSource(HOTTAIL, hottail);
    if (runaway) BuildSource(RUNAWAY, runaway);
    const auto *outer = (runaway ? runaway : hottail)->GetMomentumGrid(0);
    const real_t pMax = outer->GetP1_f(outer->GetNp1());
    production[ABOVE_GRID] = IntegrateRPSpectrum(std::max(pCutoff, pMax), infiniteMomentum);
    production[RUNAWAY_DENSITY] = production[RUNAWAY] + production[ABOVE_GRID];
    production[ALL] = production[HOTTAIL] + production[RUNAWAY_DENSITY];
}

void RosenbluthPutvinskiBoltzmann::BuildSource(Region region, const FVM::Grid *grid) {
    const len_t nPoles = includeResolvedPrimaries ? 2 : 1;
    for (len_t pole = 0; pole < nPoles; ++pole)
        source[region][pole].resize(grid->GetNCells());

    const auto *mg = grid->GetMomentumGrid(0);
    const len_t np = mg->GetNp1(), nxi = mg->GetNp2();
    std::vector<real_t> spectrum(np), cosine(np), column(nxi);
    production[region] = 0;
    for (len_t i = 0; i < np; ++i) {
        const real_t upper = mg->GetP1_f(i+1);
        const real_t lower = std::max(pCutoff, mg->GetP1_f(i));
        const real_t integral = IntegrateRPSpectrum(lower, upper);
        production[region] += integral;
        spectrum[i] = integral / mg->GetDp1(i);
        if (integral > 0)
            cosine[i] = K::EvaluateXiStar(std::max(lower, mg->GetP1(i)), infiniteMomentum);
    }

    const Q::GaussLegendreRule rule(poloidalOrder);
    for (len_t ir = 0; ir < grid->GetNr(); ++ir) {
        Q::OrbitGeometry geometry(grid->GetRadialGrid()->GetFluxSurfaceAverager(), ir);
        const real_t *volume = grid->GetVp(ir);
        for (len_t i = 0; i < np; ++i) {
            for (len_t pole = 0; pole < nPoles; ++pole) {
                if (spectrum[i] > 0) {
                    // One scattering cosine per outgoing p cell in the
                    // infinite-energy limit. Tabulate directly on that grid.
                    KnockOnUtilities::SetFieldAlignedDeltaColumnOnGrid(
                        ir, cosine[i], pole == 0 ? +1 : -1, grid, column.data(), geometry, rule
                    );
                }
                for (len_t j = 0; j < nxi; ++j) {
                    const len_t cell = j*np + i;
                    source[region][pole][ir*np*nxi + cell] = spectrum[i] > 0 && volume[cell] > 0
                        ? spectrum[i] * column[j] * grid->GetVpVol(ir) / volume[cell] : 0;
                }
            }
        }
    }
}

const real_t *RosenbluthPutvinskiBoltzmann::GetSource(Region region, bool positive) const {
    return source[region][positive ? 0 : 1].data();
}

/**
 * Use one previous-time moment for the kinetic gain and every fluid balance.
 * Remove the full f_re inventory from n_re at the same time level, so
 * the threshold applies only to the resolved primaries.
 */
void RosenbluthPutvinskiBoltzmann::RebuildPrimaryDensities(FVM::UnknownQuantityHandler *u) {
    if (u != unknowns) {
        unknowns = u;
        densitiesValid = false;
    }
    const real_t t = unknowns->GetUnknownDataPreviousTime(idNRe);
    if (densitiesValid && t == previousTime)
        return;

    for (auto &pole : resolved) std::fill(pole.begin(), pole.end(), 0);

    // n_re includes all of f_re, not just its above-threshold part. Keep
    // the raw difference for diagnostics and floor only the source weight.
    const real_t *nRe = unknowns->GetUnknownDataPrevious(idNRe);
    std::copy(nRe, nRe + externalRaw.size(), externalRaw.begin());
    if (runaway) {
        const real_t *fRe = unknowns->GetUnknownDataPrevious(idFRe);
        for (len_t ir = 0, offset = 0; ir < runaway->GetNr(); ++ir) {
            externalRaw[ir] -= runaway->IntegralMomentumAtRadius(ir, fRe + offset);
            offset += runaway->GetMomentumGrid(ir)->GetNCells();
        }
    }
    for (len_t ir = 0; ir < external.size(); ++ir)
        external[ir] = std::max(real_t(0), externalRaw[ir]);

    for (const auto &primary : primaryGrids) {
        const real_t *f = unknowns->GetUnknownDataPrevious(primary.id);
        len_t offset = 0, pitchOffset = 0, momentumOffset = 0;

        for (len_t ir = 0; ir < primary.grid->GetNr(); ++ir) {
            const auto *mg = primary.grid->GetMomentumGrid(ir);
            const len_t np = mg->GetNp1(), nxi = mg->GetNp2();
            const real_t *volume = primary.grid->GetVp(ir);
            const real_t radialVolume = primary.grid->GetVpVol(ir);
            for (len_t j = 0; j < nxi; ++j) {
                if (primary.grid->GetVpOverP2AtZero(ir)[j] == 0)
                    continue;
                const real_t positive = primary.positiveFraction[pitchOffset+j];
                const real_t pitchWeight = mg->GetDp2(j) / radialVolume;
                real_t admitted = 0;
                for (len_t k = 0; k < np; ++k) {
                    const len_t cell = j*np + k;
                    const real_t densityPerP = f[offset+cell] * volume[cell] * pitchWeight;
                    admitted += densityPerP * primary.momentumWidth[momentumOffset+k];
                }
                resolved[0][ir] += positive * admitted;
                resolved[1][ir] += (1 - positive) * admitted;
            }
            offset += np*nxi;
            pitchOffset += nxi;
            momentumOffset += np;
        }
    }
    previousTime = t;
    densitiesValid = true;
}

RosenbluthPutvinskiBoltzmannTerm::RosenbluthPutvinskiBoltzmannTerm(
    FVM::Grid *grid, FVM::UnknownQuantityHandler *unknowns,
    RosenbluthPutvinskiBoltzmann *source, RosenbluthPutvinskiBoltzmann::Region region,
    real_t scaleFactor
) : FVM::EquationTerm(grid), unknowns(unknowns), source(source), region(region),
    scaleFactor(scaleFactor),
    kinetic(region == RosenbluthPutvinskiBoltzmann::HOTTAIL || region == RosenbluthPutvinskiBoltzmann::RUNAWAY) {
    SetName("RosenbluthPutvinskiBoltzmannTerm");
}

void RosenbluthPutvinskiBoltzmannTerm::Rebuild(real_t, real_t, FVM::UnknownQuantityHandler *u) {
    if (u)
        unknowns = u;
    source->RebuildPrimaryDensities(unknowns);
}

/**
 * Use the cached previous-time primary densities. The first external
 * closure is a positive-pole beam, independent of the electric field.
 * Resolved RP primaries retain their two-pole passing/trapped split.
 */
template<typename Add>
void RosenbluthPutvinskiBoltzmannTerm::Apply(
    const real_t *nTot, Population population, Add add
) const {
    const real_t *plus = kinetic ? source->GetSource(region, true) : nullptr;
    const real_t *minus = kinetic && source->IncludesResolvedPrimaries()
        ? source->GetSource(region, false) : nullptr;
    const real_t *resolvedPlus = source->GetResolvedDensity(true);
    const real_t *resolvedMinus = source->GetResolvedDensity(false);
    const real_t *external = source->GetExternalDensity();

    for (len_t ir = 0, offset = 0; ir < grid->GetNr(); ++ir) {
        real_t positive = 0, negative = 0;
        if (population != OFF_GRID) {
            positive = resolvedPlus[ir];
            negative = resolvedMinus[ir];
        }
        if (population != RESOLVED)
            positive += external[ir];

        const real_t scale = nTot ? scaleFactor * nTot[ir] : scaleFactor;
        if (!kinetic) {
            add(ir, ir, scale * source->GetProduction(region) * (positive + negative));
            continue;
        }
        const len_t end = offset + grid->GetMomentumGrid(ir)->GetNCells();
        for (len_t cell = offset; cell < end; ++cell)
            add(cell, ir, scale * (positive * plus[cell] + (minus ? negative * minus[cell] : 0)));
        offset = end;
    }
}

void RosenbluthPutvinskiBoltzmannTerm::SetVectorElements(real_t *vec, const real_t *nTot) {
    Apply(nTot, BOTH, [&](len_t row, len_t, real_t value) { vec[row] += value; });
}

void RosenbluthPutvinskiBoltzmannTerm::SetMatrixElements(FVM::Matrix *mat, real_t *) {
    Apply(nullptr, BOTH, [&](len_t row, len_t col, real_t value) { mat->SetElement(row, col, value); });
}

bool RosenbluthPutvinskiBoltzmannTerm::SetJacobianBlock(
    len_t uqtyId, len_t derivId, FVM::Matrix *mat, const real_t *
) {
    if (derivId != uqtyId)
        return false;
    SetMatrixElements(mat, nullptr);
    return true;
}

void RosenbluthPutvinskiBoltzmannTerm::SetResolvedSource(real_t *vec, const real_t *nTot) const {
    Apply(nTot, RESOLVED, [&](len_t row, len_t, real_t value) { vec[row] += value; });
}

void RosenbluthPutvinskiBoltzmannTerm::SetExternalSource(real_t *vec, const real_t *nTot) const {
    Apply(nTot, OFF_GRID, [&](len_t row, len_t, real_t value) { vec[row] += value; });
}

void RosenbluthPutvinskiBoltzmannTerm::GetPrimaryDensity(real_t *density) const {
    const real_t *resolvedPlus = source->GetResolvedDensity(true);
    const real_t *resolvedMinus = source->GetResolvedDensity(false);
    GetExternalDensity(density);
    for (len_t ir = 0; ir < grid->GetNr(); ++ir)
        density[ir] += resolvedPlus[ir] + resolvedMinus[ir];
}

void RosenbluthPutvinskiBoltzmannTerm::GetExternalDensity(real_t *density) const {
    const real_t *external = source->GetExternalDensity();
    std::copy(external, external + grid->GetNr(), density);
}

void RosenbluthPutvinskiBoltzmannTerm::GetRawExternalDensity(real_t *density) const {
    const real_t *raw = source->GetRawExternalDensity();
    std::copy(raw, raw + grid->GetNr(), density);
}

/** Physical birth rate beyond the outermost kinetic face (not a flux). */
void RosenbluthPutvinskiBoltzmannTerm::GetAboveGridProduction(real_t *rate, const real_t *nTot) const {
    GetPrimaryDensity(rate);
    for (len_t ir = 0; ir < grid->GetNr(); ++ir)
        rate[ir] *= nTot[ir] * source->GetProduction(RosenbluthPutvinskiBoltzmann::ABOVE_GRID);
}

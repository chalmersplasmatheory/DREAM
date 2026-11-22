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
 * Prepare the positive-pole spectrum for off-grid primaries.
 * The highest kinetic grid sets the start of the off-grid population.
 */
RosenbluthPutvinskiBoltzmann::RosenbluthPutvinskiBoltzmann(
    const FVM::Grid *hottail, const FVM::Grid *runaway, FVM::UnknownQuantityHandler *unknowns,
    real_t pCutoff, len_t poloidalOrder
) : hottail(hottail), runaway(runaway), unknowns(unknowns),
    idNRe(unknowns->GetUnknownID(OptionConstants::UQTY_N_RE)),
    pCutoff(pCutoff), poloidalOrder(poloidalOrder) {
    if ((!hottail && !runaway) || poloidalOrder < 1
        || !std::isfinite(pCutoff) || pCutoff <= 5*std::numeric_limits<real_t>::epsilon())
        throw DREAMException("Boltzmann RP spectrum requires a kinetic grid, poloidal order >= 1 and a finite positive outgoing cutoff.");

    if (hottail && runaway && hottail->GetNr() != runaway->GetNr())
        throw DREAMException("Boltzmann RP requires the kinetic grids to share the same radii.");

    if (runaway)
        idFRe = unknowns->GetUnknownID(OptionConstants::UQTY_F_RE);

    GridRebuilt();
}

void RosenbluthPutvinskiBoltzmann::GridRebuilt() {
    const len_t nr = (runaway ? runaway : hottail)->GetNr();
    external.resize(nr);
    externalRaw.resize(nr);
    densitiesValid = false;

    if (hottail) BuildSource(HOTTAIL, hottail);
    if (runaway) BuildSource(RUNAWAY, runaway);
    const auto *outer = (runaway ? runaway : hottail)->GetMomentumGrid(0);
    const real_t pMax = outer->GetP1_f(outer->GetNp1());
    production[ABOVE_GRID] = IntegrateRPSpectrum(std::max(pCutoff, pMax), infiniteMomentum);
    production[RUNAWAY_DENSITY] = production[RUNAWAY] + production[ABOVE_GRID];
    production[ALL] = production[HOTTAIL] + production[RUNAWAY_DENSITY];
}

void RosenbluthPutvinskiBoltzmann::BuildSource(Region region, const FVM::Grid *grid) {
    source[region].resize(grid->GetNCells());

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
            if (spectrum[i] > 0) {
                // One scattering cosine per outgoing p cell in the
                // infinite-energy limit. Tabulate directly on that grid.
                KnockOnUtilities::SetFieldAlignedDeltaColumnOnGrid(
                    ir, cosine[i], +1, grid, column.data(), geometry, rule
                );
            }
            for (len_t j = 0; j < nxi; ++j) {
                const len_t cell = j*np + i;
                source[region][ir*np*nxi + cell] = spectrum[i] > 0 && volume[cell] > 0
                    ? spectrum[i] * column[j] * grid->GetVpVol(ir) / volume[cell] : 0;
            }
        }
    }
}

const real_t *RosenbluthPutvinskiBoltzmann::GetSource(Region region) const {
    return source[region].data();
}

/**
 * Use one previous-time moment for the kinetic gain and every fluid balance.
 * Remove the full f_re inventory from n_re at the same time level.
 */
void RosenbluthPutvinskiBoltzmann::RebuildPrimaryDensities(FVM::UnknownQuantityHandler *u) {
    if (u != unknowns) {
        unknowns = u;
        densitiesValid = false;
    }
    const real_t t = unknowns->GetUnknownDataPreviousTime(idNRe);
    if (densitiesValid && t == previousTime)
        return;

    // Keep the raw difference for diagnostics and floor the source weight.
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
 */
template<typename Add>
void RosenbluthPutvinskiBoltzmannTerm::Apply(const real_t *nTot, Add add) const {
    const real_t *shape = kinetic ? source->GetSource(region) : nullptr;
    const real_t *external = source->GetExternalDensity();

    for (len_t ir = 0, offset = 0; ir < grid->GetNr(); ++ir) {
        const real_t scale = nTot ? scaleFactor * nTot[ir] : scaleFactor;
        const real_t weight = scale * external[ir];
        if (!kinetic) {
            add(ir, ir, weight * source->GetProduction(region));
            continue;
        }
        const len_t end = offset + grid->GetMomentumGrid(ir)->GetNCells();
        for (len_t cell = offset; cell < end; ++cell)
            add(cell, ir, weight * shape[cell]);
        offset = end;
    }
}

void RosenbluthPutvinskiBoltzmannTerm::SetVectorElements(real_t *vec, const real_t *nTot) {
    SetExternalSource(vec, nTot);
}

void RosenbluthPutvinskiBoltzmannTerm::SetMatrixElements(FVM::Matrix *mat, real_t *) {
    Apply(nullptr, [&](len_t row, len_t col, real_t value) { mat->SetElement(row, col, value); });
}

bool RosenbluthPutvinskiBoltzmannTerm::SetJacobianBlock(
    len_t uqtyId, len_t derivId, FVM::Matrix *mat, const real_t *
) {
    if (derivId != uqtyId)
        return false;
    SetMatrixElements(mat, nullptr);
    return true;
}

void RosenbluthPutvinskiBoltzmannTerm::SetExternalSource(real_t *vec, const real_t *nTot) const {
    Apply(nTot, [&](len_t row, len_t, real_t value) { vec[row] += value; });
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
    GetExternalDensity(rate);
    for (len_t ir = 0; ir < grid->GetNr(); ++ir)
        rate[ir] *= nTot[ir] * source->GetProduction(RosenbluthPutvinskiBoltzmann::ABOVE_GRID);
}

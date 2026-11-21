/**
 * A local continuation at the last momentum face, with a finite reservoir.
 * This is a boundary closure, not an evolved off-grid distribution.
 */
#include <algorithm>
#include <array>
#include "FVM/Grid/FluxSurfaceAverager.hpp"
#include "DREAM/Equations/Kinetic/ExternalKineticBoundary.hpp"

using namespace DREAM;

ExternalKineticBoundary::ExternalKineticBoundary(
    FVM::Grid *g, FVM::Grid *distributionGrid, const FVM::Operator *op,
    len_t idF, std::array<len_t, 2> idExternal,
    FVM::UnknownQuantityHandler *unknowns, Side side
) : PXiAdvectionDiffusionBoundaryCondition(g, op), distributionGrid(distributionGrid),
    idF(idF), idExternal(idExternal), unknowns(unknowns), side(side) {
    SetName("ExternalKineticBoundary");
    GridRebuilt();
}

bool ExternalKineticBoundary::GridRebuilt() {
    pitchOffset.resize(distributionGrid->GetNr());
    for (len_t ir = 0, offset = 0; ir < distributionGrid->GetNr(); ++ir) {
        if (distributionGrid->GetNp1(ir) < 2)
            throw FVM::OperatorException("The external kinetic boundary requires at least two momentum cells.");
        pitchOffset[ir] = offset;
        offset += distributionGrid->GetNp2(ir);
    }
    // Split each face by its pitch orbit measure. Share the trapped
    // contribution evenly between both poles.
    BuildPoleFractions();
    return PXiAdvectionDiffusionBoundaryCondition::GridRebuilt();
}

void ExternalKineticBoundary::BuildPoleFractions() {
    auto *grid = distributionGrid;
    auto *rg = grid->GetRadialGrid();
    positiveFraction.clear();
    for (len_t ir = 0; ir < grid->GetNr(); ++ir) {
        const auto *mg = grid->GetMomentumGrid(ir);
        auto *fsa = rg->GetFluxSurfaceAverager();
        const real_t xiT = rg->GetXi0TrappedBoundary(ir);
        const auto measure = [&](real_t lower, real_t upper) -> real_t {
            if (upper <= lower)
                return 0;
            return (upper - lower) * fsa->EvaluateCellAveragedBounceIntegralOverP2(
                ir, lower, upper, FVM::FLUXGRIDTYPE_DISTRIBUTION,
                FVM::RadialGrid::BA_FUNC_UNITY, nullptr, FVM::RadialGrid::BA_PARAM_UNITY
            );
        };

        for (len_t l = 0; l < mg->GetNp2(); ++l) {
            const real_t lower = mg->GetP2_f(l);
            const real_t upper = mg->GetP2_f(l + 1);
            real_t fraction;
            if (grid->GetVpOverP2AtZero(ir)[l] == 0 || upper <= -xiT) {
                // A void row contributes no primary inventory.
                fraction = 0;
            } else if (lower >= xiT) {
                fraction = 1;
            } else if (lower >= -xiT && upper <= xiT) {
                fraction = 0.5;
            } else if (xiT == 0) {
                // Constant field and a pitch cell that straddles zero.
                fraction = upper / (upper - lower);
            } else {
                // The passing measure is even in pitch. Reflect this piece so
                // the native integral includes it up to the separatrix.
                const real_t counterPassing = measure(std::max(-upper, xiT), -lower);
                const real_t trapped = measure(std::max(lower, real_t(0)), std::min(upper, xiT));
                const real_t coPassing = measure(std::max(lower, xiT), upper);
                const real_t total = counterPassing + trapped + coPassing;
                if (!(total > 0) || !std::isfinite(total))
                    throw FVM::OperatorException("External boundary pitch cell has invalid orbit measure for the pole split.");
                fraction = (coPassing + 0.5 * trapped) / total;
            }
            positiveFraction.push_back(fraction);
        }
    }
}

bool ExternalKineticBoundary::Rebuild(real_t t, FVM::UnknownQuantityHandler *u) {
    unknowns = u;
    // BoundaryCondition::Rebuild has no dt argument. The previous accepted
    // time gives the actual step also after a rejected/adaptive step.
    dt = std::max(real_t(0), t - u->GetUnknownDataPreviousTime(idF));
    // SaveStep advances DREAM's accepted buffers before other quantities are
    // written. Keep this step's incoming trace for those flux diagnostics too.
    const auto *accepted = u->GetUnknownDataPrevious(idF);
    acceptedBoundaryValues.resize(distributionGrid->GetNr());
    for (len_t ir=0, offset=0; ir<distributionGrid->GetNr(); ++ir) {
        const len_t np = distributionGrid->GetNp1(ir), nxi = distributionGrid->GetNp2(ir);
        auto &values = acceptedBoundaryValues[ir];
        values.resize(2*nxi);
        for (len_t j=0; j<nxi; ++j) {
            const len_t cell = offset+j*np+np-1;
            values[2*j] = accepted[cell-1];
            values[2*j+1] = accepted[cell];
        }
        offset += np*nxi;
    }
    return false;
}

len_t ExternalKineticBoundary::GetNumberOfNonZerosPerRow() const {
    len_t maximum = 0;
    for (len_t ir = 0; ir < distributionGrid->GetNr(); ++ir)
        maximum = std::max(maximum, 2*distributionGrid->GetNp2(ir));
    return side == KINETIC ? std::min(len_t(6), maximum) : maximum;
}

len_t ExternalKineticBoundary::GetNumberOfNonZerosPerRow_jac() const {
    const auto *ad = oprtr->GetAdvectionDiffusion();
    return GetNumberOfNonZerosPerRow() + idExternal.size()
        + ad->GetNumberOfNonZerosPerRow_jac() - ad->GetNumberOfNonZerosPerRow();
}

/**
 * Reconstruct the diffusion trace at the face. Clip the slope so the
 * trace stays non-increasing and nonnegative, and return derivatives of
 * the selected limiter branch.
 */
ExternalKineticBoundary::Trace ExternalKineticBoundary::Reconstruct(
    real_t previous, real_t last, real_t distance, real_t halfWidth
) const {
    const real_t value = std::max(real_t(0), last);
    const real_t dLast = last > 0 ? 1 : 0;
    Trace trace{value, 0, {0, dLast}, {0, 0}};
    const real_t slope = (value - std::max(real_t(0), previous)) / distance;
    if (slope >= 0)
        return trace;
    if (slope <= -value/halfWidth) {
        trace.value = 0;
        trace.slope = -value/halfWidth;
        trace.dValue = {0, 0};
        trace.dSlope = {0, -dLast/halfWidth};
    } else {
        trace.slope = slope;
        trace.dSlope = {previous > 0 ? -1/distance : 0, dLast/distance};
        trace.value += halfWidth*slope;
        trace.dValue = {halfWidth*trace.dSlope[0], dLast + halfWidth*trace.dSlope[1]};
    }
    return trace;
}

/** Build physical face fluxes A*f - Dpp*df/dp - Dpx*df/dxi. */
std::vector<ExternalKineticBoundary::Face> ExternalKineticBoundary::EvaluateRawFaces(
    len_t ir, len_t offset, const real_t *f, bool accepted
) const {
    const auto *mg = distributionGrid->GetMomentumGrid(ir);
    const len_t np = mg->GetNp1(), nxi = mg->GetNp2();
    const real_t distance = mg->GetP1(np-1) - mg->GetP1(np-2);
    const real_t halfWidth = mg->GetP1_f(np) - mg->GetP1(np-1);
    const auto *A = oprtr->GetAdvectionCoeff1(ir);
    const auto *D = oprtr->GetDiffusionCoeff11(ir);
    const auto *Dx = oprtr->GetDiffusionCoeff12(ir);
    std::vector<Trace> traces(nxi);
    for (len_t j = 0; j < nxi; ++j) {
        const len_t cell = offset + j*np + np-1;
        const real_t previous = accepted ? acceptedBoundaryValues.at(ir)[2*j] : f[cell-1];
        const real_t last = accepted ? acceptedBoundaryValues.at(ir)[2*j+1] : f[cell];
        traces[j] = Reconstruct(previous, last, distance, halfWidth);
    }

    std::vector<Face> faces(nxi);
    for (len_t j = 0; j < nxi; ++j) {
        // These rows are constrained by the trapped-particle mirror, and
        // have no independent inventory. Omit them on BOTH sides.
        if (distributionGrid->IsNegativePitchTrappedIgnorableCell(ir, j))
            continue;
        const len_t cell = j*np + np-1, face = j*(np+1) + np;
        const real_t volume = distributionGrid->GetVp(ir)[cell];
        if (volume <= 0)
            continue;
        auto &result = faces[j];
        const auto &trace = traces[j];
        const real_t a = A ? A[face] : 0;
        const real_t d = D ? D[face] : 0;
        const real_t dx = Dx ? Dx[face] : 0;
        const real_t faceVolume = distributionGrid->GetVp_f1(ir)[face];
        result.weight = faceVolume * mg->GetDp2(j) / distributionGrid->GetVpVol(ir);
        result.kineticScale = faceVolume / (volume*mg->GetDp1(np-1));
        const real_t positive = positiveFraction[pitchOffset[ir] + j];
        result.fraction = {positive, 1-positive};
        // Use the adjacent cell for both advection directions. Keep the
        // reconstructed slope and face trace for the diffusive flux.
        result.advectionValue = std::max(real_t(0), accepted
            ? acceptedBoundaryValues.at(ir)[2*j+1] : f[offset+cell]);
        result.slope = trace.slope;
        result.flux = a*result.advectionValue - d*trace.slope;
        for (len_t k = 0; k < 2; ++k) {
            const real_t dAdvectionValue = k == 1 && result.advectionValue > 0 ? 1 : 0;
            result.derivative.emplace_back(2*j+k, a*dAdvectionValue - d*trace.dSlope[k]);
        }
        if (j > 0 && j+1 < nxi) {
            const real_t width = mg->GetP2(j+1) - mg->GetP2(j-1);
            result.pitchSlope = (traces[j+1].value - traces[j-1].value) / width;
            result.flux -= dx*result.pitchSlope;
            for (len_t k = 0; k < 2; ++k) {
                result.derivative.emplace_back(2*(j+1)+k, -dx*traces[j+1].dValue[k]/width);
                result.derivative.emplace_back(2*(j-1)+k, +dx*traces[j-1].dValue[k]/width);
            }
        }
    }
    return faces;
}

/** Split outgoing and incoming transport before applying the reservoir limit.
 * The accepted buffer is managed by DREAM's timestep save/rollback machinery.
 */
std::vector<ExternalKineticBoundary::Face> ExternalKineticBoundary::EvaluateFaces(
    len_t ir, len_t offset, const real_t *f
) const {
    auto faces = EvaluateRawFaces(ir, offset, f);
    const auto accepted = EvaluateRawFaces(ir, offset, f, true);
    for (len_t j=0; j<faces.size(); ++j) {
        faces[j].acceptedFlux = accepted[j].flux;
        faces[j].inflow = std::min(real_t(0), accepted[j].flux);
        faces[j].inflowAdvectionValue = accepted[j].advectionValue;
        faces[j].inflowSlope = accepted[j].slope;
        faces[j].inflowPitchSlope = accepted[j].pitchSlope;
    }
    return faces;
}

ExternalKineticBoundary::PoleValues ExternalKineticBoundary::RequestedInflow(
    const std::vector<Face>& faces
) {
    PoleValues rates{};
    for (const auto &face : faces)
        for (len_t pole = 0; pole < rates.size(); ++pole)
            rates[pole] -= face.weight*face.fraction[pole]*face.inflow;
    return rates;
}

/**
 * Limit the inward rate of each pole independently at each radius.
 * With I the requested rate from that pole and n its current inventory,
 *     lambda = n / (n + dt*I),  I_actual = lambda*I <= n/dt.
 * The requested rate uses the accepted trace and current coefficients.
 */
ExternalKineticBoundary::Reservoir ExternalKineticBoundary::LimitInflow(
    real_t rawDensity, real_t rate
) const {
    const real_t density = std::max(real_t(0), rawDensity);
    if (dt == 0)
        return {density > 0 ? real_t(1) : real_t(0), 0, 0};
    const real_t requestedDensity = dt*rate;
    const real_t denominator = density + requestedDensity;
    if (denominator == 0)
        return {0, 0, 0};
    const real_t scale = density / denominator;
    // Choose the right derivative at zero, so an empty reservoir does not
    // remove its coupling from the Newton system.
    const real_t dDensity = rawDensity >= 0 ? (requestedDensity/denominator)/denominator : 0;
    return {scale, dDensity, -scale*(dt/denominator)};
}

ExternalKineticBoundary::Reservoirs ExternalKineticBoundary::LimitInflow(
    len_t ir, const PoleValues& rates
) const {
    Reservoirs reservoirs;
    for (len_t pole = 0; pole < reservoirs.size(); ++pole)
        reservoirs[pole] = LimitInflow(unknowns->GetUnknownData(idExternal[pole])[ir], rates[pole]);
    return reservoirs;
}

real_t ExternalKineticBoundary::AppliedFlux(const Face& face, const Reservoirs& reservoirs) const {
    const real_t outward = std::max(real_t(0), face.flux);
    real_t flux = 0;
    for (len_t pole = 0; pole < reservoirs.size(); ++pole)
        if (IncludesPole(pole))
            flux += face.fraction[pole]*(outward + reservoirs[pole].scale*face.inflow);
    return flux;
}

len_t ExternalKineticBoundary::Row(len_t ir, len_t offset, len_t j) const {
    const len_t np = distributionGrid->GetNp1(ir);
    return side == KINETIC ? offset+j*np+np-1 : ir;
}

real_t ExternalKineticBoundary::Scale(const Face& face) const {
    return side == KINETIC ? face.kineticScale : -face.weight;
}

void ExternalKineticBoundary::AddToVectorElements(real_t *vec, const real_t *f) {
    for (len_t ir = 0, offset = 0; ir < grid->GetNr(); ++ir) {
        const auto faces = EvaluateFaces(ir, offset, f);
        const auto reservoirs = LimitInflow(ir, RequestedInflow(faces));
        for (len_t j = 0; j < faces.size(); ++j)
            vec[Row(ir, offset, j)] += Scale(faces[j])*AppliedFlux(faces[j], reservoirs);
        offset += distributionGrid->GetMomentumGrid(ir)->GetNCells();
    }
}

/**
 * At fixed coefficients and inventory, only the outgoing part depends on
 * the current distribution. Retain explicit zeros in the boundary pattern.
 */
void ExternalKineticBoundary::AddDistributionJacobian(
    len_t ir, len_t offset, const std::vector<Face>& faces, FVM::Matrix *jac
) const {
    const len_t np = distributionGrid->GetNp1(ir), nxi = faces.size();
    if (side != KINETIC) {
        std::vector<PetscInt> columns(2*nxi);
        std::vector<real_t> values(2*nxi);
        for (len_t j = 0; j < nxi; ++j) {
            columns[2*j] = offset+j*np+np-2;
            columns[2*j+1] = columns[2*j]+1;
        }
        for (const auto &face : faces) {
            const real_t fraction = side == NEGATIVE_DENSITY ? face.fraction[NEGATIVE] : 1;
            const real_t factor = face.flux > 0 ? -face.weight*fraction : 0;
            for (const auto &[column, derivative] : face.derivative)
                values[column] += factor*derivative;
        }
        jac->SetRow(ir, columns.size(), columns.data(), values.data());
        return;
    }

    for (len_t j = 0; j < nxi; ++j) {
        const auto &face = faces[j];
        const len_t count = face.derivative.size();
        if (count == 0)
            continue;

        const real_t factor = face.flux > 0 ? face.kineticScale : 0;
        std::array<PetscInt, 6> columns;
        std::array<real_t, 6> values;
        for (len_t k = 0; k < count; ++k) {
            const auto &[column, derivative] = face.derivative[k];
            columns[k] = offset+(column/2)*np+np-2+(column%2);
            values[k] = factor*derivative;
        }
        // Keep the local stencil allocated even when the flux is inward.
        jac->SetRow(Row(ir, offset, j), count, columns.data(), values.data());
    }
}

bool ExternalKineticBoundary::AddToJacobianBlock(
    len_t qtyId, len_t derivId, FVM::Matrix *jac, const real_t *f
) {
    const auto external = std::find(idExternal.begin(), idExternal.end(), derivId);
    if (derivId == idF || external != idExternal.end()) {
        const len_t pole = std::distance(idExternal.begin(), external);
        if (derivId != idF && !IncludesPole(pole))
            return false;
        for (len_t ir = 0, offset = 0; ir < grid->GetNr(); ++ir) {
            const auto faces = EvaluateFaces(ir, offset, f);
            if (derivId == idF)
                AddDistributionJacobian(ir, offset, faces, jac);
            else {
                const real_t rate = RequestedInflow(faces)[pole];
                const auto reservoir = LimitInflow(unknowns->GetUnknownData(derivId)[ir], rate);
                PetscInt column = ir;
                if (side != KINETIC) {
                    const real_t value = rate*reservoir.dDensity;
                    jac->SetRow(ir, 1, &column, &value);
                } else {
                    for (len_t j = 0; j < faces.size(); ++j) {
                        const auto &face = faces[j];
                        if (face.derivative.empty())
                            continue;
                        const real_t value = face.kineticScale*face.fraction[pole]*face.inflow*reservoir.dDensity;
                        jac->SetRow(Row(ir, offset, j), 1, &column, &value);
                    }
                }
            }
            offset += distributionGrid->GetMomentumGrid(ir)->GetNCells();
        }
        return true;
    }
    // The base helper inserts only nonzero coefficient derivatives. Reserve
    // their boundary entries as well, since fluxes can change sign during a run.
    const auto *ad = oprtr->GetAdvectionDiffusion();
    len_t multiples = 0, count;
    for (auto *term : ad->GetAdvectionTerms())
        if (term->HasJacobianContribution(derivId, &count))
            multiples = std::max(multiples, count);
    for (auto *term : ad->GetDiffusionTerms())
        if (term->HasJacobianContribution(derivId, &count))
            multiples = std::max(multiples, count);
    const len_t nr = distributionGrid->GetNr();
    const real_t zero = 0;
    for (len_t ir = 0, offset = 0; ir < nr; ++ir) {
        for (len_t multiple = 0; multiple < multiples; ++multiple) {
            PetscInt column = multiple*nr + ir;
            if (side != KINETIC)
                jac->SetRow(ir, 1, &column, &zero);
            else
                for (len_t j = 0; j < distributionGrid->GetNp2(ir); ++j)
                    if (!distributionGrid->IsNegativePitchTrappedIgnorableCell(ir, j))
                        jac->SetRow(Row(ir, offset, j), 1, &column, &zero);
        }
        offset += distributionGrid->GetMomentumGrid(ir)->GetNCells();
    }
    return AddPartialJacobianContributions(qtyId, derivId, jac, f, false);
}

/**
 * At fixed reservoir densities, coefficient changes affect both the local
 * inward flux and the integrated request of each pole.
 */
void ExternalKineticBoundary::AddToVectorElements_c(
    real_t *vec, const real_t *f, const real_t *const*, const real_t *const* dA,
    const real_t *const*, const real_t *const*, const real_t *const* dD,
    const real_t *const* dDx, const real_t *const*, const real_t *const*,
    jacobian_interp_mode
) {
    for (len_t ir = 0, offset = 0; ir < grid->GetNr(); ++ir) {
        const auto faces = EvaluateFaces(ir, offset, f);
        const auto reservoirs = LimitInflow(ir, RequestedInflow(faces));
        const len_t np = distributionGrid->GetNp1(ir);
        std::vector<real_t> outwardDerivative(faces.size()), inwardDerivative(faces.size());
        PoleValues dRequested{};
        for (len_t j = 0; j < faces.size(); ++j) {
            if (faces[j].derivative.empty())
                continue;
            const len_t index = j*(np+1)+np;
            if (faces[j].flux > 0)
                outwardDerivative[j] = (dA ? dA[ir][index]*faces[j].advectionValue : 0)
                    - (dD ? dD[ir][index]*faces[j].slope : 0)
                    - (dDx ? dDx[ir][index]*faces[j].pitchSlope : 0);
            if (faces[j].acceptedFlux < 0) {
                inwardDerivative[j] = (dA ? dA[ir][index]*faces[j].inflowAdvectionValue : 0)
                    - (dD ? dD[ir][index]*faces[j].inflowSlope : 0)
                    - (dDx ? dDx[ir][index]*faces[j].inflowPitchSlope : 0);
                for (len_t pole = 0; pole < dRequested.size(); ++pole)
                    dRequested[pole] -= faces[j].weight*faces[j].fraction[pole]*inwardDerivative[j];
            }
        }
        for (len_t j = 0; j < faces.size(); ++j) {
            const auto &face = faces[j];
            real_t value = 0;
            for (len_t pole = 0; pole < reservoirs.size(); ++pole) {
                if (!IncludesPole(pole))
                    continue;
                const auto &reservoir = reservoirs[pole];
                // The density row sums the inward flux, lambda*I, over pitch.
                const real_t inward = side == KINETIC
                    ? reservoir.scale*inwardDerivative[j] + face.inflow*reservoir.dRate*dRequested[pole]
                    : reservoir.scale*reservoir.scale*inwardDerivative[j];
                value += face.fraction[pole]*(outwardDerivative[j] + inward);
            }
            vec[Row(ir, offset, j)] += Scale(face)*value;
        }
        offset += distributionGrid->GetMomentumGrid(ir)->GetNCells();
    }
}

void ExternalKineticBoundary::AddToMatrixElements(FVM::Matrix *, real_t *) {
    throw FVM::OperatorException("The two-way external boundary requires the nonlinear solver.");
}

void ExternalKineticBoundary::GetFlux(real_t *net, real_t *outward, real_t *inward, Pole selected) const {
    const auto *f = unknowns->GetUnknownData(idF);
    for (len_t ir = 0, offset = 0; ir < grid->GetNr(); ++ir) {
        const auto faces = EvaluateFaces(ir, offset, f);
        const auto rates = RequestedInflow(faces);
        const auto reservoirs = LimitInflow(ir, rates);
        real_t in = 0;
        for (len_t pole = 0; pole < reservoirs.size(); ++pole)
            if (selected == BOTH || selected == pole)
                in += rates[pole]*reservoirs[pole].scale;
        real_t out = 0;
        for (const auto &face : faces)
            out += face.weight*std::max(real_t(0), face.flux)
                * (selected == BOTH ? 1 : face.fraction[selected]);
        if (net) net[ir] = out-in;
        if (outward) outward[ir] = out;
        if (inward) inward[ir] = in;
        offset += distributionGrid->GetMomentumGrid(ir)->GetNCells();
    }
}

/** Requested rate before inventory limiting, for output diagnostics. */
void ExternalKineticBoundary::GetRequestedInflow(real_t *rate, Pole selected) const {
    const auto *f = unknowns->GetUnknownData(idF);
    for (len_t ir = 0, offset = 0; ir < grid->GetNr(); ++ir) {
        const auto rates = RequestedInflow(EvaluateFaces(ir, offset, f));
        rate[ir] = selected == BOTH ? rates[POSITIVE]+rates[NEGATIVE] : rates[selected];
        offset += distributionGrid->GetMomentumGrid(ir)->GetNCells();
    }
}

/**
 * Contact-split quadrature for the knock-on pitch kernel.
 *
 * The primary pitch and outgoing cell faces are specified at Bmin.
 * At fixed primary pitch, contacts are expressed as field levels
 * beta = 1 - Bmin/B and then located on both poloidal branches.
 * The outer primary-pitch integral is split where those contacts
 * meet a field extremum or a bounce point. See doc/notes/theory.tex.
 */
#include "DREAM/Equations/KnockOnQuadrature.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <utility>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>

#include "DREAM/DREAMException.hpp"
#include "FVM/Grid/FluxSurfaceAverager.hpp"

namespace DREAM::KnockOnUtilities::Quadrature {

namespace {

void ValidatePitch(real_t pitch) {
    if (!std::isfinite(pitch) || pitch < -1 || pitch > 1)
        throw DREAMException("KnockOnQuadrature: pitch/cosine must be finite and in [-1,1].");
}

void ValidateBetaMax(real_t betaMax) {
    if (!std::isfinite(betaMax) || betaMax < 0 || betaMax >= 1)
        throw DREAMException("KnockOnQuadrature: betaMax must be finite and in [0,1).");
}

void ValidateInterval(Interval interval) {
    if (!std::isfinite(interval.lower) || !std::isfinite(interval.upper)
        || interval.lower > interval.upper || !std::isfinite(interval.upper - interval.lower)) {
        throw DREAMException("KnockOnQuadrature: invalid integration interval.");
    }
}

void SortUnique(std::vector<real_t> &points) {
    std::sort(points.begin(), points.end());
    points.erase(std::unique(points.begin(), points.end()), points.end());
}

}  // namespace

std::vector<real_t> SplitPoints(Interval interval, std::vector<real_t> cuts) {
    ValidateInterval(interval);
    for (real_t x : cuts) {
        if (!std::isfinite(x))
            throw DREAMException("KnockOnQuadrature: non-finite split point.");
    }

    const auto outside = [&](real_t x) {
        return x <= interval.lower || x >= interval.upper;
    };
    cuts.erase(std::remove_if(cuts.begin(), cuts.end(), outside), cuts.end());
    cuts.push_back(interval.lower);
    cuts.push_back(interval.upper);
    SortUnique(cuts);
    return cuts;
}

/**
 * Find field levels where the outgoing scattering cone touches one face.
 * A face can bounce before the primary, in which case its local pitch is
 * clipped to zero and the contact equation changes.
 */
std::vector<real_t> ContactFieldLevels(
    real_t primaryPitch, real_t secondaryFace, real_t xiStar, real_t betaMax
) {
    ValidatePitch(primaryPitch);
    ValidatePitch(secondaryFace);
    ValidatePitch(xiStar);
    ValidateBetaMax(betaMax);

    // Extra precision limits cancellation for nearly collinear scattering
    // and nearly equal invariant pitches.
    const long double primary = primaryPitch;
    const long double face = secondaryFace;
    const long double cosine = std::abs(static_cast<long double>(xiStar));
    const long double primarySquared = primary * primary;
    const long double faceSquared = face * face;
    const long double sinSquared = (1 - cosine) * (1 + cosine);
    const long double reachableLevel = std::min(primarySquared, static_cast<long double>(betaMax));

    std::vector<real_t> levels;
    const auto addLevel = [&](long double beta) {
        const real_t level = static_cast<real_t>(beta);
        // Check again after conversion: rounding can put a root on an endpoint.
        if (beta > 0 && beta < reachableLevel
            && level > 0 && level < std::min(primaryPitch * primaryPitch, betaMax)) {
            levels.push_back(level);
        }
    };

    // The face bounces at beta = v^2.
    addLevel(faceSquared);

    // Solve the squared contact equation for reachable faces. Extra roots
    // only split a smooth panel.
    if (cosine < 1) {
        const long double primarySin = std::sqrt((1 - primary) * (1 + primary));
        const long double faceSin = std::sqrt((1 - face) * (1 + face));
        const long double difference = primarySin - faceSin;
        for (long double denominator : {1 + cosine, 1 - cosine}) {
            const long double beta = 1 - difference * difference / sinSquared
                - 2 * primarySin * faceSin / denominator;
            if (beta <= std::min(primarySquared, faceSquared))
                addLevel(beta);
        }
    }

    // Contact with a face clipped to local pitch zero. Compare before dividing
    // to avoid overflow for arbitrarily small nonzero scattering cosine.
    const long double numerator = primarySquared - sinSquared;
    if (cosine > 0 && numerator > 0
        && numerator < cosine * cosine * reachableLevel
        && numerator >= cosine * cosine * faceSquared) {
        addLevel(numerator / (cosine * cosine));
    }

    SortUnique(levels);
    return levels;
}

/**
 * Find where the poloidal panel boundaries change as the primary pitch varies.
 */
std::vector<real_t> PrimaryPitchSplitPoints(
    Interval cell, real_t xiStar, real_t betaMax,
    const std::vector<real_t> &secondaryFaces
) {
    ValidateInterval(cell);
    ValidatePitch(cell.lower);
    ValidatePitch(cell.upper);
    ValidatePitch(xiStar);
    ValidateBetaMax(betaMax);

    const long double cosine = std::abs(static_cast<long double>(xiStar));
    const long double sinSquared = (1 - cosine) * (1 + cosine);
    const long double sine = std::sqrt(sinSquared);
    const real_t xiT = std::sqrt(betaMax);
    std::vector<real_t> cuts{-xiT, 0, xiT};
    const auto addSignedCuts = [&](long double pitchSquared) {
        if (pitchSquared >= 0 && pitchSquared <= 1) {
            const real_t pitch = static_cast<real_t>(std::sqrt(pitchSquared));
            cuts.push_back(-pitch);
            cuts.push_back(pitch);
        }
    };

    for (real_t secondaryFace : secondaryFaces) {
        ValidatePitch(secondaryFace);
        const long double face = secondaryFace;
        const long double faceSquared = face * face;

        // A contact enters or leaves the poloidal domain at a field extremum.
        for (long double beta : {0.L, static_cast<long double>(betaMax)}) {
            if (faceSquared < beta)
                continue;

            const long double localFace = std::copysign(
                std::sqrt((faceSquared - beta) / (1 - beta)), face
            );
            const long double halfWidth = sine * std::sqrt(
                std::max(0.L, (1 - localFace) * (1 + localFace))
            );
            for (long double localPitch : {
                cosine * localFace - halfWidth, cosine * localFace + halfWidth
            }) {
                addSignedCuts(beta + (1 - beta) * localPitch * localPitch);
            }
        }

        if (faceSquared <= betaMax) {
            // Contact at the face bounce, then coincident bounces.
            addSignedCuts(sinSquared + cosine * cosine * faceSquared);
            addSignedCuts(faceSquared);
        }
        const long double numerator = faceSquared - sinSquared;
        if (cosine > 0 && numerator >= 0 && numerator <= cosine * cosine * betaMax) {
            // The contact reaches the primary bounce.
            addSignedCuts(numerator / (cosine * cosine));
        }
    }

    return SplitPoints(cell, std::move(cuts));
}

/**
 * Construct a Gauss rule on [0,1] and its sin^2-mapped counterpart.
 */
GaussLegendreRule::GaussLegendreRule(len_t order) {
    if (order < 1)
        throw DREAMException("KnockOnQuadrature: Gauss order must be at least 1.");

    std::unique_ptr<gsl_integration_glfixed_table, decltype(&gsl_integration_glfixed_table_free)>
        table(gsl_integration_glfixed_table_alloc(order), &gsl_integration_glfixed_table_free);
    if (!table)
        throw DREAMException("KnockOnQuadrature: cannot allocate Gauss rule.");

    nodes.resize(order);
    weights.resize(order);
    mappedNodes.resize(order);
    mappedWeights.resize(order);
    real_t weightSum = 0;
    for (len_t i = 0; i < order; ++i) {
        if (gsl_integration_glfixed_point(0, 1, i, &nodes[i], &weights[i], table.get()) != GSL_SUCCESS)
            throw DREAMException("KnockOnQuadrature: cannot construct Gauss rule.");

        // Map from the nearer endpoint to retain accuracy near both ends.
        const real_t distance = std::min(nodes[i], 1 - nodes[i]);
        const real_t mappedDistance = std::pow(std::sin(M_PI * distance / 2), 2);
        mappedNodes[i] = nodes[i] < 0.5 ? mappedDistance : 1 - mappedDistance;
        mappedWeights[i] = weights[i] * (M_PI / 2) * std::sin(M_PI * nodes[i]);
        weightSum += mappedWeights[i];
    }

    // Preserve constants on each panel and make the midpoint exact.
    for (real_t &weight : mappedWeights) {
        weight /= weightSum;
    }
    if (order % 2)
        mappedNodes[order / 2] = 0.5;
}

std::vector<WeightedPoint> GaussLegendreRule::MappedPoints(Interval panel) const {
    ValidateInterval(panel);
    std::vector<WeightedPoint> points;

    // Keep rounded nodes inside the panel, away from singular endpoints.
    const real_t firstInterior = std::nextafter(panel.lower, panel.upper);
    const real_t lastInterior = std::nextafter(panel.upper, panel.lower);
    if (panel.lower == panel.upper || firstInterior >= panel.upper)
        return points;

    points.reserve(GetOrder());
    const real_t width = panel.upper - panel.lower;
    for (len_t i = 0; i < GetOrder(); ++i) {
        const real_t node = mappedNodes[i];
        const real_t x = node <= 0.5 ? panel.lower + width * node
            : panel.upper - width * (1 - node);
        points.push_back({std::max(firstInterior, std::min(lastInterior, x)),
            width * mappedWeights[i]});
    }
    return points;
}

// Each branch is stored in increasing B, from the field minimum outwards.
struct OrbitGeometry::Implementation {
    struct Branch {
        std::vector<real_t> theta;
        std::vector<real_t> B;
    };

    std::function<real_t(real_t)> evaluateB;
    std::array<Branch, 2> branches;
    FieldExtrema extrema;
    bool constant = false;
    real_t betaMax = 0;
    std::unique_ptr<gsl_root_fsolver, decltype(&gsl_root_fsolver_free)>
        solver{nullptr, &gsl_root_fsolver_free};

    real_t EvaluateField(real_t theta) const {
        const real_t B = evaluateB(theta);
        if (!std::isfinite(B) || B < 0)
            throw DREAMException("KnockOnQuadrature: invalid magnetic field evaluation.");
        return B;
    }

    struct RootParameters {
        Implementation *geometry;
        real_t targetField;
    };

    static double RootFunction(double theta, void *parameters) {
        const auto &params = *static_cast<RootParameters *>(parameters);
        return params.geometry->EvaluateField(theta) - params.targetField;
    }

    real_t InvertBranch(const Branch &branch, real_t targetField) {
        // Levels closer to an extremum than the field evaluation can resolve
        // have that extremum as their only representable field bracket.
        if (targetField <= branch.B.front())
            return branch.theta.front();
        if (targetField >= branch.B.back())
            return branch.theta.back();

        // Use the samples to locate a bracket, then invert the actual field.
        const auto sample = std::lower_bound(branch.B.begin(), branch.B.end(), targetField);
        const std::size_t i = std::max<std::size_t>(1, sample - branch.B.begin());
        real_t thetaLower = branch.theta[i - 1];
        real_t thetaUpper = branch.theta[i];
        if (thetaLower > thetaUpper)
            std::swap(thetaLower, thetaUpper);

        const real_t fLower = EvaluateField(thetaLower) - targetField;
        const real_t fUpper = EvaluateField(thetaUpper) - targetField;
        if (fLower == 0)
            return thetaLower;
        if (fUpper == 0)
            return thetaUpper;
        if ((fLower < 0) == (fUpper < 0))
            throw DREAMException("KnockOnQuadrature: field level is not bracketed.");

        RootParameters params{this, targetField};
        gsl_function function;
        function.function = &RootFunction;
        function.params = &params;
        if (gsl_root_fsolver_set(solver.get(), &function, thetaLower, thetaUpper) != GSL_SUCCESS)
            throw DREAMException("KnockOnQuadrature: cannot initialize field inversion.");

        for (len_t iteration = 0; iteration < 100; ++iteration) {
            if (gsl_root_fsolver_iterate(solver.get()) != GSL_SUCCESS)
                throw DREAMException("KnockOnQuadrature: field inversion failed.");
            thetaLower = gsl_root_fsolver_x_lower(solver.get());
            thetaUpper = gsl_root_fsolver_x_upper(solver.get());
            if (gsl_root_test_interval(thetaLower, thetaUpper, 1e-11, 1e-11) == GSL_SUCCESS)
                return gsl_root_fsolver_root(solver.get());
        }
        throw DREAMException("KnockOnQuadrature: field inversion did not converge.");
    }
};

OrbitGeometry::OrbitGeometry(
    FieldExtrema extrema, std::function<real_t(real_t)> evaluateB,
    const std::vector<real_t> &thetaSamples
) : impl(new Implementation()) {
    if (!evaluateB)
        throw DREAMException("KnockOnQuadrature: magnetic field callback is empty.");
    impl->evaluateB = std::move(evaluateB);
    impl->solver.reset(gsl_root_fsolver_alloc(gsl_root_fsolver_brent));
    if (!impl->solver)
        throw DREAMException("KnockOnQuadrature: cannot allocate field root solver.");
    Reset(extrema, thetaSamples);
}

namespace {

FieldExtrema GetExtrema(FVM::FluxSurfaceAverager *fsa, len_t ir) {
    if (!fsa)
        throw DREAMException("KnockOnQuadrature: null flux surface averager.");

    FieldExtrema extrema;
    extrema.Bmin = fsa->GetBmin(ir, FVM::FLUXGRIDTYPE_DISTRIBUTION, &extrema.thetaMin);
    extrema.Bmax = fsa->GetBmax(ir, FVM::FLUXGRIDTYPE_DISTRIBUTION, &extrema.thetaMax);
    return extrema;
}

std::vector<real_t> GetThetaSamples(FVM::FluxSurfaceAverager *fsa) {
    if (!fsa)
        throw DREAMException("KnockOnQuadrature: null flux surface averager.");
    const len_t nTheta = fsa->GetNTheta();
    if (nTheta == 0)
        return {};
    if (!fsa->GetTheta())
        throw DREAMException("KnockOnQuadrature: missing bracket sampling angles.");
    return {fsa->GetTheta(), fsa->GetTheta() + nTheta};
}

}  // namespace

OrbitGeometry::OrbitGeometry(FVM::FluxSurfaceAverager *fsa, len_t ir)
    : OrbitGeometry(
        GetExtrema(fsa, ir),
        [fsa, ir](real_t theta) {
            return fsa->BAtTheta(ir, theta, FVM::FLUXGRIDTYPE_DISTRIBUTION);
        },
        GetThetaSamples(fsa)
    ) {}

OrbitGeometry::~OrbitGeometry() = default;

/**
 * Prepare the field brackets for the two branches of the current surface.
 */
void OrbitGeometry::Reset(FieldExtrema extrema, const std::vector<real_t> &thetaSamples) {
    if (!std::isfinite(extrema.Bmin) || !std::isfinite(extrema.Bmax)
        || !std::isfinite(extrema.thetaMin) || !std::isfinite(extrema.thetaMax)
        || extrema.Bmin < 0 || extrema.Bmax < extrema.Bmin) {
        throw DREAMException("KnockOnQuadrature: invalid field extrema.");
    }

    const real_t tolerance = 100 * std::numeric_limits<real_t>::epsilon() * (1 + extrema.Bmax);
    const bool constant = extrema.Bmax - extrema.Bmin <= tolerance;
    if (!constant && extrema.Bmin == 0)
        throw DREAMException("KnockOnQuadrature: varying field requires positive Bmin.");

    // Unwrap the field minimum into [thetaMax - 2*pi, thetaMax].
    const real_t period = 2 * M_PI;
    real_t displacement = std::fmod(extrema.thetaMax - extrema.thetaMin, period);
    if (displacement < 0)
        displacement += period;
    extrema.thetaMin = extrema.thetaMax - displacement;
    if (!constant && !(displacement > 0 && displacement < period))
        throw DREAMException("KnockOnQuadrature: distinct field extrema need distinct angles.");

    const real_t betaMax = constant ? 0 : (extrema.Bmax - extrema.Bmin) / extrema.Bmax;
    ValidateBetaMax(betaMax);

    std::array<Implementation::Branch, 2> branches;
    if (constant) {
        impl->extrema = extrema;
        impl->constant = true;
        impl->betaMax = 0;
        impl->branches = std::move(branches);
        return;
    }

    const real_t thetaLeft = extrema.thetaMax - period;
    const real_t thetaMin = extrema.thetaMin;
    const real_t thetaRight = extrema.thetaMax;

    // Seed brackets on both branches even when the supplied quadrature
    // nodes cover only half the surface.
    std::vector<real_t> samples = thetaSamples;
    const len_t nBracketIntervals = 32;
    for (len_t i = 0; i <= nBracketIntervals; ++i) {
        samples.push_back(thetaLeft + (thetaMin - thetaLeft) * i / nBracketIntervals);
        samples.push_back(thetaMin + (thetaRight - thetaMin) * i / nBracketIntervals);
    }

    std::array<std::vector<real_t>, 2> branchAngles;
    for (real_t theta : samples) {
        if (!std::isfinite(theta))
            throw DREAMException("KnockOnQuadrature: non-finite bracket sample.");

        real_t unwrapped = std::fmod(theta - thetaLeft, period);
        if (unwrapped < 0)
            unwrapped += period;
        unwrapped += thetaLeft;

        if (unwrapped > thetaLeft && unwrapped < thetaMin) {
            branchAngles[0].push_back(unwrapped);
        } else if (unwrapped > thetaMin && unwrapped < thetaRight) {
            branchAngles[1].push_back(unwrapped);
        }
    }

    for (len_t side = 0; side < 2; ++side) {
        auto &angles = branchAngles[side];
        angles.push_back(thetaMin);
        angles.push_back(side == 0 ? thetaLeft : thetaRight);
        SortUnique(angles);
        if (side == 0)
            std::reverse(angles.begin(), angles.end());

        auto &branch = branches[side];
        for (real_t theta : angles) {
            const real_t B = impl->EvaluateField(theta);
            if (B < extrema.Bmin - tolerance || B > extrema.Bmax + tolerance)
                throw DREAMException("KnockOnQuadrature: field disagrees with its extrema.");

            if (!branch.B.empty()) {
                if (B < branch.B.back() - tolerance)
                    throw DREAMException("KnockOnQuadrature: expected a monotone magnetic branch.");
                // Skip roundoff plateaus so the field bracket stays ordered.
                if (B <= branch.B.back())
                    continue;
            }
            branch.theta.push_back(theta);
            branch.B.push_back(B);
        }
        if (branch.B.size() < 2)
            throw DREAMException("KnockOnQuadrature: degenerate magnetic branch.");
    }

    // Publish the brackets only after both branches have been validated.
    impl->extrema = extrema;
    impl->constant = false;
    impl->betaMax = betaMax;
    impl->branches = std::move(branches);
}

bool OrbitGeometry::IsConstantField() const { return impl->constant; }
real_t OrbitGeometry::GetBetaMax() const { return impl->betaMax; }
Interval OrbitGeometry::GetAngularDomain() const {
    return {impl->extrema.thetaMax - 2 * M_PI, impl->extrema.thetaMax};
}

Interval OrbitGeometry::InvertBeta(real_t beta) {
    if (!std::isfinite(beta) || beta < 0 || beta > impl->betaMax)
        throw DREAMException("KnockOnQuadrature: field level is outside this surface.");
    if (impl->constant)
        return GetAngularDomain();

    const real_t thetaMin = impl->extrema.thetaMin;
    if (beta == 0)
        return {thetaMin, thetaMin};
    if (beta == impl->betaMax)
        return GetAngularDomain();

    const real_t targetField = impl->extrema.Bmin / (1 - beta);
    return {
        impl->InvertBranch(impl->branches[0], targetField),
        impl->InvertBranch(impl->branches[1], targetField)
    };
}

std::vector<FieldLevel> OrbitGeometry::InvertLevels(std::vector<real_t> levels) {
    for (real_t beta : levels) {
        if (!std::isfinite(beta) || beta < 0 || beta > impl->betaMax)
            throw DREAMException("KnockOnQuadrature: field level is outside this surface.");
    }

    SortUnique(levels);
    std::vector<FieldLevel> roots;
    roots.reserve(levels.size());
    for (real_t beta : levels) {
        roots.push_back({beta, InvertBeta(beta)});
    }
    return roots;
}

Interval OrbitGeometry::GetOrbitBounds(real_t primaryPitch) {
    ValidatePitch(primaryPitch);
    if (impl->constant || primaryPitch * primaryPitch >= impl->betaMax)
        return GetAngularDomain();
    return InvertBeta(primaryPitch * primaryPitch);
}

Interval OrbitGeometry::GetReachableBounds(real_t primaryPitch, real_t lowerFace, real_t upperFace) {
    ValidatePitch(primaryPitch);
    ValidatePitch(lowerFace);
    ValidatePitch(upperFace);
    if (lowerFace > upperFace)
        throw DREAMException("KnockOnQuadrature: secondary cell faces are reversed.");

    const real_t reachableLevel = std::min(
        primaryPitch * primaryPitch, std::max(lowerFace * lowerFace, upperFace * upperFace)
    );
    if (impl->constant || reachableLevel >= impl->betaMax)
        return GetAngularDomain();
    return InvertBeta(reachableLevel);
}

/**
 * Build the poloidal panels for each outgoing cell at a given primary pitch.
 */
std::vector<std::vector<real_t>> PoloidalSplitPoints(
    OrbitGeometry &geometry, real_t primaryPitch, real_t xiStar,
    const std::vector<Interval> &cells
) {
    ValidatePitch(primaryPitch);
    ValidatePitch(xiStar);
    const real_t betaMax = geometry.GetBetaMax();
    const real_t primarySquared = primaryPitch * primaryPitch;
    const bool constant = geometry.IsConstantField();
    if (!constant && (primarySquared == betaMax || std::abs(primaryPitch) == std::sqrt(betaMax))) {
        throw DREAMException("KnockOnQuadrature: pointwise primary pitch is on the separatrix.");
    }

    // Collect contact and reachability levels per cell. Share the field
    // inversions across cells while keeping their integration cuts separate.
    std::vector<std::vector<real_t>> contacts(cells.size());
    const real_t orbitLevel = std::min(primarySquared, betaMax);
    std::vector<real_t> reachableLevels(cells.size());
    std::vector<real_t> fieldLevels{0, betaMax, orbitLevel};
    for (len_t j = 0; j < cells.size(); ++j) {
        const Interval &cell = cells[j];
        ValidatePitch(cell.lower);
        ValidatePitch(cell.upper);
        ValidateInterval(cell);

        reachableLevels[j] = std::min({
            primarySquared, betaMax,
            std::max(cell.lower * cell.lower, cell.upper * cell.upper)
        });
        fieldLevels.push_back(reachableLevels[j]);
        for (real_t face : {cell.lower, cell.upper}) {
            const auto candidates = ContactFieldLevels(primaryPitch, face, xiStar, betaMax);
            for (real_t beta : candidates) {
                if (beta < reachableLevels[j])
                    contacts[j].push_back(beta);
            }
        }
        fieldLevels.insert(fieldLevels.end(), contacts[j].begin(), contacts[j].end());
    }

    // Invert a field-level drop |primaryPitch^2 - betaMax| to resolve the
    // separatrix layer using the actual field on each branch.
    const real_t gap = std::abs(primarySquared - betaMax);
    const bool nearSeparatrix = !constant && gap > 0 && gap < 0.1 * betaMax;
    const real_t layerLevel = nearSeparatrix ? betaMax - gap : 0;
    if (nearSeparatrix)
        fieldLevels.push_back(layerLevel);

    const auto roots = geometry.InvertLevels(std::move(fieldLevels));
    const auto rootsAt = [&](real_t beta) {
        // Each requested level was included in the batch above.
        const auto root = std::lower_bound(
            roots.begin(), roots.end(), beta,
            [](const FieldLevel &level, real_t value) { return level.beta < value; }
        );
        return root->roots;
    };

    // The field minimum and graded separatrix cuts are shared by all cells.
    const auto domain = geometry.GetAngularDomain();
    const real_t thetaMin = constant ? (domain.lower + domain.upper) / 2 : rootsAt(0).lower;
    std::vector<real_t> sharedCuts{thetaMin};
    if (nearSeparatrix) {
        const auto layer = rootsAt(layerLevel);
        for (bool leftBranch : {true, false}) {
            const real_t endpoint = leftBranch ? domain.lower : domain.upper;
            const real_t width = leftBranch ? layer.lower - endpoint : endpoint - layer.upper;
            const real_t branchWidth = std::abs(thetaMin - endpoint);
            // A level indistinguishable from Bmax gives zero width.
            for (real_t distance = width; distance > 0 && distance < branchWidth; distance *= 4) {
                sharedCuts.push_back(leftBranch ? endpoint + distance : endpoint - distance);
            }
        }
    }

    std::vector<std::vector<real_t>> panels;
    panels.reserve(cells.size());
    for (len_t j = 0; j < cells.size(); ++j) {
        auto cuts = sharedCuts;
        for (real_t beta : contacts[j]) {
            const auto contact = rootsAt(beta);
            cuts.push_back(contact.lower);
            cuts.push_back(contact.upper);
        }

        // A contact close to a primary bounce puts the inverse-square-root
        // layer outside the adjacent panel. Grade distances from the bounce
        // so that the mapped rule resolves this nearby layer.
        if (!constant && primarySquared < betaMax) {
            const auto orbit = rootsAt(orbitLevel);
            const auto reachable = rootsAt(reachableLevels[j]);
            cuts.push_back(reachable.lower);
            cuts.push_back(reachable.upper);

            // Measure both branches from the same set of original cuts.
            const auto originalCuts = cuts;
            for (bool leftBranch : {true, false}) {
                const real_t bounce = leftBranch ? orbit.lower : orbit.upper;
                const real_t branchWidth = std::abs(thetaMin - bounce);
                real_t nearestCut = branchWidth;
                for (real_t cut : originalCuts) {
                    const real_t distance = leftBranch ? cut - bounce : bounce - cut;
                    if (distance > 0)
                        nearestCut = std::min(nearestCut, distance);
                }
                if (nearestCut < 0.1 * branchWidth) {
                    for (real_t distance = 4 * nearestCut; distance > 0 && distance < branchWidth;
                        distance *= 4) {
                        cuts.push_back(leftBranch ? bounce + distance : bounce - distance);
                    }
                }
            }
        }

        panels.push_back(SplitPoints(rootsAt(reachableLevels[j]), std::move(cuts)));
    }
    return panels;
}

}  // namespace DREAM::KnockOnUtilities::Quadrature

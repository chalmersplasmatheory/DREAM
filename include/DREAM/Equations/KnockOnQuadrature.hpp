#ifndef _DREAM_EQUATIONS_KNOCK_ON_QUADRATURE_HPP
#define _DREAM_EQUATIONS_KNOCK_ON_QUADRATURE_HPP

#include <functional>
#include <memory>
#include <vector>

#include "FVM/config.h"

namespace DREAM::FVM { class FluxSurfaceAverager; }

namespace DREAM::KnockOnUtilities::Quadrature {

struct Interval { real_t lower, upper; };
struct WeightedPoint { real_t x, weight; };

// Include the endpoints, discard cuts outside the interval, and remove exact
// duplicates. Nearby distinct contacts must remain distinct.
std::vector<real_t> SplitPoints(Interval interval, std::vector<real_t> cuts);

/**
 * Candidate contact levels beta=1-Bmin/B for a fixed primary pitch u and
 * secondary face v. Extra roots from squaring only oversplit the integrals.
 * Includes face bounce and clipped-face events; endpoints are supplied by
 * the caller. Returned levels lie strictly in (0,min(u*u,betaMax)).
 */
std::vector<real_t> ContactFieldLevels(
    real_t u, real_t v, real_t xiStar, real_t betaMax
);

/**
 * Split a primary-pitch interval at contact/extremum, moving-bounce and
 * bounce-coincidence events, as well as -xiT, 0 and xiT. Supply all physical
 * and virtual secondary faces (including those introduced by trapped folding).
 * Both pitch signs are retained: superfluous cuts are harmless.
 */
std::vector<real_t> PrimaryPitchSplitPoints(
    Interval cell, real_t xiStar, real_t betaMax,
    const std::vector<real_t>& secondaryFaces
);

/**
 * A rule is constructed once and reused over many panels. Ordinary Gauss
 * nodes/weights are on [0,1]. MappedPoints uses x=sin^2(pi*t/2), with positive
 * weights normalized to preserve constants on each panel. N=1 is exactly
 * midpoint, including its weight. An interval with no representable interior
 * point returns no nodes; it cannot be sampled away from its singular ends.
 */
class GaussLegendreRule {
    std::vector<real_t> nodes, weights, mappedNodes, mappedWeights;
public:
    explicit GaussLegendreRule(len_t order);
    len_t GetOrder() const { return nodes.size(); }
    real_t GetNode(len_t i) const { return nodes.at(i); }
    real_t GetWeight(len_t i) const { return weights.at(i); }
    std::vector<WeightedPoint> MappedPoints(Interval panel) const;
};

struct FieldExtrema { real_t Bmin, Bmax, thetaMin, thetaMax; };
struct FieldLevel { real_t beta; Interval roots; };

/**
 * Prepared geometry for one radius. Uses the existing single-well topology:
 * one monotone branch on either side of Bmin. Angles are unwrapped onto
 * [thetaMax-2*pi,thetaMax]. Both roots are solved independently.
 *
 * Bracket samples accelerate inversion; they are NOT a replacement field
 * interpolant. All root iterations evaluate the original B callback. Rebuild
 * this object, or Reset it with current extrema, after the geometry changes.
 * No cache of floating-point pitch samples is retained between operations.
 */
class OrbitGeometry {
    struct Implementation;
    std::unique_ptr<Implementation> impl;
public:
    OrbitGeometry(FVM::FluxSurfaceAverager *fsa, len_t ir);
    // Also useful for analytic geometry and independent numerical tests.
    OrbitGeometry(
        FieldExtrema extrema, std::function<real_t(real_t)> evaluateB,
        const std::vector<real_t>& thetaSamples = {}
    );
    ~OrbitGeometry();
    OrbitGeometry(const OrbitGeometry&) = delete;
    OrbitGeometry& operator=(const OrbitGeometry&) = delete;

    void Reset(FieldExtrema extrema, const std::vector<real_t>& thetaSamples = {});
    bool IsConstantField() const;
    real_t GetBetaMax() const;
    Interval GetAngularDomain() const;
    Interval InvertBeta(real_t beta);

    // Sort/deduplicate exact repeated levels and solve each only once in this
    // batch. The returned vector is ordered by beta, not by angle.
    std::vector<FieldLevel> InvertLevels(std::vector<real_t> levels);

    // Physical orbit bounds, with no grid-index substitutions or endpoint
    // nudges. These do not perform DREAM's trapped-cell folding/remapping.
    Interval GetOrbitBounds(real_t u);
    Interval GetReachableBounds(real_t u, real_t lowerFace, real_t upperFace);
};

} // namespace DREAM::KnockOnUtilities::Quadrature
#endif

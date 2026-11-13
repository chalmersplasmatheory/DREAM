#ifndef _DREAM_EQUATIONS_KNOCK_ON_QUADRATURE_HPP
#define _DREAM_EQUATIONS_KNOCK_ON_QUADRATURE_HPP

#include <functional>
#include <memory>
#include <vector>

#include "FVM/config.h"

namespace DREAM::FVM { class FluxSurfaceAverager; }

namespace DREAM::KnockOnUtilities::Quadrature {

struct Interval {
    real_t lower;
    real_t upper;
};

struct WeightedPoint {
    real_t x;
    real_t weight;
};

// Keep both endpoints and all distinct interior cuts. Two nearby contacts
// can bound a narrow panel with nonzero support.
std::vector<real_t> SplitPoints(Interval interval, std::vector<real_t> cuts);

/**
 * Candidate contact levels beta = 1 - Bmin/B for a primary pitch u and an
 * outgoing pitch face. Includes the face bounce and contacts with a face
 * clipped to zero local pitch. The caller supplies the interval endpoints.
 * Returned levels lie strictly in (0, min(u*u, betaMax)).
 */
std::vector<real_t> ContactFieldLevels(
    real_t primaryPitch, real_t secondaryFace, real_t xiStar, real_t betaMax
);

/**
 * Split a primary-pitch interval at changes in the poloidal panel boundaries.
 * These include contacts at field extrema, contacts at bounces, coincident
 * bounces, and the trapped-pitch boundaries -xiT, 0 and xiT. Supply the
 * outgoing pitch faces, including any faces introduced by trapped folding.
 */
std::vector<real_t> PrimaryPitchSplitPoints(
    Interval cell, real_t xiStar, real_t betaMax,
    const std::vector<real_t> &secondaryFaces
);

/**
 * Store a Gauss rule on [0,1] and its endpoint-clustered counterpart.
 * MappedPoints uses x = sin^2(pi*t/2) within each panel and rescales the
 * weights to preserve constants. Order one gives the midpoint rule.
 * A panel with no representable interior point yields no nodes.
 */
class GaussLegendreRule {
    std::vector<real_t> nodes;
    std::vector<real_t> weights;
    std::vector<real_t> mappedNodes;
    std::vector<real_t> mappedWeights;

   public:
    explicit GaussLegendreRule(len_t order);
    len_t GetOrder() const { return nodes.size(); }
    real_t GetNode(len_t i) const { return nodes.at(i); }
    real_t GetWeight(len_t i) const { return weights.at(i); }
    std::vector<WeightedPoint> MappedPoints(Interval panel) const;
};

struct FieldExtrema {
    real_t Bmin;
    real_t Bmax;
    real_t thetaMin;
    real_t thetaMax;
};

struct FieldLevel {
    real_t beta;
    Interval roots;
};

/**
 * Prepare one flux surface with a monotone branch on either side of Bmin.
 * Angles are unwrapped onto [thetaMax - 2*pi, thetaMax].
 *
 * Bracket samples locate the interval for each root. The root solver
 * evaluates the original field on each branch. Rebuild this object or
 * call Reset() with current extrema after the geometry changes.
 */
class OrbitGeometry {
    struct Implementation;
    std::unique_ptr<Implementation> impl;

   public:
    OrbitGeometry(FVM::FluxSurfaceAverager *fsa, len_t ir);
    // Also useful for analytic geometry and independent numerical tests.
    OrbitGeometry(
        FieldExtrema extrema, std::function<real_t(real_t)> evaluateB,
        const std::vector<real_t> &thetaSamples = {}
    );
    ~OrbitGeometry();
    OrbitGeometry(const OrbitGeometry &) = delete;
    OrbitGeometry &operator=(const OrbitGeometry &) = delete;

    void Reset(FieldExtrema extrema, const std::vector<real_t> &thetaSamples = {});
    bool IsConstantField() const;
    real_t GetBetaMax() const;
    Interval GetAngularDomain() const;
    Interval InvertBeta(real_t beta);

    // Solve each distinct level once and return the roots in order of beta.
    std::vector<FieldLevel> InvertLevels(std::vector<real_t> levels);

    // Physical bounds of a primary orbit and of its overlap with an outgoing
    // pitch interval. Grid-aware trapped folding is done by the caller.
    Interval GetOrbitBounds(real_t primaryPitch);
    Interval GetReachableBounds(real_t primaryPitch, real_t lowerFace, real_t upperFace);
};

/**
 * Build poloidal panels for several outgoing pitch intervals at fixed
 * primary pitch and scattering cosine. Panel boundaries follow contacts,
 * the reachable interval, the field minimum, and separatrix or bounce layers.
 * Primary-pitch averaging samples on either side of the separatrix, where
 * the pointwise orbit integral is finite.
 */
std::vector<std::vector<real_t>> PoloidalSplitPoints(
    OrbitGeometry &geometry, real_t primaryPitch, real_t xiStar,
    const std::vector<Interval> &cells
);

}  // namespace DREAM::KnockOnUtilities::Quadrature
#endif /*_DREAM_EQUATIONS_KNOCK_ON_QUADRATURE_HPP*/

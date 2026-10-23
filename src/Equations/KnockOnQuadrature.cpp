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

using namespace DREAM;
namespace Q = DREAM::KnockOnUtilities::Quadrature;

namespace {
void ValidatePitch(real_t x) {
    if (!std::isfinite(x) || x < -1 || x > 1)
        throw DREAMException("KnockOnQuadrature: pitch/cosine must be finite and in [-1,1].");
}
void ValidateBetaMax(real_t betaMax) {
    if (!std::isfinite(betaMax) || betaMax < 0 || betaMax >= 1)
        throw DREAMException("KnockOnQuadrature: betaMax must be finite and in [0,1).");
}
void ValidateInterval(Q::Interval x) {
    if (!std::isfinite(x.lower) || !std::isfinite(x.upper) || x.lower > x.upper
        || !std::isfinite(x.upper-x.lower))
        throw DREAMException("KnockOnQuadrature: invalid integration interval.");
}
void SortUnique(std::vector<real_t>& x) {
    std::sort(x.begin(), x.end());
    x.erase(std::unique(x.begin(), x.end()), x.end());
}
}

std::vector<real_t> Q::SplitPoints(Interval interval, std::vector<real_t> cuts) {
    ValidateInterval(interval);
    for (real_t x : cuts)
        if (!std::isfinite(x))
            throw DREAMException("KnockOnQuadrature: non-finite split point.");
    cuts.erase(std::remove_if(cuts.begin(), cuts.end(), [&](real_t x) {
        return x <= interval.lower || x >= interval.upper;
    }), cuts.end());
    cuts.push_back(interval.lower);
    cuts.push_back(interval.upper);
    SortUnique(cuts);
    return cuts;
}

std::vector<real_t> Q::ContactFieldLevels(
    real_t u, real_t v, real_t xiStar, real_t betaMax
) {
    ValidatePitch(u); ValidatePitch(v); ValidatePitch(xiStar); ValidateBetaMax(betaMax);
    // Long-double intermediates protect cancellation near xiStar=1 and u=v.
    const long double U=u, V=v, c=std::abs((long double)xiStar);
    const long double u2=U*U, v2=V*V, s2=(1-c)*(1+c);
    const long double reach=std::min(u2, (long double)betaMax);
    std::vector<real_t> levels;
    const auto add = [&](long double beta) {
        const real_t b=static_cast<real_t>(beta);
        // Check again after conversion: rounding can put a root on an endpoint.
        if (beta>0 && beta<reach && b>0 && b<std::min(u*u,betaMax))
            levels.push_back(b);
    };
    add(v2); // a secondary face becomes reachable
    if (c < 1) {
        const long double X=std::sqrt((1-U)*(1+U));
        const long double Y=std::sqrt((1-V)*(1+V));
        for (long double den : {1+c, 1-c}) {
            const long double beta=1-(X-Y)*(X-Y)/s2-2*X*Y/den;
            if (beta <= std::min(u2,v2))
                add(beta);
        }
    }
    // Contact with a face clipped to local pitch zero. Compare before dividing
    // to avoid overflow for arbitrarily small nonzero scattering cosine.
    const long double numerator=u2-s2;
    if (c>0 && numerator>0 && numerator<c*c*reach && numerator>=c*c*v2)
        add(numerator/(c*c));
    SortUnique(levels);
    return levels;
}

std::vector<real_t> Q::PrimaryPitchSplitPoints(
    Interval cell, real_t xiStar, real_t betaMax,
    const std::vector<real_t>& secondaryFaces
) {
    ValidateInterval(cell); ValidatePitch(cell.lower); ValidatePitch(cell.upper);
    ValidatePitch(xiStar); ValidateBetaMax(betaMax);
    const long double c=std::abs((long double)xiStar), s2=(1-c)*(1+c);
    const long double s=std::sqrt(s2);
    const real_t T=std::sqrt(betaMax);
    std::vector<real_t> cuts{-T,0,T};
    const auto addSquared = [&](long double a2) {
        if (a2>=0 && a2<=1) {
            const real_t a=static_cast<real_t>(std::sqrt(a2));
            cuts.push_back(-a); cuts.push_back(a);
        }
    };
    for (real_t v : secondaryFaces) {
        ValidatePitch(v);
        const long double V=v, v2=V*V;
        for (long double beta : {0.L, (long double)betaMax}) {
            if (v2<beta) continue;
            const long double t=std::copysign(std::sqrt((v2-beta)/(1-beta)), V);
            const long double d=s*std::sqrt(std::max(0.L,(1-t)*(1+t)));
            for (long double eta : {c*t-d,c*t+d})
                addSquared(beta+(1-beta)*eta*eta);
        }
        if (v2<=betaMax) {
            addSquared(s2+c*c*v2); // contact reaches the secondary-face bounce
            addSquared(v2);       // primary and secondary bounces coincide
        }
        const long double numerator=v2-s2;
        if (c>0 && numerator>=0 && numerator<=c*c*betaMax)
            addSquared(numerator/(c*c)); // contact reaches the primary bounce
    }
    return SplitPoints(cell,std::move(cuts));
}

Q::GaussLegendreRule::GaussLegendreRule(len_t order) {
    if (order<1)
        throw DREAMException("KnockOnQuadrature: Gauss order must be at least 1.");
    std::unique_ptr<gsl_integration_glfixed_table, decltype(&gsl_integration_glfixed_table_free)>
        table(gsl_integration_glfixed_table_alloc(order), &gsl_integration_glfixed_table_free);
    if (!table)
        throw DREAMException("KnockOnQuadrature: cannot allocate Gauss rule.");
    nodes.resize(order); weights.resize(order);
    mappedNodes.resize(order); mappedWeights.resize(order);
    real_t sum=0;
    for (len_t i=0; i<order; ++i) {
        if (gsl_integration_glfixed_point(0,1,i,&nodes[i],&weights[i],table.get()) != GSL_SUCCESS)
            throw DREAMException("KnockOnQuadrature: cannot construct Gauss rule.");
        const real_t t=std::min(nodes[i],1-nodes[i]);
        const real_t y=std::pow(std::sin(M_PI*t/2),2);
        mappedNodes[i]=nodes[i]<0.5 ? y : 1-y;
        mappedWeights[i]=weights[i]*(M_PI/2)*std::sin(M_PI*nodes[i]);
        sum+=mappedWeights[i];
    }
    for (real_t& w : mappedWeights) w/=sum;
    if (order%2) mappedNodes[order/2]=0.5;
}

std::vector<Q::WeightedPoint> Q::GaussLegendreRule::MappedPoints(Interval panel) const {
    ValidateInterval(panel);
    std::vector<WeightedPoint> out;
    const real_t first=std::nextafter(panel.lower,panel.upper);
    const real_t last=std::nextafter(panel.upper,panel.lower);
    if (panel.lower==panel.upper || first>=panel.upper) return out;
    out.reserve(GetOrder());
    const real_t width=panel.upper-panel.lower;
    for (len_t i=0; i<GetOrder(); ++i) {
        const real_t y=mappedNodes[i];
        const real_t x=y<=0.5 ? panel.lower+width*y : panel.upper-width*(1-y);
        out.push_back({std::max(first,std::min(last,x)), width*mappedWeights[i]});
    }
    return out;
}

struct Q::OrbitGeometry::Implementation {
    struct Branch { std::vector<real_t> theta, B; };
    std::function<real_t(real_t)> evaluateB;
    std::array<Branch,2> branches;
    FieldExtrema extrema;
    bool constant=false;
    real_t betaMax=0;
    std::unique_ptr<gsl_root_fsolver, decltype(&gsl_root_fsolver_free)>
        solver{nullptr,&gsl_root_fsolver_free};

    real_t Field(real_t theta) const {
        const real_t B=evaluateB(theta);
        if (!std::isfinite(B) || B<0)
            throw DREAMException("KnockOnQuadrature: invalid magnetic field evaluation.");
        return B;
    }
    struct RootParameters { Implementation *geometry; real_t target; };
    static double RootFunction(double theta, void *p) {
        const auto& params=*static_cast<RootParameters*>(p);
        return params.geometry->Field(theta)-params.target;
    }
    real_t Root(Branch& branch, real_t target) {
        // Levels closer to an extremum than the field evaluation can resolve
        // have that extremum as their only representable field bracket.
        if (target<=branch.B.front()) return branch.theta.front();
        if (target>=branch.B.back()) return branch.theta.back();
        const auto pos=std::lower_bound(branch.B.begin(),branch.B.end(),target);
        const std::size_t i=std::max<std::size_t>(1,pos-branch.B.begin());
        real_t a=branch.theta[i-1], b=branch.theta[i];
        if (a>b) std::swap(a,b);
        const real_t fa=Field(a)-target, fb=Field(b)-target;
        if (fa==0) return a;
        if (fb==0) return b;
        if ((fa<0)==(fb<0))
            throw DREAMException("KnockOnQuadrature: field level is not bracketed.");
        RootParameters params{this,target};
        gsl_function function;
        function.function=&RootFunction; function.params=&params;
        if (gsl_root_fsolver_set(solver.get(),&function,a,b) != GSL_SUCCESS)
            throw DREAMException("KnockOnQuadrature: cannot initialize field inversion.");
        for (len_t iteration=0; iteration<100; ++iteration) {
            if (gsl_root_fsolver_iterate(solver.get()) != GSL_SUCCESS)
                throw DREAMException("KnockOnQuadrature: field inversion failed.");
            a=gsl_root_fsolver_x_lower(solver.get());
            b=gsl_root_fsolver_x_upper(solver.get());
            if (gsl_root_test_interval(a,b,1e-11,1e-11)==GSL_SUCCESS)
                return gsl_root_fsolver_root(solver.get());
        }
        throw DREAMException("KnockOnQuadrature: field inversion did not converge.");
    }
};

Q::OrbitGeometry::OrbitGeometry(
    FieldExtrema extrema, std::function<real_t(real_t)> evaluateB,
    const std::vector<real_t>& thetaSamples
) : impl(new Implementation()) {
    if (!evaluateB)
        throw DREAMException("KnockOnQuadrature: magnetic field callback is empty.");
    impl->evaluateB=std::move(evaluateB);
    impl->solver.reset(gsl_root_fsolver_alloc(gsl_root_fsolver_brent));
    if (!impl->solver)
        throw DREAMException("KnockOnQuadrature: cannot allocate field root solver.");
    Reset(extrema,thetaSamples);
}

namespace {
Q::FieldExtrema GetExtrema(FVM::FluxSurfaceAverager *fsa, len_t ir) {
    if (!fsa) throw DREAMException("KnockOnQuadrature: null flux surface averager.");
    Q::FieldExtrema extrema;
    extrema.Bmin=fsa->GetBmin(ir,FVM::FLUXGRIDTYPE_DISTRIBUTION,&extrema.thetaMin);
    extrema.Bmax=fsa->GetBmax(ir,FVM::FLUXGRIDTYPE_DISTRIBUTION,&extrema.thetaMax);
    return extrema;
}
std::vector<real_t> GetThetaSamples(FVM::FluxSurfaceAverager *fsa) {
    if (!fsa) throw DREAMException("KnockOnQuadrature: null flux surface averager.");
    const len_t n=fsa->GetNTheta();
    if (n==0) return {};
    if (!fsa->GetTheta())
        throw DREAMException("KnockOnQuadrature: missing bracket sampling angles.");
    return {fsa->GetTheta(),fsa->GetTheta()+n};
}
}

Q::OrbitGeometry::OrbitGeometry(FVM::FluxSurfaceAverager *fsa, len_t ir)
    : OrbitGeometry(GetExtrema(fsa,ir),[fsa,ir](real_t theta) {
        return fsa->BAtTheta(ir,theta,FVM::FLUXGRIDTYPE_DISTRIBUTION);
    },GetThetaSamples(fsa)) {}

Q::OrbitGeometry::~OrbitGeometry() = default;

void Q::OrbitGeometry::Reset(FieldExtrema extrema, const std::vector<real_t>& thetaSamples) {
    if (!std::isfinite(extrema.Bmin) || !std::isfinite(extrema.Bmax)
        || !std::isfinite(extrema.thetaMin) || !std::isfinite(extrema.thetaMax)
        || extrema.Bmin<0 || extrema.Bmax<extrema.Bmin)
        throw DREAMException("KnockOnQuadrature: invalid field extrema.");
    const real_t tolerance=100*std::numeric_limits<real_t>::epsilon()*(1+extrema.Bmax);
    const bool constant=extrema.Bmax-extrema.Bmin<=tolerance;
    if (!constant && extrema.Bmin==0)
        throw DREAMException("KnockOnQuadrature: varying field requires positive Bmin.");
    const real_t period=2*M_PI;
    real_t displacement=std::fmod(extrema.thetaMax-extrema.thetaMin,period);
    if (displacement<0) displacement+=period;
    extrema.thetaMin=extrema.thetaMax-displacement;
    if (!constant && !(displacement>0 && displacement<period))
        throw DREAMException("KnockOnQuadrature: distinct field extrema need distinct angles.");
    const real_t betaMax=constant ? 0 : (extrema.Bmax-extrema.Bmin)/extrema.Bmax;
    ValidateBetaMax(betaMax);
    std::array<Implementation::Branch,2> branches;
    if (constant) {
        impl->extrema=extrema;
        impl->constant=true;
        impl->betaMax=0;
        impl->branches=std::move(branches);
        return;
    }

    const real_t left=extrema.thetaMax-period, middle=extrema.thetaMin, right=extrema.thetaMax;
    // Sampling only supplies brackets. A small internal seed covers both
    // branches even when FSA stores quadrature nodes on only half the surface.
    std::vector<real_t> samples=thetaSamples;
    for (len_t i=0; i<=32; ++i) {
        samples.push_back(left+(middle-left)*i/32);
        samples.push_back(middle+(right-middle)*i/32);
    }
    std::array<std::vector<real_t>,2> angles;
    for (real_t theta : samples) {
        if (!std::isfinite(theta))
            throw DREAMException("KnockOnQuadrature: non-finite bracket sample.");
        real_t x=std::fmod(theta-left,period);
        if (x<0) x+=period;
        x+=left;
        if (x>left && x<middle) angles[0].push_back(x);
        else if (x>middle && x<right) angles[1].push_back(x);
    }
    for (len_t side=0; side<2; ++side) {
        angles[side].push_back(middle);
        angles[side].push_back(side==0 ? left : right);
        SortUnique(angles[side]);
        if (side==0) std::reverse(angles[side].begin(),angles[side].end());
        auto& branch=branches[side];
        for (real_t theta : angles[side]) {
            const real_t B=impl->Field(theta);
            if (B<extrema.Bmin-tolerance || B>extrema.Bmax+tolerance)
                throw DREAMException("KnockOnQuadrature: field disagrees with its extrema.");
            if (!branch.B.empty()) {
                if (B<branch.B.back()-tolerance)
                    throw DREAMException("KnockOnQuadrature: expected a monotone magnetic branch.");
                // Tiny decreases from roundoff must not enter lower_bound.
                if (B<=branch.B.back()) continue;
            }
            branch.theta.push_back(theta); branch.B.push_back(B);
        }
        if (branch.B.size()<2)
            throw DREAMException("KnockOnQuadrature: degenerate magnetic branch.");
    }
    // Commit only fully prepared brackets, leaving the previous context usable
    // if validation or field evaluation failed during Reset().
    impl->extrema=extrema;
    impl->constant=false;
    impl->betaMax=betaMax;
    impl->branches=std::move(branches);
}

bool Q::OrbitGeometry::IsConstantField() const { return impl->constant; }
real_t Q::OrbitGeometry::GetBetaMax() const { return impl->betaMax; }
Q::Interval Q::OrbitGeometry::GetAngularDomain() const {
    return {impl->extrema.thetaMax-2*M_PI,impl->extrema.thetaMax};
}

Q::Interval Q::OrbitGeometry::InvertBeta(real_t beta) {
    if (!std::isfinite(beta) || beta<0 || beta>impl->betaMax)
        throw DREAMException("KnockOnQuadrature: field level is outside this surface.");
    if (impl->constant) return GetAngularDomain();
    const real_t minimum=impl->extrema.thetaMin;
    if (beta==0) return {minimum,minimum};
    if (beta==impl->betaMax) return GetAngularDomain();
    const real_t target=impl->extrema.Bmin/(1-beta);
    return {impl->Root(impl->branches[0],target),impl->Root(impl->branches[1],target)};
}

std::vector<Q::FieldLevel> Q::OrbitGeometry::InvertLevels(std::vector<real_t> levels) {
    for (real_t beta : levels)
        if (!std::isfinite(beta) || beta<0 || beta>impl->betaMax)
            throw DREAMException("KnockOnQuadrature: field level is outside this surface.");
    SortUnique(levels);
    std::vector<FieldLevel> out;
    out.reserve(levels.size());
    for (real_t beta : levels) out.push_back({beta,InvertBeta(beta)});
    return out;
}

Q::Interval Q::OrbitGeometry::GetOrbitBounds(real_t u) {
    ValidatePitch(u);
    if (impl->constant || u*u>=impl->betaMax) return GetAngularDomain();
    return InvertBeta(u*u);
}

Q::Interval Q::OrbitGeometry::GetReachableBounds(real_t u, real_t lowerFace, real_t upperFace) {
    ValidatePitch(u); ValidatePitch(lowerFace); ValidatePitch(upperFace);
    if (lowerFace>upperFace)
        throw DREAMException("KnockOnQuadrature: secondary cell faces are reversed.");
    const real_t reach=std::min(u*u,std::max(lowerFace*lowerFace,upperFace*upperFace));
    if (impl->constant || reach>=impl->betaMax) return GetAngularDomain();
    return InvertBeta(reach);
}

std::vector<std::vector<real_t>> Q::PoloidalSplitPoints(
    OrbitGeometry& geometry, real_t u, real_t xiStar,
    const std::vector<Interval>& cells
) {
    ValidatePitch(u); ValidatePitch(xiStar);
    const real_t betaMax=geometry.GetBetaMax(), u2=u*u;
    const bool constant=geometry.IsConstantField();
    if (!constant && (u2==betaMax || std::abs(u)==std::sqrt(betaMax)))
        throw DREAMException("KnockOnQuadrature: pointwise primary pitch is on the separatrix.");

    // Keep each row's contacts separate: applying every column contact to
    // every row would needlessly increase the number of integrand samples.
    std::vector<std::vector<real_t>> contacts(cells.size());
    const real_t orbitLevel=std::min(u2,betaMax);
    std::vector<real_t> reach(cells.size()), levels{0, betaMax, orbitLevel};
    for (len_t j=0; j<cells.size(); ++j) {
        ValidatePitch(cells[j].lower); ValidatePitch(cells[j].upper);
        ValidateInterval(cells[j]);
        reach[j]=std::min({u2, betaMax,
            std::max(cells[j].lower*cells[j].lower,cells[j].upper*cells[j].upper)});
        levels.push_back(reach[j]);
        for (real_t face : {cells[j].lower,cells[j].upper}) {
            const auto candidates=ContactFieldLevels(u,face,xiStar,betaMax);
            for (real_t beta : candidates)
                if (beta<reach[j]) contacts[j].push_back(beta);
        }
        levels.insert(levels.end(),contacts[j].begin(),contacts[j].end());
    }

    // At a smooth field maximum, betaMax-beta ~ curvature*distance^2/2.
    // Invert a field drop |u^2-betaMax| to determine the layer width on
    // each branch directly. This also works for asymmetric/flat maxima.
    const real_t gap=std::abs(u2-betaMax);
    const bool nearSeparatrix=!constant && gap>0 && gap<0.1*betaMax;
    const real_t layerLevel=nearSeparatrix ? betaMax-gap : 0;
    if (nearSeparatrix) levels.push_back(layerLevel);
    const auto roots=geometry.InvertLevels(std::move(levels));
    auto at=[&](real_t beta) {
        return std::lower_bound(roots.begin(),roots.end(),beta,
            [](const FieldLevel& level, real_t value) { return level.beta<value; })->roots;
    };
    const auto domain=geometry.GetAngularDomain();
    const real_t minimum=constant ? (domain.lower+domain.upper)/2 : at(0).lower;
    std::vector<real_t> sharedCuts{minimum};
    if (nearSeparatrix) {
        const auto layer=at(layerLevel);
        for (bool left : {true,false}) {
            const real_t endpoint=left ? domain.lower : domain.upper;
            const real_t width=left ? layer.lower-endpoint : endpoint-layer.upper;
            const real_t branch=std::abs(minimum-endpoint);
            // A field level indistinguishable from Bmax has zero width.
            // Do not loop forever or invent a finite layer in that case.
            for (real_t d=width; d>0 && d<branch; d*=4)
                sharedCuts.push_back(left ? endpoint+d : endpoint-d);
        }
    }

    std::vector<std::vector<real_t>> panels;
    panels.reserve(cells.size());
    for (len_t j=0; j<cells.size(); ++j) {
        auto cuts=sharedCuts;
        for (real_t beta : contacts[j]) {
            const auto contact=at(beta);
            cuts.push_back(contact.lower); cuts.push_back(contact.upper);
        }
        // A contact close to a primary bounce puts the inverse-square-root
        // layer just OUTSIDE the adjacent panel. Its own endpoint map cannot
        // remove that nearby singularity. Grade distances from the physical
        // bounce by factors of four before applying the panel map. Only the
        // affected row is refined; all added angles reuse the existing roots.
        if (!constant && u2<betaMax) {
            const auto orbit=at(orbitLevel), bounds=at(reach[j]);
            cuts.push_back(bounds.lower); cuts.push_back(bounds.upper);
            const auto originalCuts=cuts;
            for (bool left : {true,false}) {
                const real_t endpoint=left ? orbit.lower : orbit.upper;
                const real_t branch=std::abs(minimum-endpoint);
                real_t distance=branch;
                for (real_t cut : originalCuts) {
                    const real_t d=left ? cut-endpoint : endpoint-cut;
                    if (d>0) distance=std::min(distance,d);
                }
                if (distance<0.1*branch)
                    for (real_t d=4*distance; d>0 && d<branch; d*=4)
                        cuts.push_back(left ? endpoint+d : endpoint-d);
            }
        }
        panels.push_back(SplitPoints(at(reach[j]),std::move(cuts)));
    }
    return panels;
}

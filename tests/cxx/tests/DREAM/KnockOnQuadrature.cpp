#include "KnockOnQuadrature.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "DREAM/DREAMException.hpp"
#include "DREAM/Equations/KnockOnQuadrature.hpp"

namespace Q = DREAM::KnockOnUtilities::Quadrature;
using DREAMTESTS::_DREAM::KnockOnQuadrature;

namespace {
bool Contains(const std::vector<real_t>& x, real_t value, real_t tolerance=2e-12) {
    for (real_t y : x) if (std::abs(value-y)<=tolerance) return true;
    return false;
}
template<typename F> bool Throws(F&& f) {
    try { f(); } catch (const DREAM::DREAMException&) { return true; }
    return false;
}
// Independent local geometry: start with a primary local pitch and the two
// edges of its outgoing scattering cone, then map them to invariant pitch.
real_t InvariantPitch(real_t local, real_t beta) {
    return std::copysign(std::sqrt(beta+(1-beta)*local*local),local);
}
real_t LocalFace(real_t v, real_t beta) {
    return std::copysign(std::sqrt(std::max(0.0,(v*v-beta)/(1-beta))),v);
}
// Probability in a local secondary cell, using the gyro angle directly.
// It is used only to check support topology, independently of contact algebra.
real_t Fraction(real_t u, real_t lo, real_t hi, real_t chi, real_t beta) {
    const real_t eta=LocalFace(u,beta), lower=LocalFace(lo,beta), upper=LocalFace(hi,beta);
    const real_t a=eta*chi, b=std::sqrt((1-eta*eta)*(1-chi*chi));
    if (b==0) return lower<a && a<=upper ? 1 : 0;
    const auto acosClipped=[](real_t x) { return std::acos(std::max(-1.0,std::min(1.0,x))); };
    return (acosClipped((lower-a)/b)-acosClipped((upper-a)/b))/M_PI;
}
}

bool KnockOnQuadrature::Run(bool) {
    bool success=true;
    const auto check=[&](bool result, const char *description) {
        if (result) PrintOK(description);
        else { PrintError(description); success=false; }
    };
    check(CheckContactLevels(),"Contact splits resolve independently constructed scattering/support events.");
    check(CheckPrimaryPitchSplits(),"Primary-pitch splits include endpoint contacts and bounce coincidences.");
    check(CheckGeometry(),"Prepared B-only inversion handles both branches, explicit pitches and reset.");
    check(CheckGaussRule(),"Gauss rules and mapped panels preserve constants and midpoint weights.");
    check(CheckValidation(),"Invalid quadrature and geometry inputs are rejected.");
    return success;
}

bool KnockOnQuadrature::CheckContactLevels() {
    const real_t betaMax=0.55;
    // Construct actual contacts without solving the squared contact equation.
    for (real_t beta : {0.001,0.05,0.2,0.5})
        for (real_t eta : {-0.95,-0.4,0.1,0.7})
            for (real_t chi : {-0.99,0.0,0.3,0.9,0.9999}) {
                const real_t u=InvariantPitch(eta,beta);
                const real_t d=std::sqrt((1-eta*eta)*(1-chi*chi));
                for (real_t local : {chi*eta-d,chi*eta+d}) {
                    const real_t v=InvariantPitch(local,beta);
                    const auto levels=Q::ContactFieldLevels(u,v,chi,betaMax);
                    if (!Contains(levels,beta,3e-11)) {
                        PrintError("Missing contact: u=%.17g v=%.17g chi=%.17g beta=%.17g",u,v,chi,beta);
                        return false;
                    }
                }
            }
    // Across each panel the local probability must stay zero, partial or full:
    // a support transition inside a panel would mean an omitted contact.
    const std::vector<real_t> faces{-1,-0.9,-0.6,-0.2,0,0.15,0.65,0.95,1};
    for (real_t u : {-1.0,-0.95,-0.7,-0.1,0.2,0.65,0.99,1.0})
        for (real_t chi : {0.0,1e-8,0.1,0.55,0.99,1-1e-10,1.0})
            for (std::size_t j=0; j+1<faces.size(); ++j) {
                const real_t lo=faces[j], hi=faces[j+1];
                const real_t reach=std::min({betaMax,u*u,std::max(lo*lo,hi*hi)});
                auto levels=Q::ContactFieldLevels(u,lo,chi,betaMax);
                const auto other=Q::ContactFieldLevels(u,hi,chi,betaMax);
                levels.insert(levels.end(),other.begin(),other.end());
                levels=Q::SplitPoints({0,reach},std::move(levels));
                for (std::size_t i=1; i<levels.size(); ++i) {
                    int previous=-1;
                    for (int k=1; k<=9; ++k) {
                        const real_t beta=levels[i-1]+(levels[i]-levels[i-1])*k/10;
                        if (!(beta>levels[i-1] && beta<levels[i])) continue;
                        const real_t f=Fraction(u,lo,hi,chi,beta);
                        const int state=f<=0 ? 0 : (f>=1 ? 2 : 1);
                        if (previous>=0 && state!=previous) {
                            PrintError("Support changes inside panel: u=%.17g chi=%.17g cell=[%.17g,%.17g] beta=[%.17g,%.17g]",u,chi,lo,hi,levels[i-1],levels[i]);
                            return false;
                        }
                        previous=state;
                    }
                }
            }
    return true;
}

bool KnockOnQuadrature::CheckPrimaryPitchSplits() {
    const real_t betaMax=0.46, T=std::sqrt(betaMax);
    for (real_t chi : {0.0,1e-10,0.3,0.55,0.99,1-1e-10,1.0}) {
        const std::vector<real_t> faces{-1,-T,-0.41,0,T,0.91,1};
        const auto cuts=Q::PrimaryPitchSplitPoints({-1,1},chi,betaMax,faces);
        for (real_t x : {-1.0,-T,-0.41,0.0,0.41,T,1.0})
            if (!Contains(cuts,x)) return false;
        for (real_t beta : {0.0,betaMax})
            for (real_t v : faces) {
                if (v*v<beta) continue;
                const real_t t=LocalFace(v,beta);
                const real_t d=std::sqrt((1-t*t)*(1-chi*chi));
                for (real_t local : {chi*t-d,chi*t+d}) {
                    const real_t a=std::abs(InvariantPitch(local,beta));
                    if (!Contains(cuts,a) || !Contains(cuts,-a)) return false;
                }
            }
    }
    // In a cylinder the contact condition is u=chi*v +/- sin(theta*)*sqrt(1-v^2).
    const real_t v=0.3, chi=0.8, d=std::sqrt((1-chi*chi)*(1-v*v));
    const auto cylindrical=Q::PrimaryPitchSplitPoints({-1,1},chi,0,{v});
    if (!Contains(cylindrical,chi*v-d) || !Contains(cylindrical,chi*v+d)) return false;
    // Do not merge two genuine nearby contacts with an absolute tolerance.
    const real_t a=0.2, b=std::nextafter(a,1.0);
    const auto narrow=Q::SplitPoints({0,1},{a,b,a,-1,2});
    if (narrow.size()!=4 || narrow[1]!=a || narrow[2]!=b) return false;
    return true;
}

bool KnockOnQuadrature::CheckGeometry() {
    for (real_t epsilon : {0.05,0.3,0.7})
        for (real_t shift : {-1.3,0.0,5.4})
            for (real_t asymmetry : {0.0,0.25}) {
                const real_t Bmax=(1+epsilon)/(1-epsilon), betaMax=2*epsilon/(1+epsilon);
                std::size_t calls=0;
                const auto field=[&](real_t theta) {
                    ++calls;
                    const real_t x=theta-shift, phase=x+asymmetry*std::pow(std::sin(x),2);
                    return (1+epsilon)/(1+epsilon*std::cos(phase));
                };
                Q::OrbitGeometry geometry({1,Bmax,shift,shift+M_PI},field);
                for (real_t fraction : {1e-8,0.01,0.3,0.9,1-1e-8}) {
                    const real_t beta=fraction*geometry.GetBetaMax();
                    const auto roots=geometry.InvertBeta(beta);
                    if (!(roots.lower<shift && roots.upper>shift)) return false;
                    for (real_t theta : {roots.lower,roots.upper}) {
                        if (std::abs(1-1/field(theta)-beta)>2e-10) return false;
                        const real_t x=theta-shift, phase=x+asymmetry*std::pow(std::sin(x),2);
                        const real_t exact=2*std::asin(std::sqrt(beta/betaMax));
                        if (std::abs(std::abs(phase)-exact)>2e-8) return false;
                    }
                    // Bounds use the actual sampled u, independent of grid indices.
                    const real_t u=std::sqrt(beta);
                    const auto orbit=geometry.GetOrbitBounds(u);
                    if (std::abs(orbit.lower-roots.lower)>2e-10
                        || std::abs(orbit.upper-roots.upper)>2e-10) return false;
                }
                const auto full=geometry.GetOrbitBounds(1);
                if (std::abs(full.lower-(shift-M_PI))>1e-14
                    || std::abs(full.upper-(shift+M_PI))>1e-14) return false;
                const real_t beta=0.3*geometry.GetBetaMax();
                calls=0;
                geometry.InvertLevels({beta});
                const auto uniqueCalls=calls;
                calls=0;
                const auto batch=geometry.InvertLevels({beta,beta,beta});
                if (batch.size()!=1 || calls!=uniqueCalls) return false;
                const auto bounded=geometry.GetReachableBounds(0.9,0.1,0.2);
                const auto expected=geometry.InvertBeta(0.2*0.2);
                if (bounded.lower!=expected.lower || bounded.upper!=expected.upper) return false;
            }
    // Adjacent floating-point levels around the extrema and separatrix must
    // remain finite even when the field cannot resolve distinct root angles.
    {
        const auto field=[](real_t theta) { return 1.3/(1+0.3*std::cos(theta)); };
        Q::OrbitGeometry geometry({1,1.3/0.7,0,M_PI},field);
        const real_t betaMax=geometry.GetBetaMax(), T=std::sqrt(betaMax);
        for (real_t beta : {0.0,std::nextafter(0.0,1.0),std::nextafter(betaMax,0.0),betaMax}) {
            const auto roots=geometry.InvertBeta(beta);
            if (!std::isfinite(roots.lower) || !std::isfinite(roots.upper)
                || roots.lower>roots.upper) return false;
        }
        for (real_t u : {std::nextafter(T,0.0),T,std::nextafter(T,1.0)}) {
            const auto bounds=geometry.GetOrbitBounds(u);
            if (!std::isfinite(bounds.lower) || !std::isfinite(bounds.upper)) return false;
        }
    }
    // Exercise reset: no roots/brackets from the previous field survive.
    real_t epsilon=0.3;
    const auto field=[&](real_t theta) { return (1+epsilon)/(1+epsilon*std::cos(theta)); };
    Q::OrbitGeometry geometry({1,(1+epsilon)/(1-epsilon),0,M_PI},field);
    geometry.InvertBeta(0.2);
    epsilon=0.05;
    geometry.Reset({1,(1+epsilon)/(1-epsilon),0,M_PI});
    if (std::abs(geometry.GetBetaMax()-2*epsilon/(1+epsilon))>1e-15) return false;
    const auto root=geometry.InvertBeta(0.04);
    if (std::abs(1-1/field(root.upper)-0.04)>1e-10) return false;
    for (real_t B : {0.0,1e-6,2.0}) {
        Q::OrbitGeometry constant({B,B,0,M_PI},[B](real_t) { return B; });
        if (!constant.IsConstantField() || constant.GetBetaMax()!=0) return false;
        const auto bounds=constant.GetOrbitBounds(0);
        if (bounds.lower!=-M_PI || bounds.upper!=M_PI) return false;
    }
    // Extrema supplied on opposite sides of the periodic cut must describe
    // the same physical branches as unwrapped extrema.
    const real_t shift=5.4;
    const auto shifted=[&](real_t theta) { return 1.3/(1+0.3*std::cos(theta-shift)); };
    Q::OrbitGeometry wrapped({1,1.3/0.7,shift,std::fmod(shift+M_PI,2*M_PI)},shifted);
    const auto roots=wrapped.InvertBeta(0.2);
    if (std::abs(1-1/shifted(roots.lower)-0.2)>1e-10
        || std::abs(1-1/shifted(roots.upper)-0.2)>1e-10) return false;
    // An additional period in an extremum angle should not change the answer.
    wrapped.Reset({1,1.3/0.7,shift+2*M_PI,std::fmod(shift+M_PI,2*M_PI)});
    const auto repeated=wrapped.InvertBeta(0.2);
    if (std::abs(repeated.lower-roots.lower)>1e-10
        || std::abs(repeated.upper-roots.upper)>1e-10) return false;
    return true;
}

bool KnockOnQuadrature::CheckGaussRule() {
    for (len_t n : {1,2,4,8,16}) {
        const Q::GaussLegendreRule rule(n);
        for (len_t power=0; power<2*n; ++power) {
            real_t integral=0;
            for (len_t i=0; i<n; ++i)
                integral+=rule.GetWeight(i)*std::pow(rule.GetNode(i),power);
            if (std::abs(integral-1.0/(power+1))>3e-14) return false;
        }
        real_t mass=0, firstMoment=0;
        for (const auto& point : rule.MappedPoints({-0.7,1.1})) {
            if (!(point.x>-0.7 && point.x<1.1 && point.weight>0)) return false;
            mass+=point.weight; firstMoment+=point.weight*point.x;
        }
        if (std::abs(mass-1.8)>2e-15 || std::abs(firstMoment-0.36)>2e-15) return false;
        if (!rule.MappedPoints({0,0}).empty()) return false;
        if (!rule.MappedPoints({0.2,std::nextafter(0.2,1.0)}).empty()) return false;
    }
    const auto midpoint=Q::GaussLegendreRule(1).MappedPoints({-1,1});
    return midpoint.size()==1 && midpoint[0].x==0 && midpoint[0].weight==2;
}

bool KnockOnQuadrature::CheckValidation() {
    const real_t nan=std::numeric_limits<real_t>::quiet_NaN();
    if (!Throws([] { Q::GaussLegendreRule rule(0); })
        || !Throws([&] { Q::SplitPoints({0,1},{nan}); })
        || !Throws([] { Q::ContactFieldLevels(1.1,0,0.5,0.3); })
        || !Throws([] { Q::ContactFieldLevels(0.5,0,0.5,1); })
        || !Throws([] { Q::PrimaryPitchSplitPoints({1,-1},0.5,0.3,{}); })
        || !Throws([] { Q::OrbitGeometry invalid(nullptr,0); })
        || !Throws([] { Q::OrbitGeometry invalid({0,1,0,M_PI},[](real_t) { return 1; }); }))
        return false;
    Q::OrbitGeometry constant({1,1,0,M_PI},[](real_t) { return 1; });
    if (!Throws([&] { constant.InvertBeta(0.01); })
        || !Throws([&] { constant.InvertLevels({nan}); })
        || !Throws([&] { constant.GetReachableBounds(0.3,0.9,-0.1); })) return false;
    // A rejected reset must not leave a partially built context behind.
    Q::OrbitGeometry geometry({1,1.3/0.7,0,M_PI},[](real_t theta) {
        return 1.3/(1+0.3*std::cos(theta));
    });
    const auto before=geometry.InvertBeta(0.2);
    if (!Throws([&] { geometry.Reset({1,1.01,0,M_PI}); })) return false;
    const auto after=geometry.InvertBeta(0.2);
    if (after.lower!=before.lower || after.upper!=before.upper) return false;
    return true;
}

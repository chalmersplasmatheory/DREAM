#ifndef _DREAMTESTS_KNOCKON_TEST_REFERENCE_HPP
#define _DREAMTESTS_KNOCKON_TEST_REFERENCE_HPP

#include <algorithm>
#include <cmath>
#include <vector>
#include "FVM/Grid/RadialGrid.hpp"

namespace DREAMTESTS::_DREAM::KnockOnTest {
// Independent references: no production contacts, delta integrals, metric
// helpers, or stored cell volumes. Cache geometry, not production answers.
class SurfaceReference {
    struct Node { real_t beta, weight; };
    std::vector<Node> nodes;
    real_t trapped;
public:
    SurfaceReference(const DREAM::FVM::RadialGrid *radial, len_t ir, len_t n) {
        using namespace DREAM::FVM;
        auto *fsa = radial->GetFluxSurfaceAverager();
        const real_t bmin = fsa->GetBmin(ir,FLUXGRIDTYPE_DISTRIBUTION);
        const real_t bmax = fsa->GetBmax(ir,FLUXGRIDTYPE_DISTRIBUTION);
        trapped = bmax==0 ? 0 : std::sqrt(std::max(real_t(0),1-bmin/bmax));
        nodes.reserve(n);
        for (len_t k=0; k<n; ++k) {
            real_t B,J,R,nabla;
            fsa->GeometricQuantitiesAtTheta(ir,2*M_PI*(k+0.5)/n,B,J,R,nabla,FLUXGRIDTYPE_DISTRIBUTION);
            const real_t beta = B==0 ? 0 : std::max(real_t(0),1-bmin/B);
            nodes.push_back({beta,4*M_PI*M_PI*J*(2*M_PI/n)});
        }
    }
    real_t Volume(real_t lo, real_t hi) const {
        if (hi<=lo) return 0;
        real_t result = 0;
        for (const auto& node : nodes) {
            auto eta = [&](real_t u) {
                return std::copysign(std::sqrt(std::max(real_t(0),
                    (u*u-node.beta)/(1-node.beta))),u);
            };
            // Negative passing, positive trapped (two signs), positive passing.
            const real_t lower[] = {-1,0,trapped}, upper[] = {-trapped,trapped,1};
            for (len_t part=0; part<3; ++part) {
                const real_t a = std::max(lo,lower[part]), b = std::min(hi,upper[part]);
                if (b>a) result += node.weight*(eta(b)-eta(a))*(part==1 ? 2 : 1);
            }
        }
        return result;
    }
};

class OrbitReference {
    struct Node { real_t beta, eta, weight; };
    std::vector<Node> nodes;
public:
    OrbitReference(const DREAM::FVM::RadialGrid *radial, len_t ir, real_t u, len_t n) {
        using namespace DREAM::FVM;
        auto *fsa = radial->GetFluxSurfaceAverager();
        real_t minimum, maximum;
        const real_t bmin = fsa->GetBmin(ir,FLUXGRIDTYPE_DISTRIBUTION,&minimum);
        const real_t bmax = fsa->GetBmax(ir,FLUXGRIDTYPE_DISTRIBUTION,&maximum);
        const real_t period = 2*M_PI;
        real_t displacement = std::fmod(maximum-minimum,period);
        if (displacement<0) displacement += period;
        minimum = maximum-displacement;
        real_t left = maximum-period, right = maximum;
        if (bmax>0 && u*u<1-bmin/bmax) {
            // Independent bisection for the two bounce points. Production
            // OrbitGeometry and its inversion brackets are not used here.
            const real_t target = bmin/(1-u*u);
            auto bounce = [&](real_t endpoint) {
                real_t a=minimum, b=endpoint;
                for (len_t k=0; k<60; ++k) {
                    const real_t mid = a+(b-a)/2;
                    if (fsa->BAtTheta(ir,mid,FLUXGRIDTYPE_DISTRIBUTION)<target) a=mid;
                    else b=mid;
                }
                return a+(b-a)/2;
            };
            left = bounce(left); right = bounce(right);
        }
        nodes.reserve(n);
        for (len_t k=0; k<n; ++k) {
            const real_t t = (k+0.5)/n;
            const real_t theta = left+(right-left)*std::pow(std::sin(M_PI*t/2),2);
            const real_t weight = (right-left)*M_PI/2*std::sin(M_PI*t)/n;
            real_t B,J,R,nabla;
            fsa->GeometricQuantitiesAtTheta(ir,theta,B,J,R,nabla,FLUXGRIDTYPE_DISTRIBUTION);
            const real_t ratio = bmin==0 ? 1 : std::max(real_t(1),B/bmin);
            const real_t beta = 1-1/ratio;
            if (beta>=u*u) continue;
            const real_t eta = std::copysign(std::sqrt((u*u-beta)/(1-beta)),u);
            nodes.push_back({beta,eta,weight*4*M_PI*M_PI*J*ratio*std::abs(u/eta)});
        }
    }
    real_t Volume() const {
        real_t sum = 0;
        for (const auto& node : nodes) sum += node.weight;
        return sum;
    }
    real_t Cell(real_t lo, real_t hi, real_t chi) const {
        real_t result = 0;
        for (const auto& node : nodes) {
            auto local = [&](real_t v) {
                return std::copysign(std::sqrt(std::max(real_t(0),
                    (v*v-node.beta)/(1-node.beta))),v);
            };
            const real_t a = chi*node.eta;
            const real_t d = std::sqrt(std::max(real_t(0),(1-node.eta*node.eta)*(1-chi*chi)));
            real_t probability;
            if (d==0) probability = (local(lo)<a && a<=local(hi)) || (a==-1 && lo==-1) ? 1 : 0;
            else {
                auto angle = [](real_t x) { return std::acos(std::clamp(x,real_t(-1),real_t(1))); };
                probability = (angle((local(lo)-a)/d)-angle((local(hi)-a)/d))/M_PI;
            }
            result += node.weight*probability/(hi-lo);
        }
        return result;
    }
};
} // namespace DREAMTESTS::_DREAM::KnockOnTest
#endif

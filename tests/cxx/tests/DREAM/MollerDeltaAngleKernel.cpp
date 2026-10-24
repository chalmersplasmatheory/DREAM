#include "MollerDeltaAngleKernel.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

#include "DREAM/DREAMException.hpp"
#include "DREAM/Equations/Kinetic/MollerDeltaAngleKernel.hpp"
#include "DREAM/Equations/Kinetic/MollerEnergyKernel.hpp"
#include "KnockOnTestGrid.hpp"

using namespace DREAMTESTS::_DREAM;

namespace {
real_t ScatteringCosine(const DREAM::FVM::Grid *target, const DREAM::FVM::Grid *primary,
    len_t i, len_t k, real_t cutoff) {
    const auto *mgK = target->GetMomentumGrid(0), *mgP = primary->GetMomentumGrid(0);
    const real_t p1 = mgP->GetP1(k);
    const real_t lower = std::max(cutoff,mgK->GetP1_f(i));
    const real_t upper = std::min(mgK->GetP1_f(i+1),
        DREAM::KnockOnUtilities::Kinematics::MaximumKnockOnMomentum(p1));
    return DREAM::KnockOnUtilities::Kinematics::EvaluateXiStar(
        std::clamp(mgK->GetP1(i),lower,upper),p1);
}
}

bool MollerDeltaAngleKernel::Run(bool) {
    bool success = true;
    if (CheckMDK_DeltaInterpolationConservation())
        PrintOK("Small native toroidal tables preserve normalization and interpolation endpoints.");
    else {
        success = false;
        PrintError("Angle interpolation/normalization test failed.");
    }
    if (CheckMDK_KinematicEdgeCases())
        PrintOK("Low energies, empty/singleton ranges and selected default-resolution columns passed.");
    else {
        success = false;
        PrintError("Angle kinematic/default-resolution test failed.");
    }
    if (CheckMDK_IsotropicAndDipoleMoments())
        PrintOK("Isotropic shape/zero drift and dipole transfer passed on unlike odd/even grids.");
    else {
        success = false;
        PrintError("Angle isotropic/dipole test failed.");
    }
    if (CheckMDK_GridRebuilt())
        PrintOK("Small angle-table rebuild/lifetime regression passed.");
    else {
        success = false;
        PrintError("Angle-table rebuild regression failed.");
    }
    return success;
}

bool MollerDeltaAngleKernel::CheckMDK_IsotropicAndDipoleMoments() {
    // The moment identities hold at any table resolution. Three slices also
    // exercise interpolation between slices beyond the first interval.
    for (const auto &size : {std::pair<len_t,len_t>{12,15},{15,12}}) {
        const len_t nP = size.first, nK = size.second;
        KnockOnTest::OwnedGrid target(InitializeGridRCylPXi(1,3,nK,1,0,3));
        KnockOnTest::OwnedGrid primary(InitializeGridRCylPXi(1,4,nP,1,0,4));
        const auto *mgK = target->GetMomentumGrid(0), *mgP = primary->GetMomentumGrid(0);
        DREAM::MollerEnergyKernel energy(target.get(),primary.get(),0.1);
        DREAM::MollerDeltaAngleKernel angle(target.get(),primary.get(),0.1,3,8,4);
        std::vector<real_t> isotropic(nP), dipole(nP), output(nK), outputDipole(nK);
        real_t incidentMoment = 0;
        for (len_t l=0; l<nP; ++l) {
            const real_t xi = mgP->GetP2(l), dx = mgP->GetDp2(l);
            isotropic[l] = dx/2;
            dipole[l] = dx*xi;
            incidentMoment += dx*xi*xi;
        }
        const real_t dxK = 2.0/nK, momentTolerance = 2*dxK*dxK;
        len_t active = 0;
        for (len_t i=0; i<mgK->GetNp1(); ++i) for (len_t k=0; k<mgP->GetNp1(); ++k) {
            if (energy.DifferentialCS(i,k)==0) continue;
            ++active;
            std::fill(output.begin(),output.end(),0);
            std::fill(outputDipole.begin(),outputDipole.end(),0);
            angle.AccumulatePitch(0,i,k,isotropic.data(),1,output.data());
            angle.AccumulatePitch(0,i,k,dipole.data(),1,outputDipole.data());
            real_t norm=0, drift=0, number=0, moment=0, even=0, odd=0, shape=0;
            for (len_t j=0; j<nK; ++j) {
                const real_t dx = mgK->GetDp2(j), xi = mgK->GetP2(j);
                norm += dx*output[j];
                drift += dx*xi*output[j];
                number += dx*outputDipole[j];
                moment += dx*xi*outputDipole[j];
                even += dx*std::abs(output[j]-output[nK-1-j]);
                odd += dx*std::abs(outputDipole[j]+outputDipole[nK-1-j]);
                shape += dx*std::abs(output[j]-0.5);
            }
            const real_t chi = ScatteringCosine(target.get(),primary.get(),i,k,0.1);
            if (!std::isfinite(norm) || !std::isfinite(moment) || !std::isfinite(shape)
                || std::abs(norm-1)>5e-13 || std::abs(drift)>5e-13
                || std::abs(number)>5e-13 || even>5e-12 || odd>5e-12
                || shape>1e-3 || std::abs(moment/incidentMoment-chi)>momentTolerance) {
                PrintError("Angular moments: nP=%ld nK=%ld i=%ld k=%ld norm=%g drift=%g shape=%g moment error=%g",
                    (long)nP,(long)nK,(long)i,(long)k,norm,drift,shape,moment/incidentMoment-chi);
                return false;
            }
        }
        if (active==0) return false;
    }
    return true;
}

bool MollerDeltaAngleKernel::CheckMDK_KinematicEdgeCases() {
    struct Case { len_t npK; real_t maxK; len_t npP; real_t maxP; real_t cutoff; bool empty; };
    const Case cases[] = {
        {3,1,4,1,0.1,false},   // no primary reaches the old sqrt(3) threshold
        {1,1,1,1,2,true},      // empty support
        {1,4,1,4,0.2,false},   // singleton scattering angle
        {1,4,3,1,0.1,false},   // target centre above the primary momentum
        {19,1.9,19,1.9,0.1,false} // default table resolution; old angular range clamped valid pairs
    };
    for (len_t cidx=0; cidx<sizeof(cases)/sizeof(cases[0]); ++cidx) {
        const auto& c = cases[cidx];
        KnockOnTest::OwnedGrid target(InitializeGridRCylPXi(1,c.npK,8,1,0,c.maxK));
        KnockOnTest::OwnedGrid primary(InitializeGridRCylPXi(1,c.npP,9,1,0,c.maxP));
        bool rejected = false;
        try {
            DREAM::MollerDeltaAngleKernel invalid(target.get(),primary.get(),c.cutoff,3,8,0);
        } catch (const DREAM::DREAMException&) { rejected = true; }
        if (!rejected) return false;
        DREAM::MollerEnergyKernel energy(target.get(),primary.get(),c.cutoff);
        const bool defaults = cidx==4;
        DREAM::MollerDeltaAngleKernel angle(target.get(),primary.get(),c.cutoff,defaults ? 80 : 3,8,4);
        std::vector<real_t> weights(9,0), pitch(8,0), direct(8,0);
        weights.back() = 1;
        std::vector<std::pair<len_t,len_t>> active;
        for (len_t i=0; i<c.npK; ++i) for (len_t k=0; k<c.npP; ++k)
            if (energy.DifferentialCS(i,k)>0) active.emplace_back(i,k);
        if (c.empty) {
            angle.AccumulatePitch(0,0,0,weights.data(),1,pitch.data());
            if (!active.empty()) return false;
            for (real_t value : pitch) if (!std::isfinite(value) || value!=0) return false;
            continue;
        }
        if (active.empty()) return false;
        // Extremal and interior scattering cosines, plus clipped kinematics,
        // rather than a high-order direct integral for every momentum pair.
        auto cosine = [&](const auto& pair) {
            return ScatteringCosine(target.get(),primary.get(),pair.first,pair.second,c.cutoff);
        };
        std::sort(active.begin(),active.end(),[&](const auto& a,const auto& b) {return cosine(a)<cosine(b);});
        std::vector<std::pair<len_t,len_t>> selected{active.front(),active[active.size()/2],active.back()};
        for (const auto &pair : active)
            if (target->GetMomentumGrid(0)->GetP1(pair.first)
                > DREAM::KnockOnUtilities::Kinematics::MaximumKnockOnMomentum(
                    primary->GetMomentumGrid(0)->GetP1(pair.second))) {
                selected.push_back(pair);
                break;
            }
        std::sort(selected.begin(),selected.end());
        selected.erase(std::unique(selected.begin(),selected.end()),selected.end());
        for (const auto &pair : selected) {
            std::fill(pitch.begin(),pitch.end(),0);
            angle.AccumulatePitch(0,pair.first,pair.second,weights.data(),1,pitch.data());
            real_t norm=0;
            for (len_t j=0; j<pitch.size(); ++j) {
                if (!std::isfinite(pitch[j]) || pitch[j]<0) return false;
                norm += target->GetMomentumGrid(0)->GetDp2(j)*pitch[j];
            }
            if (std::abs(norm-1)>1e-12) return false;
            if (defaults || active.size()==1) {
                DREAM::KnockOnUtilities::SetDeltaMatrixColumnOnGrid(0,cosine(pair),8,
                    target.get(),primary.get(),direct.data(),8,4);
                real_t error=0;
                for (len_t j=0; j<pitch.size(); ++j)
                    error += target->GetMomentumGrid(0)->GetDp2(j)*std::abs(pitch[j]-direct[j]);
                if (!std::isfinite(error) || error>0.02) return false;
            }
        }
    }
    return true;
}

bool MollerDeltaAngleKernel::CheckMDK_DeltaInterpolationConservation() {
    auto target = KnockOnTest::ToroidalGrid(0.3,KnockOnTest::UniformPitchFaces(7),2,3,0,3);
    auto primary = KnockOnTest::ToroidalGrid(0.3,KnockOnTest::UniformPitchFaces(9),2,4,0,4);
    DREAM::MollerEnergyKernel energy(target.get(),primary.get(),0.2);
    DREAM::MollerDeltaAngleKernel angle(target.get(),primary.get(),0.2,3,8,4);
    std::vector<std::pair<len_t,len_t>> active;
    real_t lo=1, hi=0;
    for (len_t i=0; i<3; ++i) for (len_t k=0; k<4; ++k) {
        if (energy.DifferentialCS(i,k)==0) continue;
        active.emplace_back(i,k);
        const real_t chi = ScatteringCosine(target.get(),primary.get(),i,k,0.2);
        lo = std::min(lo,chi); hi = std::max(hi,chi);
    }
    if (!(hi>lo)) return false;
    bool offNode=false, upperSegment=false, sawVoid=false, sawActive=false;
    for (len_t ir=0; ir<2; ++ir) {
        const auto *mg = target->GetMomentumGrid(ir);
        std::vector<real_t> weights(9,0), column(7), direct(7);
        for (const auto &pair : active) {
            const real_t chi = ScatteringCosine(target.get(),primary.get(),pair.first,pair.second,0.2);
            const real_t index = 2*(chi-lo)/(hi-lo);
            offNode |= std::abs(index-std::round(index))>1e-8;
            upperSegment |= index>1;
            for (len_t l=0; l<9; ++l) {
                std::fill(weights.begin(),weights.end(),0);
                std::fill(column.begin(),column.end(),0);
                weights[l] = 1;
                angle.AccumulatePitch(ir,pair.first,pair.second,weights.data(),1,column.data());
                const bool isVoid = primary->GetVpOverP2AtZero(ir)[l]==0;
                sawVoid |= isVoid; sawActive |= !isVoid;
                real_t mass=0;
                for (len_t j=0; j<7; ++j) {
                    if (!std::isfinite(column[j]) || column[j]<0) return false;
                    if (target->GetVpOverP2AtZero(ir)[j]==0 && column[j]!=0) return false;
                    mass += mg->GetDp2(j)*column[j];
                }
                if (std::abs(mass-(isVoid ? 0 : 1))>1e-14) return false;
                if (l==8 && (chi==lo || chi==hi)) {
                    DREAM::KnockOnUtilities::SetDeltaMatrixColumnOnGrid(ir,chi,l,
                        target.get(),primary.get(),direct.data(),8,4);
                    for (len_t j=0; j<7; ++j)
                        if (std::abs(column[j]-direct[j])>1e-8) return false;
                }
            }
        }
    }
    return offNode && upperSegment && sawVoid && sawActive;
}

bool MollerDeltaAngleKernel::CheckMDK_GridRebuilt() {
    KnockOnTest::OwnedGrid target(InitializeGridRCylPXi(1,1,3,1,0,4));
    KnockOnTest::OwnedGrid primary(InitializeGridRCylPXi(1,1,5,1,0,4));
    DREAM::MollerDeltaAngleKernel angle(target.get(),primary.get(),0.2,3,8,4);
    real_t weights[] = {1,0,2,0,3}, before[3] = {}, after[3] = {};
    angle.AccumulatePitch(0,0,0,weights,1,before);
    angle.GridRebuilt();
    angle.AccumulatePitch(0,0,0,weights,1,after);
    real_t mass=0;
    for (len_t j=0; j<3; ++j) {
        if (!std::isfinite(after[j]) || std::abs(after[j]-before[j])>1e-13) return false;
        mass += target->GetMomentumGrid(0)->GetDp2(j)*after[j];
    }
    return std::abs(mass-6)<1e-12;
}

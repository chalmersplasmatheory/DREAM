#include "MollerEnergyKernel.hpp"

#include <algorithm>
#include <cmath>
#include "DREAM/Equations/Kinetic/MollerEnergyKernel.hpp"
#include "DREAM/Equations/KnockOnUtilities.hpp"
#include "KnockOnTestGrid.hpp"

using namespace DREAMTESTS::_DREAM;

bool MollerEnergyKernel::Run(bool) {
    if (CheckMEK_MollerKernelConservation()) {
        PrintOK("Momentum kernel integrates to total cross section at cutoffs inside/on/outside cells.");
        return true;
    }
    PrintError("Momentum kernel conservation failed.");
    return false;
}

bool MollerEnergyKernel::CheckMEK_MollerKernelConservation() {
    // Energy kernels have no radial or pitch dependence. Test momentum support
    // on unlike grids instead of building toroidal bounce data for every cell.
    KnockOnTest::OwnedGrid target(InitializeGridRCylPXi(1,7,1,1,0,3));
    KnockOnTest::OwnedGrid primary(InitializeGridRCylPXi(1,9,1,1,0,4));
    const real_t face = target->GetMomentumGrid(0)->GetP1_f(2);
    for (real_t cutoff : {0.05,face,face-1e-10,face+1e-10,4.0}) {
        DREAM::MollerEnergyKernel kernel(target.get(),primary.get(),cutoff);
        bool anyPositive = false;
        for (len_t k=0; k<primary->GetNp1(0); ++k) {
            real_t integral = 0;
            for (len_t i=0; i<target->GetNp1(0); ++i) {
                const real_t value = kernel.DifferentialCS(i,k);
                if (!std::isfinite(value) || value<0) return false;
                integral += target->GetMomentumGrid(0)->GetDp1(i)*value;
            }
            const real_t expected = kernel.TotalCS(k);
            const real_t reference = DREAM::KnockOnUtilities::EvaluateMollerFluxIntegratedOverKnockonGrid(
                k,target.get(),primary.get(),cutoff);
            // Cross sections are O(1e-20): a fixed tolerance of 1e-14 would
            // allow an identically zero kernel to pass this test.
            const real_t scale = std::max({real_t(1e-100),std::abs(integral),std::abs(reference)});
            if (!std::isfinite(expected) || !std::isfinite(reference) || expected<0
                || std::abs(integral-expected)>1e-12*scale
                || std::abs(expected-reference)>1e-12*scale) {
                PrintError("Energy integral: k=%ld cutoff=%g sum=%g total=%g reference=%g",
                    (long)k,cutoff,integral,expected,reference);
                return false;
            }
            anyPositive |= expected>0;
            if (cutoff==4 && (integral!=0 || expected!=0)) return false;
        }
        if (cutoff<4 && !anyPositive) return false;
    }
    return true;
}

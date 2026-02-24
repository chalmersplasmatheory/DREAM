#ifndef _DREAM_EQUATIONS_MOLLER_KERNEL_HANDLER_HPP
#define _DREAM_EQUATIONS_MOLLER_KERNEL_HANDLER_HPP

#include "DREAM/DREAMException.hpp"
#include "DREAM/Equations/Kinetic/MollerDeltaAngleKernel.hpp"
#include "DREAM/Equations/Kinetic/MollerEnergyKernel.hpp"
#include "DREAM/Equations/KnockOnUtilities.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/TimeKeeper.hpp"
#include "FVM/config.h"

namespace DREAM {
/**
 * Helper to manage the unknown-independent parts of the Boltzmann operator.
 * Its purpose is to coordinate cross sections and angle distributions, since the large-angle
 * collision operator is non-local and couples the distribution function across grids.
 */
class MollerKernelHandler {
   public:
    struct Bundle {
        MollerEnergyKernel *energy = nullptr;
        MollerDeltaAngleKernel *angle = nullptr;
        bool IsValid() const { return energy != nullptr && angle != nullptr; }
        void Deallocate() {
            if (angle != nullptr) {
                delete angle;
                angle = nullptr;
            }
            if (energy != nullptr) {
                delete energy;
                energy = nullptr;
            }
        }
    };

   private:
    real_t pCutoff;
    len_t nXiStars;
    len_t nIntPts;
    KnockOnUtilities::orbit_integration_method integrationMethod;

    const FVM::Grid *hottailGrid = nullptr;
    const FVM::Grid *runawayGrid = nullptr;

    // bundles (may be null if corresponding grids don't exist)
    Bundle hot_hot;
    Bundle hot_re;
    Bundle re_re;

    static void ValidateParams(real_t pCutoff, len_t nXiStars, len_t nIntPts) {
        if (!(pCutoff > 0)) {
            throw DREAMException("MollerKernelHandler: pCutoff must be > 0 (got %.16g).", pCutoff);
        }
        if (nXiStars < 2) {
            throw DREAMException(
                "MollerKernelHandler: nXiStars must be >= 2 (got %ld).", (long)nXiStars
            );
        }
        if (nIntPts < 1) {
            throw DREAMException(
                "MollerKernelHandler: nIntPts must be >= 1 (got %ld).", (long)nIntPts
            );
        }
    }

    static Bundle MakeBundle(
        const FVM::Grid *target, const FVM::Grid *primary, real_t pCutoff, len_t nXiStars,
        len_t nIntPts, KnockOnUtilities::orbit_integration_method integrationMethod
    ) {
        Bundle b;
        b.energy = new MollerEnergyKernel(target, primary, pCutoff);
        b.angle = new MollerDeltaAngleKernel(
            target, primary, pCutoff, nXiStars, nIntPts, integrationMethod
        );
        return b;
    }

    FVM::TimeKeeper *timeKeeper = nullptr;
    len_t timerTot;

   public:
    MollerKernelHandler(
        const FVM::Grid *hottailGrid, const FVM::Grid *runawayGrid, real_t pCutoff, len_t nXiStars,
        len_t nIntPts, KnockOnUtilities::orbit_integration_method integrationMethod
    )
        : pCutoff(pCutoff),
          nXiStars(nXiStars),
          nIntPts(nIntPts),
          integrationMethod(integrationMethod),
          hottailGrid(hottailGrid),
          runawayGrid(runawayGrid) {
        ValidateParams(pCutoff, nXiStars, nIntPts);

        timeKeeper = new FVM::TimeKeeper("MøllerKernels");
        timerTot = timeKeeper->AddTimer("total", "Total Boltzmann initialization");
        timeKeeper->StartTimer(timerTot);
        // Explicitly build only what is possible.
        if (hottailGrid != nullptr) {
            hot_hot =
                MakeBundle(hottailGrid, hottailGrid, pCutoff, nXiStars, nIntPts, integrationMethod);
            if (runawayGrid != nullptr) {
                hot_re = MakeBundle(
                    hottailGrid, runawayGrid, pCutoff, nXiStars, nIntPts, integrationMethod
                );
            }
        }
        if (runawayGrid != nullptr) {
            re_re =
                MakeBundle(runawayGrid, runawayGrid, pCutoff, nXiStars, nIntPts, integrationMethod);
        }
        timeKeeper->StopTimer(timerTot);
    }

    ~MollerKernelHandler() {
        hot_hot.Deallocate();
        hot_re.Deallocate();
        re_re.Deallocate();
        if (timeKeeper != nullptr){
            delete timeKeeper;
            timeKeeper = nullptr;
        }
    }

    // Accessors. Throw if a requested bundle is not available.
    const Bundle &HotHot() const {
        if (!hot_hot.IsValid()) {
            throw DREAMException(
                "MollerKernelHandler: hot-hot kernels requested but hottail grid is not available."
            );
        }
        return hot_hot;
    }
    const Bundle &HotRe() const {
        if (!hot_re.IsValid()) {
            throw DREAMException(
                "MollerKernelHandler: hot-re kernels requested but hottail/runaway grid "
                "combination is not available."
            );
        }
        return hot_re;
    }
    const Bundle &ReRe() const {
        if (!re_re.IsValid()) {
            throw DREAMException(
                "MollerKernelHandler: re-re kernels requested but runaway grid is not available."
            );
        }
        return re_re;
    }

    // Convenience for production-rate terms etc.
    const MollerEnergyKernel *EnergyHotHot() const { return HotHot().energy; }
    const MollerEnergyKernel *EnergyHotRe() const { return HotRe().energy; }
    const MollerEnergyKernel *EnergyReRe() const { return ReRe().energy; }
    const MollerDeltaAngleKernel *AngleHotHot() const { return HotHot().angle; }
    const MollerDeltaAngleKernel *AngleHotRe() const { return HotRe().angle; }
    const MollerDeltaAngleKernel *AngleReRe() const { return ReRe().angle; }

    void GridRebuilt() {
        timeKeeper->StartTimer(timerTot);
        if (hot_hot.energy) hot_hot.energy->GridRebuilt();
        if (hot_hot.angle) hot_hot.angle->GridRebuilt();
        if (hot_re.energy) hot_re.energy->GridRebuilt();
        if (hot_re.angle) hot_re.angle->GridRebuilt();
        if (re_re.energy) re_re.energy->GridRebuilt();
        if (re_re.angle) re_re.angle->GridRebuilt();
        timeKeeper->StopTimer(timerTot);
    }

    real_t GetPCutoff() const { return pCutoff; }
    len_t GetNXiStars() const { return nXiStars; }
    len_t GetNIntPts() const { return nIntPts; }
    KnockOnUtilities::orbit_integration_method GetIntegrationMethod() const {
        return integrationMethod;
    }

    void PrintTimings(bool normalize=false) {
        int_t normalizeto = normalize ? timerTot : -1;
        this->timeKeeper->PrintTimings(true, normalizeto);
    }
    void SaveTimings(SFile *sf, const std::string& path=""){
        this->timeKeeper->SaveTimings(sf, path);
    }
};

}  // namespace DREAM

#endif

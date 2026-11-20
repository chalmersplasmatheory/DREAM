#ifndef _DREAM_EQUATIONS_MOLLER_KERNEL_HANDLER_HPP
#define _DREAM_EQUATIONS_MOLLER_KERNEL_HANDLER_HPP

#include <string>
#include "FVM/config.h"

class SFile;

namespace DREAM::FVM {
    class Grid;
    class TimeKeeper;
}

namespace DREAM {
class MollerEnergyKernel;
class MollerDeltaAngleKernel;

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
        void Deallocate();
    };

   private:
    real_t pCutoff;
    len_t nXiStars;
    len_t poloidalOrder, pitchOrder;

    // bundles (may be null if corresponding grids don't exist)
    Bundle hot_hot;
    Bundle hot_re;
    Bundle re_re;

    static void ValidateParams(
        real_t pCutoff, len_t nXiStars, len_t poloidalOrder, len_t pitchOrder
    );
    static Bundle MakeBundle(
        const FVM::Grid *target, const FVM::Grid *primary, real_t pCutoff, len_t nXiStars,
        len_t poloidalOrder, len_t pitchOrder
    );

    FVM::TimeKeeper *timeKeeper = nullptr;
    len_t timerTot;

   public:
    MollerKernelHandler(
        const FVM::Grid *hottailGrid, const FVM::Grid *runawayGrid, real_t pCutoff, len_t nXiStars,
        len_t poloidalOrder, len_t pitchOrder
    );
    ~MollerKernelHandler();

    // Accessors. Throw if a requested bundle is not available.
    const Bundle &HotHot() const;
    const Bundle &HotRe() const;
    const Bundle &ReRe() const;

    // Convenience for production-rate terms etc.
    const MollerEnergyKernel *EnergyHotHot() const { return HotHot().energy; }
    const MollerEnergyKernel *EnergyHotRe() const { return HotRe().energy; }
    const MollerEnergyKernel *EnergyReRe() const { return ReRe().energy; }
    const MollerDeltaAngleKernel *AngleHotHot() const { return HotHot().angle; }
    const MollerDeltaAngleKernel *AngleHotRe() const { return HotRe().angle; }
    const MollerDeltaAngleKernel *AngleReRe() const { return ReRe().angle; }

    void GridRebuilt();

    real_t GetPCutoff() const { return pCutoff; }
    len_t GetNXiStars() const { return nXiStars; }
    len_t GetPitchOrder() const { return pitchOrder; }
    len_t GetPoloidalOrder() const { return poloidalOrder; }

    void PrintTimings(bool normalize=false);
    void SaveTimings(SFile *sf, const std::string& path="");
};

}  // namespace DREAM

#endif

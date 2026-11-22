#ifndef _DREAM_EQUATIONS_MOLLER_KERNEL_HANDLER_HPP
#define _DREAM_EQUATIONS_MOLLER_KERNEL_HANDLER_HPP

#include <memory>
#include <string>
#include "FVM/config.h"

class SFile;

namespace DREAM::FVM {
    class Grid;
    class TimeKeeper;
    class UnknownQuantityHandler;
}

namespace DREAM {
class MollerEnergyKernel;
class MollerDeltaAngleKernel;
class RosenbluthPutvinskiBoltzmann;

/**
 * Helper to manage the unknown-independent parts of the Boltzmann operator.
 * Its purpose is to coordinate cross sections and angle distributions, since the large-angle
 * collision operator is non-local and couples the distribution function across grids.
 */
class MollerKernelHandler {
   public:
    struct Bundle {
        std::unique_ptr<MollerEnergyKernel> energy;
        std::unique_ptr<MollerDeltaAngleKernel> angle;
        bool IsValid() const { return energy != nullptr && angle != nullptr; }
    };

   private:
    real_t pCutoff;
    len_t nXiStars;
    len_t poloidalOrder, pitchOrder;

    // bundles (may be null if corresponding grids don't exist)
    Bundle hot_hot;
    Bundle hot_re;
    Bundle re_re;
    std::unique_ptr<RosenbluthPutvinskiBoltzmann> rp;

    static void ValidateParams(
        real_t pCutoff, len_t nXiStars, len_t poloidalOrder, len_t pitchOrder
    );
    static Bundle MakeBundle(
        const FVM::Grid *target, const FVM::Grid *primary, real_t pCutoff, len_t nXiStars,
        len_t poloidalOrder, len_t pitchOrder
    );

    std::unique_ptr<FVM::TimeKeeper> timeKeeper;
    len_t timerTot;

   public:
    MollerKernelHandler(
        const FVM::Grid *hottailGrid, const FVM::Grid *runawayGrid,
        FVM::UnknownQuantityHandler *unknowns, real_t pCutoff, len_t nXiStars,
        len_t poloidalOrder, len_t pitchOrder
    );
    ~MollerKernelHandler();

    RosenbluthPutvinskiBoltzmann *GetRP() const { return rp.get(); }

    bool HasHotHot() const { return hot_hot.IsValid(); }
    bool HasHotRe() const { return hot_re.IsValid(); }
    bool HasReRe() const { return re_re.IsValid(); }

    // Accessors. Throw if a requested bundle is not available.
    const Bundle &HotHot() const;
    const Bundle &HotRe() const;
    const Bundle &ReRe() const;

    // Convenience for production-rate terms etc.
    const MollerEnergyKernel *EnergyHotHot() const { return HotHot().energy.get(); }
    const MollerEnergyKernel *EnergyHotRe() const { return HotRe().energy.get(); }
    const MollerEnergyKernel *EnergyReRe() const { return ReRe().energy.get(); }
    const MollerDeltaAngleKernel *AngleHotHot() const { return HotHot().angle.get(); }
    const MollerDeltaAngleKernel *AngleHotRe() const { return HotRe().angle.get(); }
    const MollerDeltaAngleKernel *AngleReRe() const { return ReRe().angle.get(); }

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

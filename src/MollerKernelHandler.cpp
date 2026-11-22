/**
 * Owns the finite-energy kernels and the shared infinite-energy spectrum
 * used by the kinetic gain and fluid production terms.
 */
#include "DREAM/MollerKernelHandler.hpp"
#include "DREAM/DREAMException.hpp"
#include "DREAM/Equations/Kinetic/MollerDeltaAngleKernel.hpp"
#include "DREAM/Equations/Kinetic/MollerEnergyKernel.hpp"
#include "DREAM/Equations/RosenbluthPutvinskiBoltzmann.hpp"
#include "FVM/TimeKeeper.hpp"

using namespace DREAM;

void MollerKernelHandler::ValidateParams(
    real_t pCutoff, len_t nXiStars, len_t poloidalOrder, len_t pitchOrder
) {
    if (!(pCutoff > 0)) {
        throw DREAMException("MollerKernelHandler: pCutoff must be > 0 (got %.16g).", pCutoff);
    }
    if (nXiStars < 2) {
        throw DREAMException(
            "MollerKernelHandler: nXiStars must be >= 2 (got %ld).", (long)nXiStars
        );
    }
    if (pitchOrder < 1)
        throw DREAMException("MollerKernelHandler: pitchOrder must be >= 1.");
    if (poloidalOrder < 1) {
        throw DREAMException(
            "MollerKernelHandler: poloidalOrder must be >= 1 (got %ld).", (long)poloidalOrder
        );
    }
}

MollerKernelHandler::Bundle MollerKernelHandler::MakeBundle(
    const FVM::Grid *target, const FVM::Grid *primary, real_t pCutoff, len_t nXiStars,
    len_t poloidalOrder, len_t pitchOrder, OptionConstants::boltzmann_model model
) {
    return {
        std::make_unique<MollerEnergyKernel>(target, primary, pCutoff),
        std::make_unique<MollerDeltaAngleKernel>(
            target, primary, pCutoff, nXiStars, poloidalOrder, pitchOrder, model
        )
    };
}

MollerKernelHandler::MollerKernelHandler(
    const FVM::Grid *hottailGrid, const FVM::Grid *runawayGrid,
    FVM::UnknownQuantityHandler *unknowns, real_t pCutoff, len_t nXiStars,
    len_t poloidalOrder, len_t pitchOrder, OptionConstants::boltzmann_model model
)
    : pCutoff(pCutoff),
      nXiStars(nXiStars),
      poloidalOrder(poloidalOrder),
      pitchOrder(pitchOrder) {
    ValidateParams(pCutoff, nXiStars, poloidalOrder, pitchOrder);

    timeKeeper = std::make_unique<FVM::TimeKeeper>("MøllerKernels");
    timerTot = timeKeeper->AddTimer("total", "Total Boltzmann initialization");
    timeKeeper->StartTimer(timerTot);

    rp = std::make_unique<RosenbluthPutvinskiBoltzmann>(
        hottailGrid, runawayGrid, unknowns, pCutoff, poloidalOrder
    );

    // Build kernels for each target and primary grid pair that exists.
    if (hottailGrid != nullptr) {
        hot_hot = MakeBundle(
            hottailGrid, hottailGrid, pCutoff, nXiStars, poloidalOrder, pitchOrder, model
        );
        if (runawayGrid != nullptr) {
            hot_re = MakeBundle(
                hottailGrid, runawayGrid, pCutoff, nXiStars, poloidalOrder, pitchOrder, model
            );
        }
    }
    if (runawayGrid != nullptr) {
        re_re = MakeBundle(
            runawayGrid, runawayGrid, pCutoff, nXiStars, poloidalOrder, pitchOrder, model
        );
    }
    timeKeeper->StopTimer(timerTot);
}

MollerKernelHandler::~MollerKernelHandler() = default;

const MollerKernelHandler::Bundle &MollerKernelHandler::HotHot() const {
    if (!hot_hot.IsValid()) {
        throw DREAMException(
            "MollerKernelHandler: hot-hot kernels requested but hottail grid is not available."
        );
    }
    return hot_hot;
}

const MollerKernelHandler::Bundle &MollerKernelHandler::HotRe() const {
    if (!hot_re.IsValid()) {
        throw DREAMException(
            "MollerKernelHandler: hot-re kernels requested but hottail/runaway grid "
            "combination is not available."
        );
    }
    return hot_re;
}

const MollerKernelHandler::Bundle &MollerKernelHandler::ReRe() const {
    if (!re_re.IsValid()) {
        throw DREAMException(
            "MollerKernelHandler: re-re kernels requested but runaway grid is not available."
        );
    }
    return re_re;
}

void MollerKernelHandler::GridRebuilt() {
    timeKeeper->StartTimer(timerTot);
    rp->GridRebuilt();
    if (hot_hot.energy) hot_hot.energy->GridRebuilt();
    if (hot_hot.angle) hot_hot.angle->GridRebuilt();
    if (hot_re.energy) hot_re.energy->GridRebuilt();
    if (hot_re.angle) hot_re.angle->GridRebuilt();
    if (re_re.energy) re_re.energy->GridRebuilt();
    if (re_re.angle) re_re.angle->GridRebuilt();
    timeKeeper->StopTimer(timerTot);
}

void MollerKernelHandler::PrintTimings(bool normalize) {
    int_t normalizeto = normalize ? timerTot : -1;
    this->timeKeeper->PrintTimings(true, normalizeto);
}

void MollerKernelHandler::SaveTimings(SFile *sf, const std::string& path) {
    this->timeKeeper->SaveTimings(sf, path);
}

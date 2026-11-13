/**
 * Owns and updates the Møller kernels shared by the kinetic source and
 * fluid production terms.
 */
#include "DREAM/MollerKernelHandler.hpp"
#include "DREAM/DREAMException.hpp"
#include "DREAM/Equations/Kinetic/MollerDeltaAngleKernel.hpp"
#include "DREAM/Equations/Kinetic/MollerEnergyKernel.hpp"
#include "FVM/TimeKeeper.hpp"

using namespace DREAM;

void MollerKernelHandler::Bundle::Deallocate() {
    delete angle;
    delete energy;
    angle = nullptr;
    energy = nullptr;
}

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
    len_t poloidalOrder, len_t pitchOrder
) {
    Bundle b;
    b.energy = new MollerEnergyKernel(target, primary, pCutoff);
    b.angle = new MollerDeltaAngleKernel(
        target, primary, pCutoff, nXiStars, poloidalOrder, pitchOrder
    );
    return b;
}

MollerKernelHandler::MollerKernelHandler(
    const FVM::Grid *hottailGrid, const FVM::Grid *runawayGrid, real_t pCutoff, len_t nXiStars,
    len_t poloidalOrder, len_t pitchOrder
)
    : pCutoff(pCutoff),
      nXiStars(nXiStars),
      poloidalOrder(poloidalOrder),
      pitchOrder(pitchOrder) {
    ValidateParams(pCutoff, nXiStars, poloidalOrder, pitchOrder);

    timeKeeper = new FVM::TimeKeeper("MøllerKernels");
    timerTot = timeKeeper->AddTimer("total", "Total Boltzmann initialization");
    timeKeeper->StartTimer(timerTot);

    // Build kernels for each target and primary grid pair that exists.
    if (hottailGrid != nullptr) {
        hot_hot =
            MakeBundle(hottailGrid, hottailGrid, pCutoff, nXiStars, poloidalOrder, pitchOrder);
        if (runawayGrid != nullptr) {
            hot_re = MakeBundle(
                hottailGrid, runawayGrid, pCutoff, nXiStars, poloidalOrder, pitchOrder
            );
        }
    }
    if (runawayGrid != nullptr) {
        re_re =
            MakeBundle(runawayGrid, runawayGrid, pCutoff, nXiStars, poloidalOrder, pitchOrder);
    }
    timeKeeper->StopTimer(timerTot);
}

MollerKernelHandler::~MollerKernelHandler() {
    hot_hot.Deallocate();
    hot_re.Deallocate();
    re_re.Deallocate();
    delete timeKeeper;
}

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

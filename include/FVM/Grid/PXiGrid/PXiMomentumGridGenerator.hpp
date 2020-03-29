#ifndef _DREAM_FVM_P_XI_MOMENTUM_GRID_GENERATOR_HPP
#define _DREAM_FVM_P_XI_MOMENTUM_GRID_GENERATOR_HPP

#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/PXiGrid/PGridGenerator.hpp"
#include "FVM/Grid/PXiGrid/XiGridGenerator.hpp"

namespace DREAM::FVM::PXiGrid {
    class MomentumGridGenerator : public DREAM::FVM::MomentumGridGenerator {
    private:
        PGridGenerator *pGenerator;
        XiGridGenerator *xiGenerator;

    public:
        MomentumGridGenerator(PGridGenerator *pg, XiGridGenerator *xg)
            : pGenerator(pg), xiGenerator(xg) {}

        virtual ~MomentumGridGenerator();

        virtual bool NeedsRebuild(const real_t, const bool) override;
        virtual bool Rebuild(const real_t, const len_t, DREAM::FVM::MomentumGrid*, const DREAM::FVM::RadialGrid*) override;
    };
}

#endif/*_DREAM_FVM_P_XI_MOMENTUM_GRID_GENERATOR_HPP*/

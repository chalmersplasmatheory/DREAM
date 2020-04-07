#ifndef _DREAM_FVM_P_XI_GRID_P_GRID_GENERATOR_HPP
#define _DREAM_FVM_P_XI_GRID_P_GRID_GENERATOR_HPP

#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/RadialGrid.hpp"

namespace DREAM::FVM::PXiGrid {
    class PGridGenerator {
    public:
        virtual ~PGridGenerator() {}
        virtual len_t GetNp() const = 0;

        virtual bool NeedsRebuild(const real_t, const bool) = 0;
        virtual bool Rebuild(
            const real_t, const len_t, DREAM::FVM::MomentumGrid*, const DREAM::FVM::RadialGrid *rg
        ) = 0;
    };
}

#endif/*_DREAM_FVM_P_XI_GRID_P_GRID_GENERATOR_HPP*/

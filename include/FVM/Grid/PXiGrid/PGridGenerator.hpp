#ifndef _TQS_FVM_P_XI_GRID_P_GRID_GENERATOR_HPP
#define _TQS_FVM_P_XI_GRID_P_GRID_GENERATOR_HPP

#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/RadialGrid.hpp"

namespace TQS::FVM::PXiGrid {
    class PGridGenerator {
    public:
        virtual bool NeedsRebuild(const real_t, const bool) = 0;
        virtual bool Rebuild(
            const real_t, const len_t, TQS::FVM::MomentumGrid*, const TQS::FVM::RadialGrid *rg
        ) = 0;
    };
}

#endif/*_TQS_FVM_P_XI_GRID_P_GRID_GENERATOR_HPP*/

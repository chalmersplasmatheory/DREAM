#ifndef _TQS_FVM_P_XI_GRID_P_GRID_GENERATOR_HPP
#define _TQS_FVM_P_XI_GRID_P_GRID_GENERATOR_HPP

namespace TQS::FVM::PXiGrid {
    class PGridGenerator {
    public:
        virtual bool NeedsRebuild(const real_t, const bool) = 0;
        virtual bool Rebuild(
            const real_t, const len_t, const TQS::FVM::MomentumGrid*, TQS::FVM::RadialGrid *rg
        );
    };
}

#endif/*_TQS_FVM_P_XI_GRID_P_GRID_GENERATOR_HPP*/

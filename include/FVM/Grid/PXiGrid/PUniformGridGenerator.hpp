#ifndef _TQS_FVM_P_XI_GRID_P_UNIFORM_GRID_GENERATOR_HPP
#define _TQS_FVM_P_XI_GRID_P_UNIFORM_GRID_GENERATOR_HPP

#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/PGridGenerator.hpp"

namespace TQS::FVM::PXiGrid {
    class PUniformGridGenerator : public PGridGenerator {
    private:
        len_t np;
        real_t pMin, pMax;

        bool initialized = false;
    public:
        PUniformGridGenerator(const len_t, const real_t, const real_t);

        virtual bool NeedsRebuild(const real_t, const bool) { return (!initialized); }
        virtual bool Rebuild(const real_t, const len_t, MomentumGrid*, const RadialGrid*);
    };
}

#endif/*_TQS_FVM_P_XI_GRID_P_UNIFORM_GRID_GENERATOR_HPP*/

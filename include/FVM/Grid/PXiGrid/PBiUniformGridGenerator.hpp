#ifndef _DREAM_FVM_P_XI_GRID_P_BI_UNIFORM_GRID_GENERATOR_HPP
#define _DREAM_FVM_P_XI_GRID_P_BI_UNIFORM_GRID_GENERATOR_HPP

#include "FVM/FVMException.hpp"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/PXiGrid/PGridGenerator.hpp"

namespace DREAM::FVM::PXiGrid {
    class PBiUniformGridGenerator : public PGridGenerator {
    private:
        len_t np, npSep;
        real_t pMin, pSep, pMax;

        bool initialized = false;
    public:
        PBiUniformGridGenerator(const len_t, const len_t, const real_t, const real_t, const real_t);
        PBiUniformGridGenerator(const len_t, const real_t, const real_t, const real_t, const real_t);

        virtual len_t GetNp() const { return this->np; }

        virtual bool NeedsRebuild(const real_t, const bool) { return (!initialized); }
        virtual bool Rebuild(const real_t, const len_t, MomentumGrid*, const RadialGrid*);
    };
}

#endif/*_DREAM_FVM_P_XI_GRID_P_BI_UNIFORM_GRID_GENERATOR_HPP*/

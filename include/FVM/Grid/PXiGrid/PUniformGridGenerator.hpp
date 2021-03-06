#ifndef _DREAM_FVM_P_XI_GRID_P_UNIFORM_GRID_GENERATOR_HPP
#define _DREAM_FVM_P_XI_GRID_P_UNIFORM_GRID_GENERATOR_HPP

#include "FVM/FVMException.hpp"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/PXiGrid/PGridGenerator.hpp"

namespace DREAM::FVM::PXiGrid {
    class PUniformGridGenerator : public PGridGenerator {
    private:
        len_t np;
        real_t pMin, pMax;

        bool initialized = false;
    public:
        PUniformGridGenerator(const len_t, const real_t, const real_t);

        virtual len_t GetNp() const { return this->np; }

        virtual bool NeedsRebuild(const real_t, const bool) { return (!initialized); }
        virtual bool Rebuild(const real_t, const len_t, MomentumGrid*, const RadialGrid*);
    };
}

#endif/*_DREAM_FVM_P_XI_GRID_P_UNIFORM_GRID_GENERATOR_HPP*/

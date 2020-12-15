#ifndef _DREAM_FVM_P_XI_GRID_P_CUSTOM_GRID_GENERATOR_HPP
#define _DREAM_FVM_P_XI_GRID_P_CUSTOM_GRID_GENERATOR_HPP

#include "FVM/FVMException.hpp"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/PXiGrid/PGridGenerator.hpp"

namespace DREAM::FVM::PXiGrid {
    class PCustomGridGenerator : public PGridGenerator {
    private:
        real_t *pf_provided;
        len_t np;

        bool initialized = false;
    public:
        PCustomGridGenerator(const real_t *p_f, const len_t np);

        virtual len_t GetNp() const { return this->np; }

        virtual bool NeedsRebuild(const real_t, const bool) { return (!initialized); }
        virtual bool Rebuild(const real_t, const len_t, MomentumGrid*, const RadialGrid*);
    };
}

#endif/*_DREAM_FVM_P_XI_GRID_P_CUSTOM_GRID_GENERATOR_HPP*/

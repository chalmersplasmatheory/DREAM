#ifndef _DREAM_FVM_P_XI_GRID_XI_CUSTOM_GRID_GENERATOR_HPP
#define _DREAM_FVM_P_XI_GRID_XI_CUSTOM_GRID_GENERATOR_HPP
#include "FVM/FVMException.hpp"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/PXiGrid/XiGridGenerator.hpp"

namespace DREAM::FVM::PXiGrid {
    class XiCustomGridGenerator : public XiGridGenerator {
    private:
        real_t *xif_provided;
        len_t nxi;

        bool initialized = false;
    public:
        XiCustomGridGenerator(const real_t *xi_f, const len_t nxi);

        virtual len_t GetNxi() const { return this->nxi; }

        virtual bool NeedsRebuild(const real_t, const bool) { return (!initialized); }
        virtual bool Rebuild(const real_t, const len_t, MomentumGrid*, const RadialGrid*);
    };
}

#endif/*_DREAM_FVM_P_XI_GRID_XI_CUSTOM_GRID_GENERATOR_HPP*/

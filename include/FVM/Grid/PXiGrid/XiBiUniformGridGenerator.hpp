#ifndef _DREAM_FVM_P_XI_GRID_XI_BIUNIFORM_GRID_GENERATOR_HPP
#define _DREAM_FVM_P_XI_GRID_XI_BIUNIFORM_GRID_GENERATOR_HPP

#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/PXiGrid/XiGridGenerator.hpp"

namespace DREAM::FVM::PXiGrid {
    class XiBiUniformGridGenerator : public XiGridGenerator {
    private:
        len_t nxi, nxiSep;
        real_t xiMin=-1, xiMax=1, xiSep;

        bool initialized = false;
    public:
        XiBiUniformGridGenerator(const len_t, const len_t, const real_t);
        XiBiUniformGridGenerator(const len_t, const real_t, const real_t);

        len_t GetNxi() const { return this->nxi; }

        virtual bool NeedsRebuild(const real_t, const bool) { return (!initialized); }
        virtual bool Rebuild(const real_t, const len_t, MomentumGrid*, const RadialGrid*);
    };
}

#endif/*_DREAM_FVM_P_XI_GRID_XI_BIUNIFORM_GRID_GENERATOR_HPP*/

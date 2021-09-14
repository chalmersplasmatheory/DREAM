#ifndef _DREAM_FVM_P_XI_GRID_XI_TRAPPED_PASSING_BOUNDARY_LAYER_GRID_GENERATOR_HPP
#define _DREAM_FVM_P_XI_GRID_XI_TRAPPED_PASSING_BOUNDARY_LAYER_GRID_GENERATOR_HPP

#include "FVM/FVMException.hpp"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/PXiGrid/XiGridGenerator.hpp"

namespace DREAM::FVM::PXiGrid {
    class XiTrappedPassingBoundaryLayerGridGenerator : public XiGridGenerator {
    private:
        real_t dxiMax;
        len_t nxiPass, nxiTrap;
        real_t boundaryLayerWidth;
        
        len_t nxi = 0;

        bool initialized = false;
    public:
        XiTrappedPassingBoundaryLayerGridGenerator(
            const real_t dxiMax=2, const len_t nxiPass=1, const len_t nxiTrap=1,
            const real_t boundaryLayerWidth=1e-3
        );

        virtual len_t GetNxi() const override { return this->nxi; }

        virtual bool NeedsRebuild(const real_t, const bool) override { return (!initialized); }
        virtual bool Rebuild(const real_t, const len_t, MomentumGrid*, const RadialGrid*) override;
    };
}

#endif/*_DREAM_FVM_P_XI_GRID_XI_TRAPPED_PASSING_BOUNDARY_LAYER_GRID_GENERATOR_HPP*/

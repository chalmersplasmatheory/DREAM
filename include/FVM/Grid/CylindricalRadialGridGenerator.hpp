#ifndef _DREAM_FVM_GRID_CYLINDRICAL_RADIAL_GRID_GENERATOR_HPP
#define _DREAM_FVM_GRID_CYLINDRICAL_RADIAL_GRID_GENERATOR_HPP

#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/RadialGridGenerator.hpp"
#include <functional>

namespace DREAM::FVM {
    class CylindricalRadialGridGenerator : public RadialGridGenerator {
    private:
        real_t xMin=0, xMax=1;
        real_t B0=0;
        
        // Set to true when the grid is constructed for the first time
        bool isBuilt = false;

    public:
        CylindricalRadialGridGenerator(const len_t nx, const real_t B0, const real_t x0=0, const real_t xa=1);

        virtual bool NeedsRebuild(const real_t) const override { return (!isBuilt); }
        virtual bool Rebuild(const real_t, RadialGrid*) override;
        virtual void RebuildJacobians(RadialGrid*, MomentumGrid**,MagneticQuantityHandler*) override;
        virtual void RebuildFSAvgQuantities(RadialGrid*, MomentumGrid**) override;
        virtual real_t BounceAverageQuantity(RadialGrid*, const MomentumGrid*, len_t, len_t, len_t, len_t, std::function<real_t(real_t,real_t)>) override;
        virtual real_t FluxSurfaceAverageQuantity(RadialGrid*,len_t,bool, std::function<real_t(real_t)>) override;
    };
}

#endif/*_DREAM_FVM_GRID_CYLINDRICAL_RADIAL_GRID_GENERATOR_HPP*/

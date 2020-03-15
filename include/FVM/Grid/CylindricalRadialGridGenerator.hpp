#ifndef _TQS_FVM_GRID_CYLINDRICAL_RADIAL_GRID_GENERATOR_HPP
#define _TQS_FVM_GRID_CYLINDRICAL_RADIAL_GRID_GENERATOR_HPP

#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/RadialGridGenerator.hpp"

namespace TQS::FVM {
    class CylindricalRadialGridGenerator : public RadialGridGenerator {
    private:
        len_t nx=0;
        real_t xMin=0, xMax=1;
        real_t B0=0;

        // Set to true when the grid is constructed for the first time
        bool isBuilt = false;

    public:
        CylindricalRadialGridGenerator(const len_t, const real_t, const real_t x0=0, const real_t xa=1);

        virtual bool NeedsRebuild(const real_t) override { return (!isBuilt); }
        virtual bool Rebuild(const real_t, RadialGrid*) override;
    };
}

#endif/*_TQS_FVM_GRID_CYLINDRICAL_RADIAL_GRID_GENERATOR_HPP*/

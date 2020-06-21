#ifndef _DREAM_EMPTY_RADIAL_GRID_HPP
#define _DREAM_EMPTY_RADIAL_GRID_HPP

#include "FVM/Grid/RadialGridGenerator.hpp"

namespace DREAM::FVM {
    class EmptyRadialGridGenerator : public RadialGridGenerator {
        public:
            EmptyRadialGridGenerator() : RadialGridGenerator(1) {
                ntheta_ref = 2; ntheta_interp = 1;
            }

            virtual bool NeedsRebuild(const real_t) const override { return false; }
            virtual bool Rebuild(const real_t, RadialGrid*) override;
            virtual void CreateMagneticFieldData(const real_t *x, const real_t *x_f) override;

    };

    class EmptyRadialGrid : public RadialGrid {
        public:
            EmptyRadialGrid(): RadialGrid(new EmptyRadialGridGenerator()) {}

            virtual ~EmptyRadialGrid() {}

    };
}

#endif/*_DREAM_EMPTY_RADIAL_GRID_HPP*/

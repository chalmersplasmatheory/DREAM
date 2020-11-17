#ifndef _DREAM_EMPTY_MOMENTUM_GRID_HPP
#define _DREAM_EMPTY_MOMENTUM_GRID_HPP

#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/MomentumGridGenerator.hpp"
#include "FVM/Grid/RadialGrid.hpp"

namespace DREAM::FVM {
    class EmptyMomentumGridGenerator : public DREAM::FVM::MomentumGridGenerator {
    private:
        bool isBuilt = false;
    public:
        EmptyMomentumGridGenerator() {}

        virtual bool NeedsRebuild(const real_t, const bool) override { return !isBuilt; }
        virtual bool Rebuild(const real_t, const len_t, DREAM::FVM::MomentumGrid*, const DREAM::FVM::RadialGrid*) override;
    };

    class EmptyMomentumGrid : public DREAM::FVM::MomentumGrid {
    public:
        EmptyMomentumGrid(const RadialGrid *rGrid, const real_t t0=0)
            : DREAM::FVM::MomentumGrid(new EmptyMomentumGridGenerator(), 0, rGrid, t0) {}

        virtual ~EmptyMomentumGrid() {}


        virtual void EvaluateMetricOverP2(
           const len_t i, const len_t j ,
            fluxGridType fluxGridType, 
            const len_t ntheta, const real_t* theta,
            const real_t* BOverBmin, real_t *&sqrtg
        ) const override;
    };
}

#endif/*_DREAM_EMPTY_MOMENTUM_GRID_HPP*/

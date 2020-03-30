#ifndef _DREAM_EMPTY_MOMENTUM_GRID_HPP
#define _DREAM_EMPTY_MOMENTUM_GRID_HPP

#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/MomentumGridGenerator.hpp"
#include "FVM/Grid/RadialGrid.hpp"

namespace DREAM::FVM {
    class EmptyMomentumGridGenerator : public DREAM::FVM::MomentumGridGenerator {
    private:
    public:
        EmptyMomentumGridGenerator() {}

        virtual bool NeedsRebuild(const real_t, const bool) override { return false; }
        virtual bool Rebuild(const real_t, const len_t, DREAM::FVM::MomentumGrid*, const DREAM::FVM::RadialGrid*) override;
    };

    class EmptyMomentumGrid : public DREAM::FVM::MomentumGrid {
    public:
        EmptyMomentumGrid(const RadialGrid *rGrid, const real_t t0=0)
            : DREAM::FVM::MomentumGrid(new EmptyMomentumGridGenerator(), 0, rGrid, t0) {}

        virtual ~EmptyMomentumGrid() {}


        virtual void EvaluateMetric(
            const real_t p1, const real_t p2,
            const len_t ir, const RadialGrid *rGrid,
            const len_t ntheta, const real_t *theta,
            bool rFluxGrid, real_t *sqrtg
        ) const override;
    };
}

#endif/*_DREAM_EMPTY_MOMENTUM_GRID_HPP*/

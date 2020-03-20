#ifndef _TQS_FVM_P_XI_MOMENTUM_GRID_GENERATOR_HPP
#define _TQS_FVM_P_XI_MOMENTUM_GRID_GENERATOR_HPP

#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/PXiGrid/PGridGenerator.hpp"
#include "FVM/Grid/PXiGrid/XiGridGenerator.hpp"

namespace TQS::FVM::PXiGrid {
    class MomentumGridGenerator : public TQS::FVM::MomentumGridGenerator {
    private:
        PGridGenerator *pGenerator;
        XiGridGenerator *xiGenerator;

    protected:
        void GenerateLameCoeffs(
            const len_t, const len_t, const real_t*, const real_t*,
            real_t**, real_t**, real_t**
        );

    public:
        MomentumGridGenerator(PGridGenerator *pg, XiGridGenerator *xg)
            : pGenerator(pg), xiGenerator(xg) {}

        virtual ~MomentumGridGenerator();

        virtual bool NeedsRebuild(const real_t, const bool) override;
        virtual bool Rebuild(const real_t, const len_t, TQS::FVM::MomentumGrid*, const TQS::FVM::RadialGrid*) override;
    };
}

#endif/*_TQS_FVM_P_XI_MOMENTUM_GRID_GENERATOR_HPP*/

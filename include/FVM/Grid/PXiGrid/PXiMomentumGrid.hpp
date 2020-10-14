#ifndef _DREAM_FVM_PXI_MOMENTUM_GRID_HPP
#define _DREAM_FVM_PXI_MOMENTUM_GRID_HPP

#include "FVM/Grid/PXiGrid/PXiMomentumGridGenerator.hpp"

namespace DREAM::FVM::PXiGrid {
    class PXiMomentumGrid : public DREAM::FVM::MomentumGrid {
    public:
        PXiMomentumGrid(MomentumGridGenerator *g, const len_t ir, const RadialGrid *rGrid, const real_t t0=0)
            : DREAM::FVM::MomentumGrid(g, ir, rGrid, t0) {}

        virtual ~PXiMomentumGrid() {}

        // Translation from P1/P2  -->  P/XI
        const real_t *GetP() const { return this->GetP1(); }
        const real_t *GetXi() const { return this->GetP2(); }
        const real_t *GetP_f() const { return this->GetP1_f(); }
        const real_t *GetXi_f() const { return this->GetP2_f(); }
        const real_t *GetDp() const { return this->GetDp1(); }
        const real_t *GetDxi() const { return this->GetDp2(); }
        const real_t *GetDp_f() const { return this->GetDp1_f(); }
        const real_t *GetDxi_f() const { return this->GetDp2_f(); }

        virtual void EvaluateMetric(
            const len_t i, const len_t j ,
            fluxGridType fluxGridType, 
            const len_t ntheta, const real_t* theta,
            const real_t* BOverBmin, real_t *&sqrtg
        ) const override;

    };
}

#endif/*_DREAM_FVM_PXI_MOMENTUM_GRID_HPP*/

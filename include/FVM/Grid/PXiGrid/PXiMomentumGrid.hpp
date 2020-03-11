#ifndef _TQS_FVM_PXI_MOMENTUM_GRID_HPP
#define _TQS_FVM_PXI_MOMENTUM_GRID_HPP

#include "FVM/Grid/PXiMomentumGridGenerator.hpp"

namespace TQS::FVM::PXiGrid {
    class PXiMomentumGrid : public MomentumGrid {
    public:
        PXiMomentumGrid(PXiMomentumGridGenerator *g, const real_t t0=0)
            : MomentumGrid(g, t0) {}

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
    };
}

#endif/*_TQS_FVM_PXI_MOMENTUM_GRID_HPP*/

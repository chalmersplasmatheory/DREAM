#ifndef _DREAM_FVM_GRID_HPP
#define _DREAM_FVM_GRID_HPP

namespace DREAM::FVM { class Grid; }

#include "FVM/config.h"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
//#include "FVM/Grid/fluxGridType.enum.hpp"

namespace DREAM::FVM {
    class Grid {
    private:
    protected:
        RadialGrid *rgrid;
		MomentumGrid **momentumGrids;

    public:
        Grid(RadialGrid*, MomentumGrid*, const real_t t0=0);
        ~Grid();

        // Returns pointer to the momentum grid with the specified index
        MomentumGrid *GetMomentumGrid(const len_t i) const { return this->momentumGrids[i]; }
        RadialGrid *GetRadialGrid() const { return this->rgrid; }

        const len_t GetNCells() const;
        const len_t GetNr() const { return this->rgrid->GetNr(); }

        real_t *const* GetVp() const { return this->rgrid->GetVp(); }
        const real_t *GetVp(const len_t ir) const { return this->rgrid->GetVp(ir); }
        real_t *const* GetVp_fr() const { return this->rgrid->GetVp_fr(); }
        const real_t *GetVp_fr(const len_t ir) const { return this->rgrid->GetVp_fr(ir); }
        real_t *const* GetVp_f1() const { return this->rgrid->GetVp_f1(); }
        const real_t *GetVp_f1(const len_t ir) const { return this->rgrid->GetVp_f1(ir); }
        real_t *const* GetVp_f2() const { return this->rgrid->GetVp_f2(); }
        const real_t *GetVp_f2(const len_t ir) const { return this->rgrid->GetVp_f2(ir); }

        const real_t *GetVpVol() const {return this->rgrid->GetVpVol(); }
        const real_t  GetVpVol(const len_t ir) const {return this->rgrid->GetVpVol(ir); }
        const real_t *GetVpVol_f() const {return this->rgrid->GetVpVol_f(); }
        const real_t  GetVpVol_f(const len_t ir) const {return this->rgrid->GetVpVol_f(ir); }

        const real_t *const* GetVpOverP2AtZero() const { return this->rgrid->GetVpOverP2AtZero(); }
        const real_t *GetVpOverP2AtZero(const len_t ir) const { return this->rgrid->GetVpOverP2AtZero(ir); }

        bool Rebuild(const real_t);
        void RebuildJacobians() { this->rgrid->RebuildJacobians(momentumGrids); }

        real_t Integral(const real_t*) const;
        real_t *IntegralMomentum(const real_t*, real_t *I=nullptr) const;
        real_t IntegralMomentumAtRadius(const len_t, const real_t*) const;
    };
}

#endif/*_DREAM_FVM_GRID_HPP*/

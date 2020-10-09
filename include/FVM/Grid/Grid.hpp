#ifndef _DREAM_FVM_GRID_HPP
#define _DREAM_FVM_GRID_HPP

namespace DREAM::FVM { class Grid; }

#include "FVM/config.h"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/BounceAverager.hpp"

namespace DREAM::FVM {
    class Grid {
    private:

        // Bounce averaged quantities.
        real_t 
            **BA_xi_f1                  = nullptr, // {xi}/xi0 
            **BA_xi_f2                  = nullptr, // {xi}/xi0
            **BA_xi2OverB_f1            = nullptr, // {xi^2(1-xi^2)*Bmin^2/B^2}/(xi0^2(1-xi0^2))
            **BA_xi2OverB_f2            = nullptr, // {xi^2(1-xi^2)*Bmin^2/B^2}/(xi0^2(1-xi0^2))
            **BA_BOverBOverXi_f1        = nullptr, // Theta * sqrt(<B^2>) / (xi0<B/xi>)
            **BA_BOverBOverXi_f2        = nullptr, // Theta * sqrt(<B^2>) / (xi0<B/xi>)
            **BA_B3_f1                  = nullptr, // {B^3}/Bmin^3
            **BA_B3_f2                  = nullptr, // {B^3}/Bmin^3
            **BA_xi2B2_f1               = nullptr, // {xi^2*B^2}/Bmin^2xi0^2
            **BA_xi2B2_f2               = nullptr; // {xi^2*B^2}/Bmin^2xi0^2
        
        // bounce averaged pitch delta function for RP avalanche source
        real_t 
            **avalancheDeltaHat = nullptr,
            **avalancheDeltaHatNegativePitch = nullptr;
        
        // Orbit-averaged metric V'. Size Nr+ x (Np1+ x Np2+).
        real_t
            **Vp    = nullptr,    // Size NR x (N1*N2)
            **Vp_fr = nullptr,    // Size (NR+1) x (N1*N2)
            **Vp_f1 = nullptr,    // Size NR x ((N1+1)*N2)
            **Vp_f2 = nullptr,    // Size NR x (N1*N2)
            **VpOverP2AtZero = nullptr; // Size NR x N2

        // True if phase-space coordinate represents trapped orbit.
        // Size Nr+ x (Np1+ x Np2+).
        bool 
            **isTrapped = nullptr, 
            **isTrapped_fr = nullptr, 
            **isTrapped_f1 = nullptr,
            **isTrapped_f2 = nullptr;

        // If isTrapped, contains poloidal-angle bounce point 
        // theta_b1 or theta_b2, otherwise empty.
        // Size Nr+ x (Np1+ x Np2+).
        real_t  **theta_b1    = nullptr, // on distribution grid 
                **theta_b1_fr = nullptr, // on radial flux grid 
                **theta_b1_f1 = nullptr, // on p1 flux grid
                **theta_b1_f2 = nullptr, // on p2 flux grid
                **theta_b2    = nullptr, // on distribution grid 
                **theta_b2_fr = nullptr, // on radial flux grid 
                **theta_b2_f1 = nullptr, // on p1 flux grid
                **theta_b2_f2 = nullptr; // on p2 flux grid

        void DeallocateVprime();
        void DeallocateBounceParameters();

        void RebuildBounceAveragedQuantities();
        void SetBounceAverage(real_t **&BA_quantity, std::function<real_t(real_t,real_t,real_t,real_t)> F, fluxGridType fluxGridType);
        void DeallocateBAvg();
        void InitializeBAvg(
            real_t **xiAvg_f1, real_t **xiAvg_f2,
            real_t **xi2B2Avg_f1, real_t **xi2B2Avg_f2,
            real_t **B3_f1, real_t **B3_f2,
            real_t **xi2B2_f1, real_t **xi2B2_f2);

    protected:
        BounceAverager *bounceAverager;
        RadialGrid *rgrid;
		MomentumGrid **momentumGrids;

    public:
        Grid(RadialGrid*, MomentumGrid*, const real_t t0=0, 
            FluxSurfaceAverager::quadrature_method qm = FluxSurfaceAverager::QUAD_FIXED_CHEBYSHEV, len_t ntheta = 0);
        ~Grid();

        bool Rebuild(const real_t);
        void RebuildJacobians();

        real_t Integral(const real_t*) const;
        real_t *IntegralMomentum(const real_t*, real_t *I=nullptr) const;
        real_t IntegralMomentumAtRadius(const len_t, const real_t*) const;

        // Returns pointer to the momentum grid with the specified index
        MomentumGrid *GetMomentumGrid(const len_t i) const { return this->momentumGrids[i]; }
        RadialGrid *GetRadialGrid() const { return this->rgrid; }

        /**
         * Grid size getters
         */
        const len_t GetNCells() const;
        const len_t GetNCells_fr() const;
        const len_t GetNCells_f1() const;
        const len_t GetNCells_f2() const;
        const len_t GetNr() const { return this->rgrid->GetNr(); }
        const len_t GetNp1(const len_t ir) const 
            { return this->momentumGrids[ir]->GetNp1(); }
        const len_t GetNp2(const len_t ir) const 
            { return this->momentumGrids[ir]->GetNp2(); }

        /**
         * Getters of bounce-averaged metric
         */
        real_t *const* GetVp() const { return this->Vp; }
        const real_t  *GetVp(const len_t ir) const { return this->Vp[ir]; }
        const real_t GetVp(const len_t ir, const len_t i, const len_t j) const 
            { return Vp[ir][GetNp1(ir)*j+i]; }
        real_t *const* GetVp_fr() const { return this->Vp_fr; }
        const real_t  *GetVp_fr(const len_t ir) const { return this->Vp_fr[ir]; }
        const real_t GetVp_fr(const len_t ir, const len_t i, const len_t j) const 
            { return Vp_fr[ir][GetNp1(0)*j+i]; }
        real_t *const* GetVp_f1() const { return this->Vp_f1; }
        const real_t  *GetVp_f1(const len_t ir) const { return this->Vp_f1[ir]; }
        const real_t GetVp_f1(const len_t ir, const len_t i, const len_t j) const 
            { return Vp_f1[ir][(GetNp1(ir)+1)*j+i]; }
        real_t *const* GetVp_f2() const { return this->Vp_f2; }
        const real_t  *GetVp_f2(const len_t ir) const { return this->Vp_f2[ir]; }
        const real_t GetVp_f2(const len_t ir, const len_t i, const len_t j) const 
            { return Vp_f2[ir][GetNp1(ir)*j+i]; }
        void SetVp(real_t **Vp, real_t **Vp_fr, real_t **Vp_f1, real_t **Vp_f2, real_t **VpOverP2AtZero);
        // (metric * p^2) evaluated at p=0
        const real_t *const* GetVpOverP2AtZero() const { return this->VpOverP2AtZero; }
        const real_t *GetVpOverP2AtZero(const len_t ir) const { return this->VpOverP2AtZero[ir]; }

        /**
         * Getters of flux surface averaged Jacobian
         */
        const real_t *GetVpVol() const {return this->rgrid->GetVpVol(); }
        const real_t  GetVpVol(const len_t ir) const {return this->rgrid->GetVpVol(ir); }
        const real_t *GetVpVol_f() const {return this->rgrid->GetVpVol_f(); }
        const real_t  GetVpVol_f(const len_t ir) const {return this->rgrid->GetVpVol_f(ir); }

        
        /**
         * Getters of isTrapped: true if phase-space point represents a trapped orbit
         */
        const bool IsTrapped(const len_t ir, const len_t i, const len_t j) const 
            {return isTrapped[ir][GetNp1(ir)*j+i];}
        // XXX: Assumes the same momentum grid at all radii 
        const bool IsTrapped_fr(const len_t ir, const len_t i, const len_t j) const 
            {return isTrapped_fr[ir][GetNp1(0)*j+i];}
        const bool IsTrapped_f1(const len_t ir, const len_t i, const len_t j) const 
            {return isTrapped_f1[ir][(GetNp1(ir)+1)*j+i];}
        const bool IsTrapped_f2(const len_t ir, const len_t i, const len_t j) const 
            {return isTrapped_f2[ir][GetNp1(ir)*j+i];}

        /**
         * Getters of lower poloidal-angle bounce points
         */
        const real_t GetThetaBounce1(const len_t ir, const len_t i, const len_t j) const
            {return theta_b1[ir][GetNp1(ir)*j+i];}
        // XXX: Assumes the same momentum grid at all radii 
        const real_t GetThetaBounce1_fr(const len_t ir, const len_t i, const len_t j) const
            {return theta_b1_fr[ir][GetNp1(0)*j+i];}
        const real_t GetThetaBounce1_f1(const len_t ir, const len_t i, const len_t j) const
            {return theta_b1_f1[ir][(GetNp1(ir)+1)*j+i];}
        const real_t GetThetaBounce1_f2(const len_t ir, const len_t i, const len_t j) const
            {return theta_b1_f2[ir][GetNp1(ir)*j+i];}

        /**
         * Getters of upper poloidal-angle bounce points
         */
        const real_t GetThetaBounce2(const len_t ir, const len_t i, const len_t j) const
            {return theta_b2[ir][GetNp1(ir)*j+i];}
        // XXX: Assumes the same momentum grid at all radii 
        const real_t GetThetaBounce2_fr(const len_t ir, const len_t i, const len_t j) const
            {return theta_b2_fr[ir][GetNp1(0)*j+i];}
        const real_t GetThetaBounce2_f1(const len_t ir, const len_t i, const len_t j) const
            {return theta_b2_f1[ir][(GetNp1(ir)+1)*j+i];}
        const real_t GetThetaBounce2_f2(const len_t ir, const len_t i, const len_t j) const
            {return theta_b2_f2[ir][GetNp1(ir)*j+i];}

        /**
         *  Getters for bounce-averaged quantities
         */
        real_t *const* GetBA_xi_f1() const { return this->BA_xi_f1; }
        const real_t  *GetBA_xi_f1(const len_t ir) const { return this->BA_xi_f1[ir]; }
        real_t *const* GetBA_xi_f2() const { return this->BA_xi_f2; }
        const real_t  *GetBA_xi_f2(const len_t ir) const { return this->BA_xi_f2[ir]; }
        real_t *const* GetBA_xi2OverB_f1() const { return this->BA_xi2OverB_f1; }
        const real_t  *GetBA_xi2OverB_f1(const len_t ir) const { return this->BA_xi2OverB_f1[ir]; }
        real_t *const* GetBA_xi2OverB_f2() const { return this->BA_xi2OverB_f2; }
        const real_t  *GetBA_xi2OverB_f2(const len_t ir) const { return this->BA_xi2OverB_f2[ir]; }
        real_t *const* GetBA_BOverBOverXi_f1() const { return this->BA_BOverBOverXi_f1; }
        const real_t  *GetBA_BOverBOverXi_f1(const len_t ir) const { return this->BA_BOverBOverXi_f1[ir]; }
        real_t *const* GetBA_BOverBOverXi_f2() const { return this->BA_BOverBOverXi_f2; }
        const real_t  *GetBA_BOverBOverXi_f2(const len_t ir) const { return this->BA_BOverBOverXi_f2[ir]; }
        real_t *const* GetBA_B3_f1() const { return this->BA_B3_f1; }
        const real_t  *GetBA_B3_f1(const len_t ir) const { return this->BA_B3_f1[ir]; }
        real_t *const* GetBA_B3_f2() const { return this->BA_B3_f2; }
        const real_t  *GetBA_B3_f2(const len_t ir) const { return this->BA_B3_f2[ir]; }
        real_t *const* GetBA_xi2B2_f1() const { return this->BA_xi2B2_f1; }
        const real_t  *GetBA_xi2B2_f1(const len_t ir) const { return this->BA_xi2B2_f1[ir]; }
        real_t *const* GetBA_xi2B2_f2() const { return this->BA_xi2B2_f2; }
        const real_t  *GetBA_xi2B2_f2(const len_t ir) const { return this->BA_xi2B2_f2[ir]; }
        
        const real_t GetAvalancheDeltaHat(const len_t ir, const len_t i, const len_t j, int_t RESign=1)
        {
            len_t pind = GetNp1(ir)*j+i;
            if(RESign>=0)
                return avalancheDeltaHat[ir][pind]; // placeholder for avalanche calculation
            else
                return avalancheDeltaHatNegativePitch[ir][pind];
    }
        void CalculateAvalancheDeltaHat();

        real_t CalculateBounceAverage(len_t ir, len_t i, len_t j, fluxGridType fluxGridType, std::function<real_t(real_t,real_t,real_t,real_t)> F);
        real_t CalculateFluxSurfaceAverage(len_t ir, fluxGridType fluxGridType, std::function<real_t(real_t,real_t,real_t)> F);


        void SetBounceParameters(bool **isTrapped, bool **isTrapped_fr, 
            bool **isTrapped_f1, bool **isTrapped_f2, 
            real_t **theta_b1, real_t **theta_b1_fr, real_t **theta_b1_f1, real_t **theta_b1_f2, 
            real_t **theta_b2, real_t **theta_b2_fr, real_t **theta_b2_f1, real_t **theta_b2_f2 );
    };
}

#endif/*_DREAM_FVM_GRID_HPP*/

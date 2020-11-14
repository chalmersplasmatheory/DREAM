#ifndef _DREAM_FVM_MOMENTUM_GRID_HPP
#define _DREAM_FVM_MOMENTUM_GRID_HPP

namespace DREAM::FVM { class MomentumGrid; }

#include <string>
#include "FVM/Grid/MomentumGridGenerator.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/fluxGridType.enum.hpp"
#include <cmath>

namespace DREAM::FVM {
    class MomentumGrid {
    private:
        len_t np1=0, np2=0;

        // Cell grid coordinate vectors
        real_t *p1, *p2;
        // Flux grid coordinate vectors
        real_t *p1_f, *p2_f;
        // Grid step vectors
        real_t *dp1, *dp2, *dp1_f, *dp2_f;

        real_t 
            *xi0    = nullptr,
            *xi0_f1 = nullptr,
            *xi0_f2 = nullptr,
            *p      = nullptr,
            *p_f1   = nullptr,
            *p_f2   = nullptr,
            *gamma    = nullptr,
            *gamma_f1 = nullptr,
            *gamma_f2 = nullptr;
            

        // Names of coordinates
        std::string p1name, p2name;

        MomentumGridGenerator *generator;

    protected:
        void DeallocateP1();
        void DeallocateP2();
        void DeallocatePAndXi0();

    public:
        MomentumGrid(MomentumGridGenerator *generator, const len_t ir, const RadialGrid *rgrid, const real_t t0=0);
        virtual ~MomentumGrid();

        // Returns the number of cells in this momentum grid
        len_t GetNCells() const { return (np1*np2); }
        len_t GetNCells_f1() const { return ((np1+1)*np2); }
        len_t GetNCells_f2() const { return (np1*(np2+1)); }
        len_t GetNp1() const { return np1; }
        len_t GetNp2() const { return np2; }

        // Grid coordinate getters
        const real_t *GetP1() const { return this->p1; }
        const real_t  GetP1(const len_t i) const { return this->p1[i]; }
        const real_t *GetP2() const { return this->p2; }
        const real_t  GetP2(const len_t i) const { return this->p2[i]; }
        const real_t *GetP1_f() const { return this->p1_f; }
        const real_t  GetP1_f(const len_t i) const { return this->p1_f[i]; }
        const real_t *GetP2_f() const { return this->p2_f; }
        const real_t  GetP2_f(const len_t i) const { return this->p2_f[i]; }
        const real_t *GetDp1() const { return this->dp1; }
        const real_t  GetDp1(const len_t i) const { return this->dp1[i]; }
        const real_t *GetDp2() const { return this->dp2; }
        const real_t  GetDp2(const len_t i) const { return this->dp2[i]; }
        const real_t *GetDp1_f() const { return this->dp1_f; }
        const real_t  GetDp1_f(const len_t i) const { return this->dp1_f[i]; }
        const real_t *GetDp2_f() const { return this->dp2_f; }
        const real_t  GetDp2_f(const len_t i) const { return this->dp2_f[i]; }
        const real_t *GetXi0() const { return this->xi0; }
        const real_t  GetXi0(const len_t i, const len_t j) const { return this->xi0[j*GetNp1()+i]; }
        const real_t *GetXi0_f1() const { return this->xi0_f1; }
        const real_t  GetXi0_f1(const len_t i, const len_t j) const { return this->xi0_f1[j*(GetNp1()+1)+i]; }
        const real_t *GetXi0_f2() const { return this->xi0_f2; }
        const real_t  GetXi0_f2(const len_t i, const len_t j) const { return this->xi0_f2[j*GetNp1()+i]; }
        const real_t *GetP() const { return this->p; }
        const real_t  GetP(const len_t i, const len_t j) const { return this->p[j*GetNp1()+i]; }
        const real_t *GetP_f1() const { return this->p_f1; }
        const real_t  GetP_f1(const len_t i, const len_t j) const { return this->p_f1[j*(GetNp1()+1)+i]; }
        const real_t *GetP_f2() const { return this->p_f2; }
        const real_t  GetP_f2(const len_t i, const len_t j) const { return this->p_f2[j*GetNp1()+i]; }
        const real_t *GetGamma() const { return this->gamma; }
        const real_t  GetGamma(const len_t i, const len_t j) const { return this->gamma[j*GetNp1()+i]; }
        const real_t *GetGamma_f1() const { return this->gamma_f1; }
        const real_t  GetGamma_f1(const len_t i, const len_t j) const { return this->gamma_f1[j*(GetNp1()+1)+i]; }
        const real_t *GetGamma_f2() const { return this->gamma_f2; }
        const real_t  GetGamma_f2(const len_t i, const len_t j) const { return this->gamma_f2[j*GetNp1()+i]; }
        
        const std::string& GetP1Name() const { return this->p1name; }
        const std::string& GetP2Name() const { return this->p2name; }

        virtual bool NeedsRebuild(const real_t t, const bool rGridRebuilt)
        { return this->generator->NeedsRebuild(t, rGridRebuilt); }
        virtual bool Rebuild(const real_t t, const len_t ri, const RadialGrid *rGrid);

        /**
         * Evaluate the metric sqrt(g) on the given poloidal
         * angle grid 'theta' (which contains 'ntheta' grid points)
         * in the given momentum space point.
         *
         * p1:     Value of first momentum coordinate to evaluate metric for.
         * p2:     Value of second momentum coordinate to evaluate metric for.
         * ntheta: Number of points in poloidal angle grid.
         * theta:  Poloidal angle grid.
         * irad:   Index of radial grid point to evaluate metric in.
         * rgrid:  Radial grid to evaluate metric on.
         *
         * RETURNS
         * sqrtg:  Contains the metric upon return (or, rather, sqrt(g)/J)
         */
        
        virtual void EvaluateMetric(
            const len_t i, const len_t j ,
            fluxGridType fluxGridType, 
            const len_t ntheta, const real_t* theta,
            const real_t* BOverBmin, real_t *&sqrtg
        ) const = 0;

        static real_t evaluatePXiMetricOverP2(real_t xi0, real_t BOverBmin){
            if(BOverBmin==1)
                return 2*M_PI;
            real_t xi2_particle = 1 - BOverBmin*(1-xi0*xi0);    
            if (xi2_particle <= 0)
                return 0;
            else {
                return 2*M_PI* BOverBmin * sqrt(xi0*xi0/xi2_particle); 
            }
        }


        // Initialize this momentum grid
        void InitializeP1(
            const std::string& name,
            len_t np1, real_t *p1, real_t *p1_f,
            real_t *dp1, real_t *dp1_f
        ) {
            DeallocateP1();

            this->p1name = name;
            this->np1    = np1;
            this->p1     = p1;
            this->p1_f   = p1_f;
            this->dp1    = dp1;
            this->dp1_f  = dp1_f;
        }

        void InitializeP2(
            const std::string& name,
            len_t np2, real_t *p2, real_t *p2_f,
            real_t *dp2, real_t *dp2_f
        ) {
            DeallocateP2();

            this->p2name = name;
            this->np2    = np2;
            this->p2     = p2;
            this->p2_f   = p2_f;
            this->dp2    = dp2;
            this->dp2_f  = dp2_f;
        }

        void InitializePAndXi0(
            real_t *p, real_t *p_1, real_t *p_2,
            real_t *gamma, real_t *gamma_1, real_t *gamma_2,
            real_t *xi0, real_t *xi01, real_t *xi02
        ) {
            DeallocatePAndXi0();
            this->p      = p;
            this->p_f1   = p_1;
            this->p_f2   = p_2;
            this->gamma  = gamma;
            this->gamma_f1 = gamma_1;
            this->gamma_f2 = gamma_2;
            this->xi0    = xi0;
            this->xi0_f1 = xi01;
            this->xi0_f2 = xi02;
            
        }

    };
}

#endif/*_DREAM_FVM_MOMENTUM_GRID_HPP*/

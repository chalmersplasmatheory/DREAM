#ifndef _DREAM_FVM_RADIAL_GRID_HPP
#define _DREAM_FVM_RADIAL_GRID_HPP

namespace DREAM::FVM { class RadialGrid; }

#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/RadialGridGenerator.hpp"

namespace DREAM::FVM {
	class RadialGrid {
	private:
        len_t nr;

        // Radial grid
        // NOTE that 'r' has 'nr' elements, while
        // 'r_f' has (nr+1) elements.
        real_t *r=nullptr, *r_f=nullptr;
        // Radial grid steps
        //   dr[i]   = r_f[i+1] - r_f[i]   (nr elements)
        //   dr_f[i] = r[i+1] - r[i]       (nr-1 elements)
        real_t *dr=nullptr, *dr_f=nullptr;
        // Jacobian factors
        real_t
            **Vp=nullptr,       // Size NR x (N1*N2)
            **Vp_fr=nullptr,    // Size (NR+1) x (N1*N2)
            **Vp_f1=nullptr,    // Size NR x ((N1+1)*N2)
            **Vp_f2=nullptr;    // Size NR x (N1*N2)

        // Magnetic field quantities
        len_t ntheta;          // Number of poloidal angle points
        real_t *theta=nullptr; // Poloidal angle grid
        real_t
            *B=nullptr,        // Magnetic field strength on r/theta grid (size nr*ntheta)
            *B_f=nullptr;      // Magnetic field strength on r_f/theta grid (size (nr+1)*ntheta)

	protected:
		MomentumGrid **momentumGrids;
        RadialGridGenerator *generator;

    public:
        RadialGrid(RadialGridGenerator*, const real_t t0=0);
        RadialGrid(RadialGridGenerator*, MomentumGrid*, const real_t t0=0);
        ~RadialGrid();

        void DeallocateGrid();
        void DeallocateMagneticField();
        void DeallocateVprime();

        void Initialize(
            real_t *r, real_t *r_f,
            real_t *dr, real_t *dr_f
        ) {
            DeallocateGrid();

            this->r = r;
            this->r_f = r_f;
            this->dr = dr;
            this->dr_f = dr_f;
        }
        void InitializeMagneticField(
            len_t ntheta, real_t *theta,
            real_t *B, real_t *B_f
        ) {
            DeallocateMagneticField();

            this->ntheta = ntheta;
            this->theta  = theta;
            this->B      = B;
            this->B_f    = B_f;
        }
        void InitializeVprime(
            real_t **Vp, real_t **Vp_fr,
            real_t **Vp_f1, real_t **Vp_f2
        ) {
            DeallocateVprime();

            this->Vp = Vp;
            this->Vp_fr = Vp_fr;
            this->Vp_f1 = Vp_f1;
            this->Vp_f2 = Vp_f2;
        }

        bool Rebuild(const real_t);

        // Get number of poloidal angle points
        // Get theta (poloidal angle) grid
        const real_t *GetTheta() const { return this->theta; }
        // Evaluate magnetic field strength at all poloidal angles (on specified flux surface)
        const real_t *BOfTheta() const { return this->B; }
        const real_t *BOfTheta(const len_t ir) const { return this->B+(ir*ntheta); }
        const real_t *BOfTheta_f() const { return this->B_f; }
        const real_t *BOfTheta_f(const len_t ir) const { return this->B_f+(ir*ntheta); }

        // Returns pointer to the momentum grid with the specified index
        MomentumGrid *GetMomentumGrid(const len_t i) { return this->momentumGrids[i]; }
        // Returns the number of radial grid points in this grid
        len_t GetNr() const { return this->nr; }
        // Returns the vector containing all radial grid points
        const real_t *GetR() const { return this->r; }
        const real_t GetR(const len_t i) const { return this->r[i]; }
        const real_t *GetR_f() const { return this->r_f; }
        const real_t GetR_f(const len_t i) const { return this->r_f[i]; }
        
        real_t *const* GetVp() const { return this->Vp; }
        const real_t *GetVp(const len_t ir) const { return this->Vp[ir]; }
        real_t *const* GetVp_fr() const { return this->Vp_fr; }
        const real_t *GetVp_fr(const len_t ir) const { return this->Vp_fr[ir]; }
        real_t *const* GetVp_f1() const { return this->Vp_f1; }
        const real_t *GetVp_f1(const len_t ir) const { return this->Vp_f1[ir]; }
        real_t *const* GetVp_f2() const { return this->Vp_f2; }
        const real_t *GetVp_f2(const len_t ir) const { return this->Vp_f2[ir]; }

        len_t GetNCells() const;
        void SetMomentumGrid(const len_t i, MomentumGrid *m, const real_t t0=0);
        void SetAllMomentumGrids(MomentumGrid*, const real_t t0=0);
	};
}

#endif/*_DREAM_FVM_RADIAL_GRID_HPP*/

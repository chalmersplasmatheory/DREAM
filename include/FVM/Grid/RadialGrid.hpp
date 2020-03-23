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
        real_t *r, *r_f;
        // Radial grid steps
        //   dr[i]   = r_f[i+1] - r_f[i]   (nr elements)
        //   dr_f[i] = r[i+1] - r[i]       (nr-1 elements)
        real_t *dr, *dr_f;
        // Volume enclosed by flux surface
        real_t *volumes;
        // Jacobian factors
        real_t *avGradr2;
        real_t *avGradr2_R2_f;

	protected:
		MomentumGrid **momentumGrids;
        RadialGridGenerator *generator;

    public:
        RadialGrid(RadialGridGenerator*, const real_t t0=0);
        RadialGrid(RadialGridGenerator*, MomentumGrid*, const real_t t0=0);
        ~RadialGrid();

        void Initialize(
            real_t *r, real_t *r_f,
            real_t *dr, real_t *dr_f,
            real_t *V,
            real_t *avGradr2, real_t *avGradr2_R2_f
        ) {
            this->r = r;
            this->r_f = r_f;
            this->dr = dr;
            this->dr_f = dr_f;
            this->volumes = V;
            this->avGradr2 = avGradr2;
            this->avGradr2_R2_f = avGradr2_R2_f;
        }

        bool Rebuild(const real_t);

        // Returns pointer to the momentum grid with the specified index
        MomentumGrid *GetMomentumGrid(const len_t i) { return this->momentumGrids[i]; }
        // Returns the number of radial grid points in this grid
        len_t GetNr() const { return this->nr; }
        // Returns the vector containing all radial grid points
        const real_t *GetR() const { return this->r; }
        const real_t GetR(const len_t i) const { return this->r[i]; }
        const real_t *GetR_f() const { return this->r_f; }
        const real_t GetR_f(const len_t i) const { return this->r_f[i]; }

        const real_t *GetVolumes() const { return this->volumes; }


        len_t GetNCells() const;
        void SetMomentumGrid(const len_t i, MomentumGrid *m, const real_t t0=0);
        void SetAllMomentumGrids(MomentumGrid*, const real_t t0=0);
	};
}

#endif/*_DREAM_FVM_RADIAL_GRID_HPP*/

#ifndef _TQS_FVM_MOMENTUM_GRID_HPP
#define _TQS_FVM_MOMENTUM_GRID_HPP

#include "FVM/Grid/MomentumGridGenerator.hpp"
#include "FVM/Grid/RadialGrid.hpp"

namespace TQS::FVM {
    class MomentumGrid {
    private:
        len_t np1, np2;

        // Cell grid coordinate vectors
        real_t *p1, *p2;
        // Flux grid coordinate vectors
        real_t *p1_f, *p2_f;
        // Grid step vectors
        real_t *dp1, *dp2, *dp1_f, *dp2_f;

        // Coordinate system Jacobian (size np1*np2)
        real_t *J, *J_f;
        // Lam\'{e} coefficients (aka scale factors)
        real_t *h1, *h2, *h1_f, *h2_f;

        MomentumGridGenerator *generator;

    protected:
        // Initialize this momentum grid
        void InitializeP1(
            len_t np1, real_t *p1, real_t *p1_f,
            real_t *dp1, real_t *dp1_f
        ) {
            this->np1   = np1;
            this->p1    = p1;
            this->p1_f  = p1_f;
            this->dp1   = dp1;
            this->dp1_f = dp1_f;
        }

        void InitializeP2(
            len_t np2, real_t *p2, real_t *p2_f,
            real_t *dp2, real_t *dp2_f
        ) {
            this->np2   = np2;
            this->p2    = p2;
            this->p2_f  = p2_f;
            this->dp2   = dp2;
            this->dp2_f = dp2_f;
        }

        void InitializeMetric(
            real_t *volumes, real_t *J, real_t *J_f,
            real_t *h1, real_t *h2, real_t *h1_f, real_t *h2_f
        ) {
            this->volumes = volumes;
            this->J       = J;
            this->J_f     = J_f;
            this->h1      = h1;
            this->h2      = h2;
            this->h1_f    = h1_f;
            this->h2_f    = h2_f;
        }

    public:
        MomentumGrid(MomentumGridGenerator*, const real_t t0=0);
        ~MomentumGrid();

        // Returns the number of cells in this momentum grid
        len_t GetNCells() const { return (np1*np2); }

        // Grid coordinate getters
        const real_t *GetP1() const { return this->p1; }
        const real_t *GetP2() const { return this->p2; }
        const real_t *GetP1_f() const { return this->p1_f; }
        const real_t *GetP2_f() const { return this->p2_f; }
        const real_t *GetDp1() const { return this->dp1; }
        const real_t *GetDp2() const { return this->dp2; }
        const real_t *GetDp1_f() const { return this->dp1_f; }
        const real_t *GetDp2_f() const { return this->dp2_f; }

        virtual bool NeedsRebuild(const real_t t, const bool rGridRebuilt)
        { return this->generator->NeedsRebuild(t, rGridRebuilt); }
        virtual bool Rebuild(const real_t t, const len_t ri, const RadialGrid *rGrid);
    };
}

#endif/*_TQS_FVM_MOMENTUM_GRID_HPP*/

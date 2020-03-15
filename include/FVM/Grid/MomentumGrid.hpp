#ifndef _TQS_FVM_MOMENTUM_GRID_HPP
#define _TQS_FVM_MOMENTUM_GRID_HPP

namespace TQS::FVM { class MomentumGrid; }

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

        // Momentum space "volume"
        real_t *volumes;
        // Lam\'{e} coefficients (aka scale factors)
        real_t *h1, *h2, *h3,
            *h1_f1, *h2_f1, *h3_f1,
            *h1_f2, *h2_f2, *h3_f2;

        MomentumGridGenerator *generator;

    public:
        MomentumGrid(MomentumGridGenerator*, const len_t, const RadialGrid*, const real_t t0=0);
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

        const real_t *GetH1() const { return this->h1; }
        const real_t *GetH2() const { return this->h2; }
        const real_t *GetH3() const { return this->h3; }
        const real_t *GetH1_f1 const { return this->h1_f1; }
        const real_t *GetH2_f1 const { return this->h2_f1; }
        const real_t *GetH3_f1 const { return this->h3_f1; }
        const real_t *GetH1_f2 const { return this->h1_f2; }
        const real_t *GetH2_f2 const { return this->h2_f2; }
        const real_t *GetH3_f2 const { return this->h3_f2; }

        virtual bool NeedsRebuild(const real_t t, const bool rGridRebuilt)
        { return this->generator->NeedsRebuild(t, rGridRebuilt); }
        virtual bool Rebuild(const real_t t, const len_t ri, const RadialGrid *rGrid);

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
            real_t *volumes, real_t *J,
            real_t *h1, real_t *h2, real_t *h3,
            real_t *h1_f1, real_t *h2_f1, real_t *h3_f1,
            real_t *h1_f2, real_t *h2_f2, real_t *h3_f2
        ) {
            this->volumes = volumes;
            this->h1      = h1;
            this->h2      = h2;
            this->h3      = h3;
            this->h1_f1   = h1_f1;
            this->h2_f1   = h2_f1;
            this->h3_f1   = h3_f1;
            this->h1_f2   = h1_f2;
            this->h2_f2   = h2_f2;
            this->h3_f2   = h3_f2;
        }

    };
}

#endif/*_TQS_FVM_MOMENTUM_GRID_HPP*/

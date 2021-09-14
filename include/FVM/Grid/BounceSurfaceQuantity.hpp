
namespace DREAM::FVM { class BounceSurfaceQuantity; }

#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/FluxSurfaceQuantity.hpp"

#ifndef _DREAM_FVM_BOUNCE_SURFACE_QUANTITY_HPP
#define _DREAM_FVM_BOUNCE_SURFACE_QUANTITY_HPP

namespace DREAM::FVM {
    class BounceSurfaceQuantity {
    
    protected:
        real_t
            ***bounceData      = nullptr,
            ***bounceData_fr   = nullptr,
            ***bounceData_f1   = nullptr,
            ***bounceData_f2   = nullptr;

        // Size NR+ x (NP1+ x NP2+).
        // If isTrapped, contains bounce point theta_b1 or theta_b2,
        // otherwise empty.
        real_t  **theta_b1    = nullptr, // on distribution grid 
                **theta_b1_fr = nullptr, // on radial flux grid 
                **theta_b1_f1 = nullptr, // on p1 flux grid
                **theta_b1_f2 = nullptr, // on p2 flux grid
                **theta_b2    = nullptr, // on distribution grid 
                **theta_b2_fr = nullptr, // on radial flux grid 
                **theta_b2_f1 = nullptr, // on p1 flux grid
                **theta_b2_f2 = nullptr; // on p2 flux grid

        FluxSurfaceQuantity *fluxSurfaceQuantity;

        bool passingAllocated = false;
        bool trappedAllocated = false;

        Grid *grid; 
        len_t nr;        
        len_t *np1;
        len_t *np2;

        len_t ntheta_interp_trapped;
        real_t *quad_x_ref;

        gsl_interp_accel *gsl_acc;
        virtual void InterpolateToBounceGrid(real_t ***&bounceData, fluxGridType fluxGridType);
        
        real_t ThetaBounceAtIt(len_t ir, len_t i, len_t j, len_t it, fluxGridType fluxGridType);

        const real_t *GetBounceData(len_t ir, len_t i, len_t j, fluxGridType) const;
        virtual void DeleteData(real_t ***&data, len_t nr, len_t np1, len_t np2, fluxGridType fluxGrid);
    public:
        BounceSurfaceQuantity(Grid *grid, FluxSurfaceQuantity *fluxSurfaceQuantity);
        virtual ~BounceSurfaceQuantity();


        void DeallocateData();
        void AllocateData();
        void AllocateSingle(real_t ***&bounceData, len_t nr, len_t n1, len_t n2);
        
        virtual const real_t *GetData(len_t ir, len_t i, len_t j, fluxGridType) const;
        virtual const real_t evaluateAtTheta(len_t ir, real_t theta, fluxGridType) const;
        
        void SetDataForTrapped(len_t ntheta_interp_trapped, real_t *quad_x_ref);

        static bool IsTrapped(len_t ir, len_t i, len_t j, fluxGridType, Grid*);
        static real_t Theta_B1(len_t ir, len_t i, len_t j, fluxGridType, Grid*);
        static real_t Theta_B2(len_t ir, len_t i, len_t j, fluxGridType, Grid*);
    
        const FluxSurfaceQuantity *GetFluxSurfaceQuantity() const
            {return fluxSurfaceQuantity;}
    
        void SetGridResolution(len_t nr, len_t *np1, len_t *np2)
           {this->nr  = nr;
            this->np1 = np1;
            this->np2 = np2;}
        
    };


}


#endif/*_DREAM_FVM_BOUNCE_SURFACE_QUANTITY_HPP*/


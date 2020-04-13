#ifndef _DREAM_FVM_GRID_MAGNETIC_QUANTITY_HANDLER_HPP
#define _DREAM_FVM_GRID_MAGNETIC_QUANTITY_HANDLER_HPP


namespace DREAM::FVM { class MagneticQuantityHandler; }

#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/config.h"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_integration.h"


namespace DREAM::FVM {
    class MagneticQuantityHandler{

    private:
//        MomentumGrid **momentumGrids;
        len_t nr;

        // ideally B_provided is defined on (r_provided,theta_provided), and we perform 2D-interpolation.
        // for now I assume that they are defined on the radial grid (distribution and flux seperately), as 
        // is output from AnalyticBRadialGridGenerator. In the future, perhaps NumericalBRadialGridGenerator
        // should interpolate to the r and r_f grids and provide to this class.
        // The interpolation method used here should be monotonicity-preserving, so that the min and max of the data are conserved.
        len_t ntheta_interp = 11;
        len_t ntheta_provided;
        real_t *theta_provided;
        real_t **B_provided,
               **B_provided_f,
               **Jacobian_provided,
               **Jacobian_provided_f,
               **ROverR0_provided,
               **ROverR0_provided_f, 
               **NablaR2_provided,
               **NablaR2_provided_f; 

       // Interpolation objects which can evaluate the provided B and J data
       // at any (r,) theta.
       gsl_spline **B_interpolator,
                  **B_interpolator_fr,
                  **Jacobian_interpolator,
                  **Jacobian_interpolator_fr,
                  **ROverR0_interpolator,
                  **ROverR0_interpolator_fr,
                  **NablaR2_interpolator,
                  **NablaR2_interpolator_fr;
                  

       const gsl_integration_fixed_type *thetaGridType = gsl_integration_fixed_legendre; // type of quadrature, default here gauss-legendre
       real_t *x_GL_ref;       // grid points (Gauss-Legendre) on x \in [0,1]
       real_t *weights_GL_ref; // corresponding (Gauss-Legendre) quadrature weights


       real_t *Bmin   = nullptr, // Minimum B on flux surface
              *Bmin_f = nullptr, // Bmin on radial flux grid
              *Bmax   = nullptr, // Maximum B on flux surface
              *Bmax_f = nullptr; // Bmax on radial flux grid


       // Size NR+ x (NP1+ x NP2+). 
       // True if particle is on a trapped orbit
       bool **isTrapped    = nullptr,    // on distribution grid 
            **isTrapped_fr = nullptr,    // on radial flux grid 
            **isTrapped_f1 = nullptr,    // on p1 flux grid 
            **isTrapped_f2 = nullptr;    // on p2 flux grid

       // Size NR+ x (NP1+ x NP2+).
       // If isTrapped, contains bounce point theta_b1 or theta_b2,
       // otherwise empty.
       real_t **theta_b1    = nullptr, // on distribution grid 
              **theta_b1_fr = nullptr, // on radial flux grid 
              **theta_b1_f1 = nullptr, // on p1 flux grid
              **theta_b1_f2 = nullptr, // on p2 flux grid
              **theta_b2    = nullptr, // on distribution grid 
              **theta_b2_fr = nullptr, // on radial flux grid 
              **theta_b2_f1 = nullptr, // on p1 flux grid
              **theta_b2_f2 = nullptr; // on p2 flux grid

       // poloidal angle points and corresponding weights on [0,2*pi] on which 
       // we store theta-dependent magnetic quantities for passing particles
       real_t *theta      = nullptr, // grid points
              *weights    = nullptr; // corresponding quadrature weights

       // Size NR+ x ntheta_interp
       // Magnetic field, Jacobian, R/R0 and |nabla r|^2 on theta grid.
       real_t **B          = nullptr, // on distribution grid
              **B_f        = nullptr, // on radial flux grid
              **Jacobian   = nullptr, // on distribution grid
              **Jacobian_f = nullptr, // on radial flux grid
              **ROverR0    = nullptr, // on distribution grid
              **ROverR0_f  = nullptr, // on radial flux grid
              **NablaR2    = nullptr, // on distribution grid
              **NablaR2_f  = nullptr; // on radial flux grid

       // Major radius of magnetic axis
       //real_t R0;

       // Size NR+ x (NP1+ x NP2+) x ntheta_interp. 
       // If isTrapped, contains theta grid between theta_b1 and theta_b2,
       // otherwise empty.
       real_t ***theta_bounceGrid    = nullptr, // on distribution grid 
              ***theta_bounceGrid_fr = nullptr, // on radial flux grid 
              ***theta_bounceGrid_f1 = nullptr, // on p1 flux grid
              ***theta_bounceGrid_f2 = nullptr; // on p2 flux grid
       // If isTrapped, contains quadrature weights corresponding to theta_bounceGrid.
       real_t ***weights_bounceGrid    = nullptr, // on distribution grid 
              ***weights_bounceGrid_fr = nullptr, // on radial flux grid 
              ***weights_bounceGrid_f1 = nullptr, // on p1 flux grid
              ***weights_bounceGrid_f2 = nullptr; // on p2 flux grid

       // If isTrapped, contains magnetic field evaluated on theta_bounceGrid.
       real_t ***B_bounceGrid    = nullptr, // on distribution grid 
              ***B_bounceGrid_fr = nullptr, // on radial flux grid 
              ***B_bounceGrid_f1 = nullptr, // on p1 flux grid
              ***B_bounceGrid_f2 = nullptr; // on p2 flux grid
       

       // If isTrapped, contains Jacobian evaluated on theta_bounceGrid.
       real_t ***Jacobian_bounceGrid    = nullptr, // on distribution grid 
              ***Jacobian_bounceGrid_fr = nullptr, // on radial flux grid 
              ***Jacobian_bounceGrid_f1 = nullptr, // on p1 flux grid
              ***Jacobian_bounceGrid_f2 = nullptr; // on p2 flux grid

       // Size NR+ x (NP1+ x NP2+) x ntheta_interp.
       // Contains the full metric "J*sqrt(g)" on entire grid.
       real_t ***metricSqrtG    = nullptr,
              ***metricSqrtG_fr = nullptr,
              ***metricSqrtG_f1 = nullptr,
              ***metricSqrtG_f2 = nullptr;

       // Size NR+ x (NP1+ x NP2+) .
       // Contains the bounce-integrated metric VPrime on entire grid.
       real_t **Vp    = nullptr,
              **Vp_fr = nullptr,
              **Vp_f1 = nullptr,
              **Vp_f2 = nullptr;

       // Size NR+
       // Contains the flux surface-integrated spatial Jacobian VPrime_Vol on radial grid.
       real_t *VpVol    = nullptr,
              *VpVol_fr = nullptr;


       

        gsl_interp_accel *gsl_acc  = gsl_interp_accel_alloc();
        virtual bool GetIsTrapped(MomentumGrid *mg, len_t ir, len_t i, len_t j, len_t fluxGridType);
        //virtual void EvaluateVps(MomentumGrid **momentumGrids);
        virtual void EvaluateGrids(MomentumGrid **momentumGrids);
        virtual void SetGrids(MomentumGrid *mg, len_t ir, len_t fluxGridType, bool *isTrapped, 
           real_t *theta_b1, real_t *theta_b2, real_t **theta_bounceGrid, real_t **weights_bounceGrid, real_t **B_bounceGrid, real_t **Jacobian_bounceGrid, real_t **metricSqrtG, real_t *VPrime);
        virtual void SetBounceGrid(MomentumGrid *mg, len_t ir, len_t i, len_t j, len_t fluxGridType,real_t *theta_b1, real_t *theta_b2, real_t *theta_bounceGrid, 
                     real_t *weights_bounceGrid, real_t *B_bounceGrid, real_t *Jacobian_bounceGrid, real_t *metric_bounceGrid);
        virtual void FindBouncePoints(len_t ir, real_t xi0, bool rFluxGrid, real_t *thetab_1, real_t *thetab_2);
        static real_t xiParticleFunction(real_t, void*);
        virtual void FindThetaBounceRoots(real_t *x_lo, real_t *x_up, real_t *root, gsl_function);

        virtual real_t EvaluateBounceIntegral(MomentumGrid *mg, len_t ir, len_t i, len_t j, len_t fluxGridType, std::function<real_t(real_t,real_t)> F);
        virtual real_t EvaluateFluxSurfaceIntegral(len_t ir, bool rFluxGrid, std::function<real_t(real_t,real_t,real_t)> F);

        virtual real_t* GetB(MomentumGrid *mg, len_t ir, len_t i, len_t j, len_t fluxGridType);
        virtual real_t* GetTheta(MomentumGrid *mg, len_t ir, len_t i, len_t j, len_t fluxGridType);
        virtual real_t* GetWeights(MomentumGrid *mg, len_t ir, len_t i, len_t j, len_t fluxGridType);
        virtual real_t* GetMetric(MomentumGrid *mg, len_t ir, len_t i, len_t j, len_t fluxGridType);

        virtual void InitializeGridQuantities(MomentumGrid **momentumGrids);
        virtual void DeallocateGridQuantities(MomentumGrid **momentumGrids);
        
        virtual void InitializeMagneticQuantities();
        virtual void DeallocateMagneticQuantities();

        virtual void InitializeInterpolators();
        virtual void DeallocateInterpolators();
    public:

        MagneticQuantityHandler(len_t ntheta_interp);
        ~MagneticQuantityHandler(); 

        virtual void Initialize(MomentumGrid **momentumGrids,
                      len_t    ntheta_provided,    real_t  *theta_provided, 
                      real_t **B_provided,         real_t **Jacobian_provided,
                      real_t **ROverR0_provided,   real_t **NablaR2_provided,
                      real_t **B_provided_f,       real_t **Jacobian_provided_f,
                      real_t **ROverR0_provided_f, real_t **NablaR2_provided_f);

        virtual real_t CalculateBounceAverage(MomentumGrid *mg, len_t ir, len_t i, len_t j, len_t fluxGridType, std::function<real_t(real_t,real_t)> F);

        virtual real_t CalculateFluxSurfaceAverage(len_t ir, bool rFluxGrid, std::function<real_t(real_t,real_t,real_t)> F);

        virtual real_t GetVp(MomentumGrid *mg, len_t ir, len_t i, len_t j, len_t fluxGridType);
        virtual real_t *GetVp(len_t ir, len_t fluxGridType);
        virtual real_t **GetVp(len_t fluxGridType);

        virtual real_t GetVpVol(len_t ir,bool rFluxGrid);
        virtual real_t *GetVpVol(bool rFluxGrid);
       

    };
}

#endif/*_DREAM_FVM_GRID_MAGNETIC_QUANTITY_HANDLER_HPP*/
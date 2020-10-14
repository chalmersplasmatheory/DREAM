#ifndef _DREAM_FVM_GRID_CYLINDRICAL_RADIAL_GRID_GENERATOR_HPP
#define _DREAM_FVM_GRID_CYLINDRICAL_RADIAL_GRID_GENERATOR_HPP

//#include "FVM/Grid/Grid.hpp"
//#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/RadialGridGenerator.hpp"
#include <functional>

namespace DREAM::FVM {
    class CylindricalRadialGridGenerator : public RadialGridGenerator {
    private:
        real_t xMin=0, xMax=1;
        real_t B0=0;
        real_t *x, *x_f;

        // Set to true when the grid is constructed for the first time
        bool isBuilt = false;

    public:
        CylindricalRadialGridGenerator(const len_t nx, const real_t B0, const real_t x0=0, const real_t xa=1);

        virtual bool NeedsRebuild(const real_t) const override { return (!isBuilt); }
        virtual bool Rebuild(const real_t, RadialGrid*) override;

        virtual real_t JacobianAtTheta(const len_t ir, const real_t) override
            {return 2*M_PI*x[ir];}
        virtual real_t JacobianAtTheta(const len_t ir, const real_t,const real_t,const real_t) override
            {return 2*M_PI*x[ir];}
        virtual real_t ROverR0AtTheta(const len_t, const real_t) override 
            {return 1.0;}
        virtual real_t ROverR0AtTheta(const len_t, const real_t, const real_t, const real_t) override 
            {return 1.0;}
        virtual real_t NablaR2AtTheta(const len_t, const real_t) override 
            {return 1.0;}
        virtual real_t NablaR2AtTheta(const len_t, const real_t, const real_t, const real_t) override 
            {return 1.0;}
        virtual real_t JacobianAtTheta_f(const len_t ir, const real_t) override
            {return 2*M_PI*x_f[ir];}
        virtual real_t JacobianAtTheta_f(const len_t ir, const real_t, const real_t, const real_t) override
            {return 2*M_PI*x_f[ir];}
        virtual real_t ROverR0AtTheta_f(const len_t, const real_t) override 
            {return 1.0;}
        virtual real_t ROverR0AtTheta_f(const len_t, const real_t, const real_t, const real_t) override 
            {return 1.0;}
        virtual real_t NablaR2AtTheta_f(const len_t, const real_t) override 
            {return 1.0;}
        virtual real_t NablaR2AtTheta_f(const len_t, const real_t, const real_t, const real_t) override 
            {return 1.0;}
        virtual void EvaluateGeometricQuantities(const len_t ir, const real_t, real_t &B, real_t &Jacobian, real_t &ROverR0, real_t &NablaR2) override
            {Jacobian=2*M_PI*x[ir]; B=B0; NablaR2 = 1; ROverR0 = 1;}
        virtual void EvaluateGeometricQuantities_fr(const len_t ir, const real_t, real_t &B, real_t &Jacobian, real_t &ROverR0, real_t &NablaR2) override
            {Jacobian=2*M_PI*x_f[ir]; B=B0; NablaR2 = 1; ROverR0 = 1;}
    };
}

#endif/*_DREAM_FVM_GRID_CYLINDRICAL_RADIAL_GRID_GENERATOR_HPP*/

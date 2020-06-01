#ifndef _DREAM_EQUATIONS_COLLISION_QUANTITY_HANDLER_HPP
#define _DREAM_EQUATIONS_COLLISION_QUANTITY_HANDLER_HPP

namespace DREAM { class CollisionQuantityHandler; }
#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/Equations/ParallelDiffusionFrequency.hpp"
#include <gsl/gsl_interp2d.h>

namespace DREAM {
    class CollisionQuantityHandler{
    private:
        const real_t constPreFactor = 4*M_PI
                                *Constants::r0*Constants::r0
                                *Constants::c;
        len_t nr;  // number of radial grid points 
        len_t nZ;  // number of atomic species
        len_t nzs; // number of ion species (including charge states)
        FVM::Grid *grid;
        //EquationSystem *eqSys = nullptr;
        FVM::UnknownQuantityHandler *unknowns = nullptr;
        IonHandler *ionHandler = nullptr;
        enum OptionConstants::momentumgrid_type gridtype;

        CoulombLogarithm *lnLambdaEE;
        CoulombLogarithm *lnLambdaEI;
        SlowingDownFrequency *nuS;
        PitchScatterFrequency *nuD;
        ParallelDiffusionFrequency *nuPar;

        static const len_t  conductivityLenT;
        static const len_t  conductivityLenZ;
        static const real_t conductivityBraams[];
        static const real_t conductivityTmc2[];   // list of T/mc2 
        static const real_t conductivityX[];      // where X = 1/(1+Zeff) 
        

        struct CollisionQuantity::collqty_settings *settings;

        gsl_interp2d *gsl_cond;
        gsl_interp_accel *gsl_xacc;
        gsl_interp_accel *gsl_yacc;

        real_t evaluateElectricalConductivity(len_t i);

    public:

        CollisionQuantityHandler(FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ih,  
                enum OptionConstants::momentumgrid_type mgtype,  struct CollisionQuantity::collqty_settings *cqset);
        ~CollisionQuantityHandler();

        void gridRebuilt();

        void Rebuild();

        SlowingDownFrequency* GetNuS(){return nuS;}
        PitchScatterFrequency* GetNuD(){return nuD;}
        ParallelDiffusionFrequency* GetNuPar(){return nuPar;}
        CoulombLogarithm* GetLnLambdaEE(){return lnLambdaEE;}
        CoulombLogarithm* GetLnLambdaEI(){return lnLambdaEI;}

        const real_t GetLnLambdaC(len_t ir){return lnLambdaEE->GetLnLambdaC(ir);}
        const real_t *GetLnLambdaC(){return lnLambdaEE->GetLnLambdaC();}
        const real_t GetLnLambdaT(len_t ir){return lnLambdaEE->GetLnLambdaT(ir);}
        const real_t *GetLnLambdaT(){return lnLambdaEE->GetLnLambdaT();}
    };

}


#endif/*_DREAM_EQUATIONS_COLLISION_QUANTITY_HANDLER_HPP*/

    



#ifndef _DREAM_EQUATIONS_COLLISION_QUANTITY_HANDLER_HPP
#define _DREAM_EQUATIONS_COLLISION_QUANTITY_HANDLER_HPP

namespace DREAM { class CollisionQuantityHandler; }
#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/Equations/ParallelDiffusionFrequency.hpp"
#include <gsl/gsl_interp2d.h>

namespace DREAM {
    class CollisionQuantityHandler{
    private:
        FVM::UnknownQuantityHandler *unknowns = nullptr;
        IonHandler *ionHandler = nullptr;

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
        
        gsl_interp2d *gsl_cond;
        gsl_interp_accel *gsl_xacc;
        gsl_interp_accel *gsl_yacc;
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

        real_t evaluateElectricalConductivity(len_t i);
    };

}


#endif/*_DREAM_EQUATIONS_COLLISION_QUANTITY_HANDLER_HPP*/

    



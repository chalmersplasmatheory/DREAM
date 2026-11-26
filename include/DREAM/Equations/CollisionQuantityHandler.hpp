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

        CollisionQuantity::collqty_settings *collQtySettings;

        CoulombLogarithm *lnLambdaEE=nullptr;
		CoulombLogarithm *lnLambdaEEhot=nullptr;
        CoulombLogarithm *lnLambdaEI=nullptr;
        SlowingDownFrequency *nuS;
        PitchScatterFrequency *nuD;
        ParallelDiffusionFrequency *nuPar;
    public:
        CollisionQuantityHandler(FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ih,  
                enum OptionConstants::momentumgrid_type mgtype,  CollisionQuantity::collqty_settings *cqset);
        ~CollisionQuantityHandler();

        void gridRebuilt();

        void Rebuild();

        SlowingDownFrequency* GetNuS(){return nuS;}
        PitchScatterFrequency* GetNuD(){return nuD;}
        ParallelDiffusionFrequency* GetNuPar(){return nuPar;}
        CoulombLogarithm* GetLnLambdaEE(){return lnLambdaEE;}
        CoulombLogarithm* GetLnLambdaEEhot(){return lnLambdaEEhot;}
        CoulombLogarithm* GetLnLambdaEI(){return lnLambdaEI;}

        CollisionQuantity::collqty_settings* GetCollisionQuantitySettings()
            {return collQtySettings;}
    };

}


#endif/*_DREAM_EQUATIONS_COLLISION_QUANTITY_HANDLER_HPP*/

    



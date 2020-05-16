#ifndef _DREAM_EQUATIONS_SLOWING_DOWN_FREQUENCY_HPP
#define _DREAM_EQUATIONS_SLOWING_DOWN_FREQUENCY_HPP


#include "FVM/config.h"
#include "CollisionQuantity.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Constants.hpp"
#include "DREAM/Equations/CoulombLogarithm.hpp"
#include <gsl/gsl_math.h>
#include "gsl/gsl_spline.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_interp2d.h>
#include <string>


namespace DREAM {
    class SlowingDownFrequency : public CollisionQuantity {
    private:


        static const len_t  meanExcI_len;
        static const real_t meanExcI_data[];
        static const real_t meanExcI_Zs[];
        static const real_t meanExcI_Z0s[];
        
        CoulombLogarithm *lnLambdaEE;
        real_t *nbound = nullptr;
//        real_t **coulombLogarithmEE; // size nr x np1 (x np2) 
        real_t *hiBethe        = nullptr; // size np1 (x np2). Constant term
        real_t *hiBethe_fr     = nullptr; // size np1 (x np2). Constant term
        real_t *hiBethe_f1     = nullptr; // size np1 (x np2). Constant term
        real_t *hiBethe_f2     = nullptr; // size np1 (x np2). Constant term
        
        real_t **hCold    = nullptr; // size np1 (x np2). Coulomb logarithm dependent + temperature dependent
        real_t **hCold_fr = nullptr; // size np1 (x np2). Coulomb logarithm dependent + temperature dependent
        real_t **hCold_f1 = nullptr; // size np1 (x np2). Coulomb logarithm dependent + temperature dependent
        real_t **hCold_f2 = nullptr; // size np1 (x np2). Coulomb logarithm dependent + temperature dependent
        

        real_t *preFactor     = nullptr;
        real_t *preFactor_fr  = nullptr;
        real_t *preFactor_f1  = nullptr;
        real_t *preFactor_f2  = nullptr;
        real_t *meanExcitationEnergy = nullptr; // size nzs. Constant term
        //len_t *Zs = nullptr;

        real_t ReallyLargeNumber = 1e50; // something that is practically infinite but not quite

        void setPreFactor(real_t *&preFactor, const real_t *pIn, const real_t *gammaIn, len_t np1, len_t np2);
        void setHBethe(real_t *&hiBethe, const real_t *pIn, len_t np1, len_t np2);
        void setHCold(real_t **&hCold, const real_t *pIn, const real_t *gammaIn, len_t nr, len_t np1, len_t np2);

        void GetPartialContributionNonlinear(const real_t lnLc, len_t pind, len_t np1,real_t *&partQty);
        void GetPartialContributionNi(real_t preFactor, real_t *hiBethe, real_t hCold, const real_t lnLee, len_t pind,len_t np1, real_t *&partQty);        
        void evaluatePartialContribution(len_t id_unknown, real_t *preFactor, real_t *hiBethe, real_t **hCold, real_t lnLee, len_t ir, len_t pind, len_t np1, real_t *&partQty);

        void calculateIsotropicNonlinearOperatorMatrix();
        void DeallocateNonlinearMatrix();

        real_t evaluateHColdAtP(len_t ir, real_t p);
        real_t evaluateHBetheAtP(len_t iz, len_t Z0, real_t p);
   protected:
        virtual void AllocateHColdFunc(real_t **&hColdFunc,len_t nr, len_t np1, len_t np2);
        virtual void DeallocateHColdFunc(real_t **&hColdFunc,len_t nr);

        virtual void DeallocatePartialQuantities();
        virtual void AllocatePartialQuantities() override;
        virtual void RebuildPlasmaDependentTerms() override;
        virtual void RebuildConstantTerms() override;


    public:
        SlowingDownFrequency(FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ih,  
                CoulombLogarithm *lnLee,
                enum OptionConstants::momentumgrid_type mgtype,  struct CollisionQuantityHandler::collqtyhand_settings *cqset);
        ~SlowingDownFrequency();

        real_t GetMeanExcitationEnergy(len_t iz, len_t Z0);

        virtual real_t evaluateAtP(len_t ir, real_t p) override;
        virtual void GetPartialContribution(len_t id_unknown, len_t ir, len_t i, len_t j, real_t *&partQty) override;
        virtual void GetPartialContribution_fr(len_t id_unknown, len_t ir, len_t i, len_t j, real_t *&partQty) override;
        virtual void GetPartialContribution_f1(len_t id_unknown, len_t ir, len_t i, len_t j, real_t *&partQty) override;
        virtual void GetPartialContribution_f2(len_t id_unknown, len_t ir, len_t i, len_t j, real_t *&partQty) override;

    };

}

#endif/*_DREAM_EQUATIONS_SLOWING_DOWN_FREQUENCY_HPP*/

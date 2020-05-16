#ifndef _DREAM_EQUATIONS_PITCH_SCATTER_FREQUENCY_HPP
#define _DREAM_EQUATIONS_PITCH_SCATTER_FREQUENCY_HPP


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
    class PitchScatterFrequency : public CollisionQuantity {
    private:


        static const len_t  ionSizeAj_len; 
        static const real_t ionSizeAj_data[];       
        static const real_t ionSizeAj_Zs[];
        static const real_t ionSizeAj_Z0s[];
        
        CoulombLogarithm *lnLambdaEI;
        CoulombLogarithm *lnLambdaEE;
        real_t *nbound = nullptr;
//        real_t **coulombLogarithmEE; // size nr x np1 (x np2) 
        real_t *gKirillov        = nullptr; // size np1 (x np2). Constant term
        real_t *gKirillov_fr     = nullptr; // size np1 (x np2). Constant term
        real_t *gKirillov_f1     = nullptr; // size np1 (x np2). Constant term
        real_t *gKirillov_f2     = nullptr; // size np1 (x np2). Constant term
        
        real_t **gCold    = nullptr; // size np1 (x np2). Coulomb logarithm dependent + temperature dependent
        real_t **gCold_fr = nullptr; // size np1 (x np2). Coulomb logarithm dependent + temperature dependent
        real_t **gCold_f1 = nullptr; // size np1 (x np2). Coulomb logarithm dependent + temperature dependent
        real_t **gCold_f2 = nullptr; // size np1 (x np2). Coulomb logarithm dependent + temperature dependent
        

        real_t *preFactor     = nullptr;
        real_t *preFactor_fr  = nullptr;
        real_t *preFactor_f1  = nullptr;
        real_t *preFactor_f2  = nullptr;
        real_t *ionEffectiveSize = nullptr; // size nzs. Constant term
        //len_t *Zs = nullptr;

        real_t ReallyLargeNumber = 1e50; // something that is practically infinite but not quite

        void setPreFactor(real_t *&preFactor, const real_t *pIn, const real_t *gammaIn, len_t np1, len_t np2);
        void setGKirillov(real_t *&gKirillov, const real_t *pIn, len_t np1, len_t np2);
        void setGCold(real_t **&gCold, const real_t *pIn, len_t nr, len_t np1, len_t np2);

        
        void GetPartialContributionNonlinear(const real_t lnLc, len_t pind,len_t np1, real_t *&partQty);
        void GetPartialContributionNi(real_t preFactor, real_t *gKirillov, real_t gCold, const real_t lnLei, const real_t lnLee, len_t pind, len_t np1,real_t *&partQty);        
        
        void GetPartialContributionNCold(real_t gCold, real_t preFactor, real_t lnLee, real_t *&partQty);
        void evaluatePartialContribution(len_t id_unknown, real_t *preFactor, real_t *gKirillov, real_t **gCold, real_t lnLei, real_t lnLee,len_t ir, len_t pind, len_t np1,real_t *&partQty);

        void calculateIsotropicNonlinearOperatorMatrix();
        void DeallocateNonlinearMatrix();
        real_t evaluateGColdAtP(len_t ir, real_t p);

        real_t evaluateGKirillovAtP(len_t iz, len_t Z0, real_t p);
   protected:
        virtual void AllocateGColdFunc(real_t **&gColdFunc,len_t nr, len_t np1, len_t np2);
        virtual void DeallocateGColdFunc(real_t **&gColdFunc,len_t nr);

        virtual void DeallocatePartialQuantities();
        virtual void AllocatePartialQuantities() override;
        virtual void RebuildPlasmaDependentTerms() override;
        virtual void RebuildConstantTerms() override;

        virtual void GetPartialContribution(len_t id_unknown, len_t ir, len_t i, len_t j, real_t *&partQty) override;
        virtual void GetPartialContribution_fr(len_t id_unknown, len_t ir, len_t i, len_t j, real_t *&partQty) override;
        virtual void GetPartialContribution_f1(len_t id_unknown, len_t ir, len_t i, len_t j, real_t *&partQty) override;
        virtual void GetPartialContribution_f2(len_t id_unknown, len_t ir, len_t i, len_t j, real_t *&partQty) override;

    public:
        PitchScatterFrequency(FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ih,  
                CoulombLogarithm *lnLei,CoulombLogarithm *lnLee,
                enum OptionConstants::momentumgrid_type mgtype,  struct CollisionQuantityHandler::collqtyhand_settings *cqset);
        ~PitchScatterFrequency();

        real_t GetIonEffectiveSizeAj(len_t iz, len_t Z0);

        virtual real_t evaluateAtP(len_t ir, real_t p) override;

    };

}

#endif/*_DREAM_EQUATIONS_PITCH_SCATTER_FREQUENCY_HPP*/

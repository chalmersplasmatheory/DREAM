#ifndef _DREAM_EQUATIONS_COULOMB_LOGARITHM_HPP
#define _DREAM_EQUATIONS_COULOMB_LOGARITHM_HPP


#include "FVM/config.h"
#include "CollisionQuantity.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Constants.hpp"
#include <gsl/gsl_math.h>
#include "gsl/gsl_spline.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_interp2d.h>
#include <string>


namespace DREAM {
    class CoulombLogarithm : public CollisionQuantity {
    private:
        real_t *lnLambda_c = nullptr;
        real_t *lnLambda_T = nullptr;

        bool isLnEE = false;
        bool isLnEI = false;
        
        void AssembleConstantLnLambda(); 
        void AssembleConstantLnLambda(real_t **&lnLambda, len_t nr, len_t np1, len_t np2);
        void AssembleWithPXiGrid();
        void AssembleWithPXiGrid(real_t **&lnLambda,const real_t *pVec, len_t nr, len_t np1, len_t np2);
        void AssembleWithGeneralGrid();
        void AssembleWithGeneralGrid(real_t **&lnLambda,const real_t *pVec, len_t nr, len_t np1, len_t np2);
        
        void DeallocatePartialQuantities();
    protected:
        void AssembleQuantity() override;
        virtual void AllocatePartialQuantities() override;
        virtual void RebuildPlasmaDependentTerms() override;
        virtual void RebuildConstantTerms() override{return;};


    public:

        CoulombLogarithm(FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ih,  
                enum OptionConstants::momentumgrid_type mgtype,  struct CollisionQuantityHandler::collqtyhand_settings *cqset,
                len_t lnLambdaType);
        ~CoulombLogarithm();
        
        const real_t GetLnLambdaC(const len_t ir) const {return lnLambda_c[ir];}
        const real_t  *GetLnLambdaC() const{return lnLambda_c;}

        const real_t GetLnLambdaT(const len_t ir) const {return lnLambda_T[ir];}
        const real_t  *GetLnLambdaT() const{return lnLambda_T;}

        virtual real_t evaluateAtP(len_t ir, real_t p) override;

        virtual void GetPartialContribution(len_t,len_t,len_t,len_t,real_t*&) override {return;}
        virtual void GetPartialContribution_fr(len_t,len_t,len_t,len_t, real_t*&) override {return;}
        virtual void GetPartialContribution_f1(len_t,len_t,len_t,len_t, real_t*&) override {return;}
        virtual void GetPartialContribution_f2(len_t,len_t,len_t,len_t, real_t*&) override {return;}

    };
}


#endif/*_DREAM_EQUATIONS_COULOMB_LOGARITHM_HPP*/

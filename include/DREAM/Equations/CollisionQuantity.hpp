#ifndef _DREAM_EQUATIONS_COLLISION_QUANTITY_HPP
#define _DREAM_EQUATIONS_COLLISION_QUANTITY_HPP

namespace DREAM { class CollisionQuantity; }

#include "FVM/config.h"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Constants.hpp"
#include "DREAM/Equations/CollisionQuantityHandler.hpp"

namespace DREAM {
    class CollisionQuantity{
        
    private:
        void AssembleQuantity(real_t **&collisionQuantity, len_t nr, len_t np1, len_t np2, len_t fluxGridType);

    protected:
        // XXX we assume explicitly that CollisionQuantities have
        // the same MomentumGrid at all radii 
        FVM::MomentumGrid *mg;
        FVM::RadialGrid *rGrid;
        bool isPXiGrid;
        bool isNonlinear;
        bool isNonScreened;
        bool isPartiallyScreened;
        struct CollisionQuantityHandler::collqtyhand_settings *collQtySettings; 
        IonHandler *ionHandler;
        FVM::UnknownQuantityHandler *unknowns;

        len_t id_ncold, id_ni, id_Tcold, id_fhot;

        len_t np1, np2, nr, nzs, nZ;
        len_t np2_store;
        real_t kInterpolate = 5;

        const real_t constPreFactor = 4*M_PI
                                *Constants::r0*Constants::r0
                                *Constants::c;
        
        bool gridRebuilt = true;

        // Set buildOnlyF1F2 to false if we need to build collision frequencies 
        // on the distribution grid or radial flux grid
        bool buildOnlyF1F2 = true;
        

        virtual void AllocatePartialQuantities()=0;
        
        virtual void RebuildPlasmaDependentTerms()=0;
        virtual void RebuildConstantTerms()=0;
        virtual void AssembleQuantity();
        
        real_t **collisionQuantity    = nullptr;
        real_t **collisionQuantity_fr = nullptr;
        real_t **collisionQuantity_f1 = nullptr;
        real_t **collisionQuantity_f2 = nullptr;
        
        real_t **nonlinearMat = nullptr;
        real_t *trapzWeights = nullptr;

        void AllocateCollisionQuantity(real_t **&cty, len_t nr, len_t np1, len_t np2);
        void AllocateCollisionQuantities();
        void DeallocateCollisionQuantity(real_t **&collisionQuantity, len_t nr);
        void DeallocateCollisionQuantities();

        bool parametersHaveChanged();
        real_t evaluatePrefactorAtP(real_t p, real_t gamma) {return constPreFactor * gamma*gamma/(p*p*p);}
        
        virtual real_t evaluatePsi0(len_t ir, real_t p);
        virtual real_t evaluatePsi1(len_t ir, real_t p);
        static real_t psi0Integrand(real_t s, void *params);
        static real_t psi1Integrand(real_t s, void *params);
        virtual real_t evaluateExp1OverThetaK(real_t Theta, real_t n);
        void InitializeGSLWorkspace();
        void DeallocateGSL();
        gsl_integration_fixed_workspace **gsl_w = nullptr;
        virtual void GetPartialContribution(len_t id_unknown, len_t ir, len_t i, len_t j, real_t *&partQty) = 0;
        virtual void GetPartialContribution_fr(len_t id_unknown, len_t ir, len_t i, len_t j, real_t *&partQty) = 0;
        virtual void GetPartialContribution_f1(len_t id_unknown, len_t ir, len_t i, len_t j, real_t *&partQty) = 0;
        virtual void GetPartialContribution_f2(len_t id_unknown, len_t ir, len_t i, len_t j, real_t *&partQty) = 0;

    public: 

        CollisionQuantity(FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ih,  
                enum OptionConstants::momentumgrid_type mgtype,  struct CollisionQuantityHandler::collqtyhand_settings *cqset);
        ~CollisionQuantity();

        void Rebuild();
        void GridRebuilt(){gridRebuilt=true;}
        void AddNonlinearContribution(const real_t *lnL);
       
        const real_t  GetValue(const len_t ir, const len_t i, const len_t j) const 
            {return this->collisionQuantity[ir][np1*j+i]; }
        const real_t  *GetValue(const len_t ir) const 
            { return this->collisionQuantity[ir]; }
        real_t *const* GetValue() const 
        { return this->collisionQuantity; }

        const real_t  GetValue_fr(const len_t ir, const len_t i, const len_t j) const 
            {return this->collisionQuantity_fr[ir][np1*j+i]; }
        const real_t  *GetValue_fr(const len_t ir) const 
            { return this->collisionQuantity_fr[ir]; }
        real_t *const* GetValue_fr() const 
        { return this->collisionQuantity_fr; }

        const real_t  GetValue_f1(const len_t ir, const len_t i, const len_t j) const 
            {return this->collisionQuantity_f1[ir][(np1+1)*j+i]; }
        const real_t  *GetValue_f1(const len_t ir) const 
            { return this->collisionQuantity_f1[ir]; }
        real_t *const* GetValue_f1() const 
        { return this->collisionQuantity_f1; }

        const real_t  GetValue_f2(const len_t ir, const len_t i, const len_t j) const 
            {return this->collisionQuantity_f2[ir][np1*j+i]; }
        const real_t  *GetValue_f2(const len_t ir) const 
            { return this->collisionQuantity_f2[ir]; }
        real_t *const* GetValue_f2() const 
        { return this->collisionQuantity_f2; }

        virtual real_t evaluateAtP(len_t ir, real_t p) = 0;
        real_t *GetUnknownPartialContribution(len_t id_unknown, len_t ir, len_t i, len_t j, len_t fluxGridMode, real_t *&partQty);
        
        

    };

}


#endif/*_DREAM_EQUATIONS_COLLISION_QUANTITY_HPP*/

    

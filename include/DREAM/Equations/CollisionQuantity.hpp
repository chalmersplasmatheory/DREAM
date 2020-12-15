#ifndef _DREAM_EQUATIONS_COLLISION_QUANTITY_HPP
#define _DREAM_EQUATIONS_COLLISION_QUANTITY_HPP

namespace DREAM { class CollisionQuantity; }

#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Constants.hpp"

namespace DREAM {
    class CollisionQuantity{
    public:
        enum LnLambdaType {
            LNLAMBDATYPE_EE,
            LNLAMBDATYPE_EI
        };

        struct collqty_settings {
            enum OptionConstants::collqty_collfreq_type 
                        collfreq_type   = OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_PARTIALLY_SCREENED;
            enum OptionConstants::collqty_collfreq_mode 
                        collfreq_mode   = OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_SUPERTHERMAL;
            enum OptionConstants::collqty_lnLambda_type 
                        lnL_type        = OptionConstants::COLLQTY_LNLAMBDA_ENERGY_DEPENDENT;
            enum OptionConstants::uqty_n_cold_eqn       
                        ncold_type      = OptionConstants::UQTY_N_COLD_EQN_PRESCRIBED;
            enum OptionConstants::eqterm_nonlinear_mode
                        nonlinear_mode  = OptionConstants::EQTERM_NONLINEAR_MODE_NEGLECT;
            enum OptionConstants::eqterm_bremsstrahlung_mode 
                        bremsstrahlung_mode = OptionConstants::EQTERM_BREMSSTRAHLUNG_MODE_NEGLECT;
            enum OptionConstants::collqty_pstar_mode
                        pstar_mode = OptionConstants::COLLQTY_PSTAR_MODE_COLLISIONLESS;
            enum OptionConstants::collqty_screened_diffusion_mode
                        screened_diffusion = OptionConstants::COLLQTY_SCREENED_DIFFUSION_MODE_MAXWELLIAN;
        };
        
    private:
        void AssembleQuantity();
        void AllocateCollisionQuantity(real_t **&cty, len_t nr, len_t np1, len_t np2);
        void AllocateCollisionQuantities();
        void DeallocateCollisionQuantity(real_t **&collisionQuantity, len_t nr);
        bool parametersHaveChanged();
        
    protected:
        bool gridRebuilt = true;
        // XXX we assume explicitly that CollisionQuantities have
        // the same MomentumGrid at all radii 
        FVM::MomentumGrid *mg;
        FVM::RadialGrid *rGrid;
        bool isPXiGrid;
        bool isNonlinear;
        bool isNonScreened;
        bool isPartiallyScreened;
        bool isBrems;
        IonHandler *ionHandler;
        FVM::UnknownQuantityHandler *unknowns;
        collqty_settings *collQtySettings;

        len_t id_ncold, id_ni, id_Tcold;
        len_t np1, np2, nr, nzs, nZ, np2_store;
        real_t kInterpolate;

        const real_t constPreFactor = 4*M_PI
                                *Constants::r0*Constants::r0
                                *Constants::c;
        bool buildOnlyF1F2;
        
        virtual void AllocatePartialQuantities()=0;
        
        virtual void RebuildPlasmaDependentTerms()=0;
        virtual void RebuildConstantTerms()=0;
        virtual void AssembleQuantity(real_t **&collisionQuantity, len_t nr, len_t np1, len_t np2, FVM::fluxGridType) = 0;
        void DeallocateCollisionQuantities();
        
        real_t **collisionQuantity    = nullptr;
        real_t **collisionQuantity_fr = nullptr;
        real_t **collisionQuantity_f1 = nullptr;
        real_t **collisionQuantity_f2 = nullptr;

    public: 

        CollisionQuantity(FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ih,  
                enum OptionConstants::momentumgrid_type mgtype,  struct collqty_settings *cqset);
        virtual ~CollisionQuantity();
        void Rebuild();
        void GridRebuilt(){gridRebuilt=true;}
       
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

        real_t *const* GetValue(FVM::fluxGridType fluxGridType) const 
        { 
            switch(fluxGridType){
                case FVM::FLUXGRIDTYPE_DISTRIBUTION:
                    return this->collisionQuantity; 
                case FVM::FLUXGRIDTYPE_RADIAL:
                    return this->collisionQuantity_fr; 
                case FVM::FLUXGRIDTYPE_P1:
                    return this->collisionQuantity_f1; 
                case FVM::FLUXGRIDTYPE_P2:
                    return this->collisionQuantity_f2;
                default:
                    return nullptr; 
            }
        }

        real_t evaluateAtP(len_t ir, real_t p);
        virtual real_t evaluateAtP(len_t ir, real_t p, struct collqty_settings *inSettings) = 0;

        real_t evaluatePartialAtP(len_t ir, real_t p, len_t derivId, len_t n);
        virtual real_t evaluatePartialAtP(len_t ir, real_t p, len_t derivId, len_t n,struct collqty_settings *inSettings) = 0;

        const collqty_settings *GetSettings() const{return collQtySettings;}
    };
}

#endif/*_DREAM_EQUATIONS_COLLISION_QUANTITY_HPP*/   

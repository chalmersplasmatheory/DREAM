#ifndef _DREAM_EQUATIONS_RUNAWAY_FLUID_HPP
#define _DREAM_EQUATIONS_RUNAWAY_FLUID_HPP

namespace DREAM { class RunawayFluid; }

#include "DREAM/Equations/SlowingDownFrequency.hpp"
#include "DREAM/Equations/PitchScatterFrequency.hpp"
#include <algorithm>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>

namespace DREAM {
    class RunawayFluid{
    private:
        const real_t constPreFactor = 4*M_PI
                                *Constants::r0*Constants::r0
                                *Constants::c;
        static const real_t tritiumHalfLife;
        static const real_t tritiumDecayEnergyEV;
        
        FVM::RadialGrid *rGrid;
        FVM::UnknownQuantityHandler *unknowns;
        SlowingDownFrequency *nuS;
        PitchScatterFrequency *nuD;
        CoulombLogarithm *lnLambdaEE;
        CoulombLogarithm *lnLambdaEI;
        len_t nr;
        CollisionQuantity::collqty_settings *collQtySettings;

        gsl_integration_workspace *gsl_ad_w;
        const gsl_root_fsolver_type *GSL_rootsolver_type;
        gsl_root_fsolver *fsolve;
        const gsl_min_fminimizer_type *fmin_type;
        gsl_min_fminimizer *fmin;

        CollisionQuantity::collqty_settings *collSettingsForEc;
        CollisionQuantity::collqty_settings *collSettingsForPc;

        len_t id_ncold;
        len_t id_ntot;
        len_t id_ni;
        len_t id_Tcold;
        len_t id_Eterm;

        real_t *ncold;
        real_t *ntot;
        real_t *Tcold;
        real_t *Eterm;

        real_t *Ec_free=nullptr;                 // Connor-Hastie field with only bound
        real_t *Ec_tot=nullptr;                  // Connor-Hastie field with free+bound
        real_t *tauEERel=nullptr;                // Relativistic electron collision time
        real_t *tauEETh=nullptr;                 // Thermal electron collision time
        
        real_t *EDreic=nullptr;                  // Dreicer field
        real_t *criticalREMomentum=nullptr;      // Critical momentum for runaway p_star 
        real_t *criticalREMomentumInvSq=nullptr; // Inverse square p_star
        real_t *pc_COMPLETESCREENING = nullptr;
        real_t *pc_NOSCREENING = nullptr;
        real_t *avalancheGrowthRate=nullptr;     // (dnRE/dt)_ava = nRE*Gamma_ava
        real_t *tritiumRate=nullptr;             // (dnRE/dt)_Tritium = nTritium * ...
        real_t *comptonRate=nullptr;             // (dnRE/dt)_Compton = n_tot * ...
        real_t *effectiveCriticalField=nullptr;  // Eceff: Gamma_ava(Eceff) = 0

        bool gridRebuilt;
        bool parametersHaveChanged();

        void AllocateQuantities();
        void DeallocateQuantities();
        
        void CalculateDerivedQuantities();
        void CalculateEffectiveCriticalField(bool useApproximateMethod);
        void CalculateCriticalMomentum();
        void CalculateGrowthRates();


        static void FindECritInterval(len_t ir, real_t *E_lower, real_t *E_upper, void *par);
        static void FindPExInterval(real_t *p_ex_guess, real_t *p_ex_lower, real_t *p_ex_upper, void *par, real_t p_upper_threshold);
        static void FindRoot(real_t x_lower, real_t x_upper, real_t *root, gsl_function gsl_func, gsl_root_fsolver *s);
        static void FindInterval(real_t *x_lower, real_t *x_upper, gsl_function gsl_func );

        real_t BounceAverageFunc(len_t ir, std::function<real_t(real_t,real_t)> Func);

        static real_t FindUExtremumAtE(real_t Eterm, void *par);
        static real_t UAtPFunc(real_t p, void *par);
        
        
        static real_t pStarFunction(real_t, void *);
        real_t evaluateBarNuSNuDAtP(len_t ir, real_t p, CollisionQuantity::collqty_settings *inSettings);        
        real_t evaluateNuDHat(len_t ir, real_t p, CollisionQuantity::collqty_settings *inSettings);
        real_t evaluateNuSHat(len_t ir, real_t p, CollisionQuantity::collqty_settings *inSettings);
        
//        real_t evaluateNuSNuDTerm(len_t ir, real_t p, OptionConstants::collqty_collfreq_type collfreq_type);

        static const len_t  conductivityLenT;
        static const len_t  conductivityLenZ;
        static const real_t conductivityBraams[];
        static const real_t conductivityTmc2[];   // list of T/mc2 
        static const real_t conductivityX[];      // where X = 1/(1+Zeff) 
        
        gsl_interp2d *gsl_cond;
        gsl_interp_accel *gsl_xacc;
        gsl_interp_accel *gsl_yacc;


    protected:
    public:
        RunawayFluid(FVM::Grid *g, FVM::UnknownQuantityHandler *u, SlowingDownFrequency *nuS, 
        PitchScatterFrequency *nuD, CoulombLogarithm *lnLEE,CoulombLogarithm *lnLEI, CollisionQuantity::collqty_settings *cqs);
        ~RunawayFluid();

        real_t testEvalU(len_t ir, real_t p, real_t Eterm, bool useApproximateMethod, CollisionQuantity::collqty_settings *inSettings);

        real_t evaluateAnalyticPitchDistribution(len_t ir, real_t xi0, real_t p, real_t Eterm,CollisionQuantity::collqty_settings *inSettings, gsl_integration_workspace *gsl_ad_w);
        real_t evaluateApproximatePitchDistribution(len_t ir, real_t xi0, real_t p, real_t Eterm,CollisionQuantity::collqty_settings *inSettings);
        real_t evaluatePitchDistribution(len_t ir, real_t xi0, real_t p, real_t Eterm, CollisionQuantity::collqty_settings *inSettings, gsl_integration_workspace *gsl_ad_w, bool useApproximatePitchDistribution);

        static real_t evaluateTritiumRate(real_t gamma_c);
        static real_t evaluateComptonRate(real_t pc,gsl_integration_workspace *gsl_ad_w);
        static real_t evaluateComptonPhotonFluxSpectrum(real_t Eg);
        static real_t evaluateComptonTotalCrossSectionAtP(real_t Eg, real_t pc);


        void Rebuild(bool useApproximateMethod);
        void GridRebuilt();
        const real_t GetEffectiveCriticalField(len_t ir) const
            {return effectiveCriticalField[ir];}
        const real_t* GetEffectiveCriticalField() const
            {return effectiveCriticalField;}
        
        const real_t GetDreicerElectricField(len_t ir) const
            {return EDreic[ir];}
        const real_t* GetDreicerElectricField() const
            {return EDreic;}
        
        const real_t GetConnorHastieField_COMPLETESCREENING(len_t ir) const
            {return Ec_free[ir];}
        const real_t* GetConnorHastieField_COMPLETESCREENING() const
            {return Ec_free;}
        
        const real_t GetConnorHastieField_NOSCREENING(len_t ir) const
            {return Ec_tot[ir];}
        const real_t* GetConnorHastieField_NOSCREENING() const
            {return Ec_tot;}
        
        const real_t GetElectronCollisionTimeRelativistic(len_t ir) const
            {return tauEERel[ir];}
        const real_t* GetElectronCollisionTimeRelativistic() const
            {return tauEERel;}
        
        const real_t GetElectronCollisionTimeThermal(len_t ir) const
            {return tauEETh[ir];}
        const real_t* GetElectronCollisionTimeThermal() const
            {return tauEETh;}
        
        

        const real_t GetAvalancheGrowthRate(len_t ir) const
            {return avalancheGrowthRate[ir];}
        const real_t* GetAvalancheGrowthRate() const
            {return avalancheGrowthRate;}
        
        const real_t GetTritiumRunawayRate(len_t ir) const
            {return tritiumRate[ir];}
        const real_t* GetTritiumRunawayRate() const
            {return tritiumRate;}
        
        const real_t GetComptonRunawayRate(len_t ir) const
            {return comptonRate[ir];}
        const real_t* GetComptonRunawayRate() const
            {return comptonRate;}

        const real_t GetEffectiveCriticalRunawayMomentum(len_t ir) const
            {return criticalREMomentum[ir];}
        const real_t* GetEffectiveCriticalRunawayMomentum() const
            {return criticalREMomentum;}
        

        const CollisionQuantity::collqty_settings *GetSettings() const{return collQtySettings;}
        CoulombLogarithm* GetLnLambda(){return lnLambdaEE;}

        real_t evaluateNeoclassicalConductivityCorrection(len_t ir, real_t Zeff, bool collisionless = false);

        real_t evaluateSauterElectricConductivity(len_t ir, real_t Zeff);
        real_t evaluateBraamsElectricConductivity(len_t ir, real_t Zeff);

        /**
         * Placeholder calculation of the partial derivative of conductivity
         * with respect to temperature; assumes for now that it has 
         * a pure 1/T^1.5 dependence.
         */  
        real_t* evaluatePartialContributionSauterConductivity(real_t *Zeff, len_t derivId) {
            real_t *dSigma = new real_t[nr];
            if(derivId!=id_Tcold)
                for(len_t ir = 0; ir<nr; ir++)
                    dSigma[ir] = 0;
            else { 
                real_t *Tcold = unknowns->GetUnknownData(id_Tcold);
                for(len_t ir = 0; ir<nr; ir++)
                    dSigma[ir] = -1.5 * evaluateSauterElectricConductivity(ir,Zeff[ir]) / Tcold[ir];
            }
            return dSigma; 
        }

        /**
         * Placeholder (?) calculation of the partial derivative of the 
         * conductivity with respect to temperature; assumes for now that 
         * it has a pure T^1.5 dependence.
         */  
        real_t* evaluatePartialContributionBraamsConductivity(real_t *Zeff, len_t derivId) {
            real_t *dSigma = new real_t[nr];
            if(derivId!=id_Tcold)
                for(len_t ir = 0; ir<nr; ir++)
                    dSigma[ir] = 0;
            else { 
                real_t *Tcold = unknowns->GetUnknownData(id_Tcold);
                for(len_t ir = 0; ir<nr; ir++)
                    dSigma[ir] = 1.5 * evaluateBraamsElectricConductivity(ir,Zeff[ir]) / Tcold[ir];
            }
            return dSigma; 
        }

    };

}


#endif/*_DREAM_EQUATIONS_RUNAWAY_FLUID_HPP*/

    



#ifndef _DREAM_EQUATIONS_RUNAWAY_FLUID_HPP
#define _DREAM_EQUATIONS_RUNAWAY_FLUID_HPP

namespace DREAM { class RunawayFluid; }

#include <algorithm>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <string>
#include "DREAM/Equations/AnalyticDistributionRE.hpp"
#include "DREAM/Equations/ConnorHastie.hpp"
#include "DREAM/Equations/DreicerNeuralNetwork.hpp"
#include "DREAM/Equations/EffectiveCriticalField.hpp"
#include "DREAM/Equations/PitchScatterFrequency.hpp"
#include "DREAM/Equations/SlowingDownFrequency.hpp"
#include "DREAM/IonHandler.hpp"
#include "FVM/Interpolator1D.hpp"
#include "FVM/TimeKeeper.hpp"

namespace DREAM {
    class RunawayFluid {
    private:
        const real_t constPreFactor = 
            4*M_PI*Constants::r0*Constants::r0*Constants::c;
        static const real_t tritiumHalfLife;
        static const real_t tritiumDecayEnergyEV;
        
        FVM::RadialGrid *rGrid;
        SlowingDownFrequency *nuS;
        PitchScatterFrequency *nuD;
        CoulombLogarithm *lnLambdaEE;
        bool extrapolateDreicer;
        CoulombLogarithm *lnLambdaEI;
        len_t nr;
        FVM::UnknownQuantityHandler *unknowns;
        IonHandler *ions;
        AnalyticDistributionRE *analyticRE;      // analytic distribution of runaway electrons 

        CollisionQuantity::collqty_settings *collSettingsForPc;
        CollisionQuantity::collqty_settings *collSettingsForEc;

        // Dreicer runaway rate objects
        ConnorHastie *dreicer_ConnorHastie=nullptr;
        DreicerNeuralNetwork *dreicer_nn=nullptr;

        gsl_integration_workspace *gsl_ad_w;
        gsl_root_fsolver *fsolve;
        gsl_min_fminimizer *fmin;

        OptionConstants::conductivity_mode cond_mode;
        OptionConstants::eqterm_dreicer_mode dreicer_mode;
        OptionConstants::collqty_Eceff_mode Eceff_mode;
        OptionConstants::eqterm_avalanche_mode ava_mode;
        OptionConstants::eqterm_compton_mode compton_mode;
        FVM::Interpolator1D *compton_photon_flux;

        len_t id_ncold;
        len_t id_ntot;
        len_t id_ni;
        len_t id_Tcold;
        len_t id_Eterm;
        len_t id_jtot;

        real_t *ncold;
        real_t *ntot;
        real_t *Tcold;
        real_t *Eterm;

        real_t *Ec_free=nullptr;                 // Connor-Hastie field with only free
        real_t *Ec_tot=nullptr;                  // Connor-Hastie field with free+bound
        real_t *tauEERel=nullptr;                // Relativistic electron collision time
        real_t *tauEETh=nullptr;                 // Thermal electron collision time
        
        real_t *EDreic=nullptr;                  // Dreicer field
        real_t *criticalREMomentum=nullptr;      // Critical momentum for runaway p_star 
        real_t *criticalREMomentumInvSq=nullptr; // Inverse square p_star
        real_t *pc_COMPLETESCREENING = nullptr;
        real_t *pc_NOSCREENING = nullptr;
        real_t *avalancheGrowthRate=nullptr;     // (dnRE/dt)_ava = nRE*Gamma_ava
        real_t *dreicerRunawayRate=nullptr;      // (dnRE/dt)_Dreicer = gamma_Dreicer
        real_t *tritiumRate=nullptr;             // (dnRE/dt)_Tritium = nTritium * ...
        real_t *comptonRate=nullptr;             // (dnRE/dt)_Compton = n_tot * ...
        real_t *DComptonRateDpc=nullptr;         // d/dpc((dnRE/dt)_Compton)
        real_t *effectiveCriticalField=nullptr;  // Eceff: Gamma_ava(Eceff) = 0
        real_t *electricConductivity=nullptr;
		real_t *pStar=nullptr;					 // Effective critical momentum pStar
		real_t *nusnuDatPStar=nullptr;			 // normalized nu_s*nu_D evaluated at pStar

        EffectiveCriticalField *effectiveCriticalFieldObject = nullptr; 
        
        FVM::TimeKeeper *timeKeeper;
        len_t
            timerTot,
            timerLnLambdaEE, timerLnLambdaEI,
            timerNuS, timerNuD, timerDerived,
            timerEcEff, timerPCrit, timerGrowthrates;

        bool gridRebuilt;
        bool parametersHaveChanged();

        void AllocateQuantities();
        void DeallocateQuantities();
        
        void CalculateDerivedQuantities();
        void CalculateCriticalMomentum();
        void CalculateGrowthRates(const real_t t);

        static void FindECritInterval(len_t ir, real_t *E_lower, real_t *E_upper, void *par);

        real_t BounceAverageFunc(len_t ir, std::function<real_t(real_t,real_t)> Func);

        
        static real_t pStarFunction(real_t, void *);
        static real_t pStarFunctionAlt(real_t, void *);
        real_t evaluatePStar(len_t ir, real_t E, gsl_function gsl_func, real_t *nuSHat_COMPSCREEN);
        real_t evaluateNuDHat(len_t ir, real_t p, CollisionQuantity::collqty_settings *inSettings);
        real_t evaluateNuSHat(len_t ir, real_t p, CollisionQuantity::collqty_settings *inSettings);
        
        static const len_t  conductivityLenT;
        static const len_t  conductivityLenZ;
        static const real_t conductivityBraams[];
        static const real_t conductivityTmc2[];   // list of T/mc2 
        static const real_t conductivityX[];      // where X = 1/(1+Zeff) 
        
        gsl_interp2d *gsl_cond;
        gsl_interp_accel *gsl_xacc;
        gsl_interp_accel *gsl_yacc;
        int QAG_KEY = GSL_INTEG_GAUSS31;


    protected:
    public:
        RunawayFluid(
            FVM::Grid *g, FVM::UnknownQuantityHandler *u, SlowingDownFrequency *nuS, 
            PitchScatterFrequency *nuD, CoulombLogarithm *lnLEE,
            bool extrapolateDreicer,
            CoulombLogarithm *lnLEI, IonHandler *ions, AnalyticDistributionRE *distRE,
            CollisionQuantity::collqty_settings *cqForPc, CollisionQuantity::collqty_settings *cqForEc,
            OptionConstants::conductivity_mode cond_mode,
            OptionConstants::eqterm_dreicer_mode,
            OptionConstants::collqty_Eceff_mode,
            OptionConstants::eqterm_avalanche_mode,
            OptionConstants::eqterm_compton_mode,
            FVM::Interpolator1D *compton_flux
        );
        ~RunawayFluid();

        static void FindRoot(real_t x_lower, real_t x_upper, real_t *root, gsl_function gsl_func, gsl_root_fsolver *s, real_t epsrel=1e-3, real_t epsabs=0);
        static void FindRoot_fdf(real_t &root, gsl_function_fdf gsl_func, gsl_root_fdfsolver *s, real_t epsrel=3e-3, real_t epsabs=0);
        static void FindRoot_fdf_bounded(real_t x_lower, real_t x_upper, real_t &root, gsl_function_fdf gsl_func, gsl_root_fdfsolver *s, real_t epsrel=3e-3, real_t epsabs=0);
        static void FindInterval(real_t *x_lower, real_t *x_upper, gsl_function gsl_func );

        static real_t evaluateTritiumRate(real_t gamma_c);
        static real_t evaluateComptonRate(real_t pc, real_t photonFlux, gsl_integration_workspace *gsl_ad_w);
        static real_t evaluateDComptonRateDpc(real_t pc, real_t photonFlux, gsl_integration_workspace *gsl_ad_w);
        static real_t evaluateComptonPhotonFluxSpectrum(real_t Eg, real_t photonFlux);
        static real_t evaluateComptonTotalCrossSectionAtP(real_t Eg, real_t pc);
        static real_t evaluateDSigmaComptonDpcAtP(real_t Eg, real_t pc);


        void Rebuild(const real_t);
        void GridRebuilt();
        bool GetExtrapolateDreicer()
            {return extrapolateDreicer;}
        const real_t GetEffectiveCriticalField(len_t ir) const
            {return effectiveCriticalField[ir];}
        const real_t* GetEffectiveCriticalField() const
            {return effectiveCriticalField;}
        
        const real_t GetElectricConductivity(len_t ir) const
            {return electricConductivity[ir];}
        const real_t* GetElectricConductivity() const
            {return electricConductivity;}

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

        const real_t GetDreicerRunawayRate(len_t ir) const
            { return dreicerRunawayRate[ir]; }
        const real_t *GetDreicerRunawayRate() const
            { return dreicerRunawayRate; }
        
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

		const real_t GetPStar(len_t ir) const
			{return pStar[ir];}
		const real_t *GetPStar() const
			{return pStar;}

		const real_t GetNusNuDatPStar(len_t ir) const
			{return nusnuDatPStar[ir];}
		const real_t *GetNusNuDatPStar() const
			{return nusnuDatPStar;}
        
        ConnorHastie *GetConnorHastieRunawayRate() { return this->dreicer_ConnorHastie; }
        DreicerNeuralNetwork *GetDreicerNeuralNetwork() { return this->dreicer_nn; }
        IonHandler *GetIonHandler() { return this->ions; }
        FVM::UnknownQuantityHandler *GetUnknowns() { return this->unknowns; }
        AnalyticDistributionRE *GetAnalyticDistributionRE() { return this->analyticRE; }

        CoulombLogarithm* GetLnLambda(){return lnLambdaEE;}
        SlowingDownFrequency* GetNuS(){return nuS;}
        PitchScatterFrequency* GetNuD(){return nuD;}

        real_t evaluateNeoclassicalConductivityCorrection(len_t ir, bool collisionless);
        real_t evaluateNeoclassicalConductivityCorrection(len_t ir, real_t Tcold, real_t Zeff, real_t ncold, bool collisionless);

        real_t evaluateSauterElectricConductivity(len_t ir, bool collisionless);
        real_t evaluateSauterElectricConductivity(len_t ir, real_t Tcold, real_t Zeff, real_t ncold, bool collisionless);
        real_t evaluateBraamsElectricConductivity(len_t ir);
        real_t evaluateBraamsElectricConductivity(len_t ir, real_t Tcold, real_t Zeff);

        real_t evaluatePartialContributionConductivity(len_t ir, len_t derivId, len_t n); 
        real_t evaluatePartialContributionSauterConductivity(len_t ir, len_t derivId, len_t n, bool collisionless);
        real_t evaluatePartialContributionBraamsConductivity(len_t ir, len_t derivId, len_t n);
        void evaluatePartialContributionAvalancheGrowthRate(real_t *dGamma, len_t derivId);
        void evaluatePartialContributionComptonGrowthRate(real_t *dGamma, len_t derivId);


        void PrintTimings();
        void SaveTimings(SFile*, const std::string& path="");
    };

}


#endif/*_DREAM_EQUATIONS_RUNAWAY_FLUID_HPP*/

    



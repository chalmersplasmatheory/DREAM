#ifndef _DREAM_EQUATIONS_SLOWING_DOWN_FREQUENCY_HPP
#define _DREAM_EQUATIONS_SLOWING_DOWN_FREQUENCY_HPP

#include "CollisionFrequency.hpp"
namespace DREAM {
    class SlowingDownFrequency : public CollisionFrequency {
    private:
        gsl_integration_workspace *gsl_ad_w;
        gsl_spline *bremsSpline;
        gsl_interp_accel *gsl_acc;

        static const len_t  MAX_Z = 18; // tabulated mean excitation energies up to Z = 18
        static const len_t  MAX_NE = 14; // tabulated constants for analytic formula up to Ne = 14
        static const real_t MEAN_EXCITATION_ENERGY_DATA[MAX_Z][MAX_Z];
        static const real_t MEAN_EXCITATION_ENERGY_FUNCTION_D[];
        static const real_t MEAN_EXCITATION_ENERGY_FUNCTION_S_0[];
        static const real_t HIGH_Z_EXCITATION_ENERGY_PER_Z; 
        static const real_t HYDROGEN_MEAN_EXCITATION_ENERGY;
        
        static const len_t BREMS_INTEGRAL_N = 100;
        static const real_t BREMS_INTEGRAL_X[BREMS_INTEGRAL_N];
        static const real_t BREMS_INTEGRAL_I[BREMS_INTEGRAL_N];

        void GetPartialContributionNi(real_t preFactor, real_t *hiBethe, real_t hCold, const real_t lnLee, len_t pind,len_t np1, real_t *&partQty);        

        void calculateIsotropicNonlinearOperatorMatrix() override;

        virtual real_t evaluateElectronTermAtP(len_t ir, real_t p, OptionConstants::collqty_collfreq_mode collfreq_mode) override;
        virtual real_t evaluateDDTElectronTermAtP(len_t ir, real_t p, OptionConstants::collqty_collfreq_mode collfreq_mode) override;
        virtual real_t evaluateScreenedTermAtP(len_t iz, len_t Z0, real_t p, OptionConstants::collqty_collfreq_mode collfreq_mode) override;
        virtual real_t evaluateIonTermAtP(len_t /*iz*/, len_t /*Z0*/, real_t /*p*/) override {return 0;}
        virtual real_t evaluateBremsstrahlungTermAtP(len_t iz, len_t Z0, real_t p, OptionConstants::eqterm_bremsstrahlung_mode brems_mode, OptionConstants::collqty_collfreq_type collfreq_type) override;
    protected:
        virtual real_t GetAtomicParameter(len_t iz, len_t Z0) override;        
    public:
        SlowingDownFrequency(FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ih,  
                CoulombLogarithm *lnLee,CoulombLogarithm *lnLei,
                enum OptionConstants::momentumgrid_type mgtype,  struct collqty_settings *cqset);
        ~SlowingDownFrequency();
        virtual real_t evaluatePreFactorAtP(real_t p, OptionConstants::collqty_collfreq_mode collfreq_mode) override;
        real_t GetMeanExcitationEnergy(len_t iz, len_t Z0)
            {len_t ind = ionIndex[iz][Z0]; return atomicParameter[ind];}

        real_t GetP3NuSAtZero(len_t ir);
        real_t *GetPartialP3NuSAtZero(len_t derivId); 
    };

}

#endif/*_DREAM_EQUATIONS_SLOWING_DOWN_FREQUENCY_HPP*/

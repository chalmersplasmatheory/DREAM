#ifndef _DREAM_EQUATIONS_SLOWING_DOWN_FREQUENCY_HPP
#define _DREAM_EQUATIONS_SLOWING_DOWN_FREQUENCY_HPP

#include "CollisionFrequency.hpp"

namespace DREAM {
    class SlowingDownFrequency : public CollisionFrequency {
    private:
        gsl_integration_workspace *gsl_ad_w;
        static const len_t  meanExcI_len;
        static const real_t meanExcI_data[];
        static const real_t meanExcI_Zs[];
        static const real_t meanExcI_Z0s[];
        
        void GetPartialContributionNi(real_t preFactor, real_t *hiBethe, real_t hCold, const real_t lnLee, len_t pind,len_t np1, real_t *&partQty);        

        void calculateIsotropicNonlinearOperatorMatrix();

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

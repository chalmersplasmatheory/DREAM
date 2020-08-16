#ifndef _DREAM_EQUATIONS_PITCH_SCATTER_FREQUENCY_HPP
#define _DREAM_EQUATIONS_PITCH_SCATTER_FREQUENCY_HPP

#include "CollisionFrequency.hpp"

namespace DREAM {
    class PitchScatterFrequency : public CollisionFrequency {
    private:

        static const len_t  ionSizeAj_len;
        static const real_t ionSizeAj_data[];
        static const real_t ionSizeAj_Zs[];
        static const real_t ionSizeAj_Z0s[];
        
        

        void GetPartialContributionNi(real_t preFactor, real_t *hiBethe, real_t hCold, const real_t lnLee, len_t pind,len_t np1, real_t *&partQty);        

        void calculateIsotropicNonlinearOperatorMatrix();

        virtual real_t evaluateElectronTermAtP(len_t ir, real_t p, OptionConstants::collqty_collfreq_mode collfreq_mode) override;
        virtual real_t evaluateDDTElectronTermAtP(len_t ir, real_t p, OptionConstants::collqty_collfreq_mode collfreq_mode) override;
        virtual real_t evaluateScreenedTermAtP(len_t iz, len_t Z0, real_t p,OptionConstants::collqty_collfreq_mode collfreq_mode) override;
        virtual real_t evaluateIonTermAtP(len_t iz, len_t Z0, real_t p) override;
        virtual real_t evaluateBremsstrahlungTermAtP(len_t /*iz*/, len_t /*Z0*/, real_t /*p*/, OptionConstants::eqterm_bremsstrahlung_mode /*brems_mode*/, OptionConstants::collqty_collfreq_type /*collfreq_type*/) override {return 0;}
   protected:
        virtual real_t GetAtomicParameter(len_t iz, len_t Z0) override;

    public:
        PitchScatterFrequency(FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ih,  
                CoulombLogarithm *lnLei,CoulombLogarithm *lnLee,
                enum OptionConstants::momentumgrid_type mgtype,  struct collqty_settings *cqset);
        ~PitchScatterFrequency();

        real_t GetIonEffectiveSizeAj(len_t iz, len_t Z0)
            {len_t ind = ionIndex[iz][Z0]; return atomicParameter[ind];}
        virtual real_t evaluatePreFactorAtP(real_t p, OptionConstants::collqty_collfreq_mode collfreq_mode) override;


    };

}

#endif/*_DREAM_EQUATIONS_PITCH_SCATTER_FREQUENCY_HPP*/

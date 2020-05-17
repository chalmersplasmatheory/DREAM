#ifndef _DREAM_EQUATIONS_PITCH_SCATTER_FREQUENCY_HPP
#define _DREAM_EQUATIONS_PITCH_SCATTER_FREQUENCY_HPP


#include "FVM/config.h"
#include "CollisionFrequency.hpp"
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
    class PitchScatterFrequency : public CollisionFrequency {
    private:

        static const len_t  ionSizeAj_len;
        static const real_t ionSizeAj_data[];
        static const real_t ionSizeAj_Zs[];
        static const real_t ionSizeAj_Z0s[];
        
        

        void GetPartialContributionNi(real_t preFactor, real_t *hiBethe, real_t hCold, const real_t lnLee, len_t pind,len_t np1, real_t *&partQty);        

        void calculateIsotropicNonlinearOperatorMatrix();

        virtual real_t evaluateElectronTermAtP(len_t ir, real_t p) override;
        virtual real_t evaluateScreenedTermAtP(len_t iz, len_t Z0, real_t p) override;
        virtual real_t evaluateIonTermAtP(len_t iz, len_t Z0, real_t p) override;
        virtual real_t evaluatePreFactorAtP(real_t p) override {return constPreFactor * sqrt(1+p*p)/(p*p*p);}
   protected:
        virtual real_t GetAtomicParameter(len_t iz, len_t Z0) override;

    public:
        PitchScatterFrequency(FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ih,  
                CoulombLogarithm *lnLei,CoulombLogarithm *lnLee,
                enum OptionConstants::momentumgrid_type mgtype,  struct collqty_settings *cqset);
        ~PitchScatterFrequency();

        real_t GetIonEffectiveSizeAj(len_t iz, len_t Z0)
            {len_t ind = ionIndex[iz][Z0]; return atomicParameter[ind];}


    };

}

#endif/*_DREAM_EQUATIONS_PITCH_SCATTER_FREQUENCY_HPP*/

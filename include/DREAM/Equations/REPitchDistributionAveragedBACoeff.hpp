#ifndef _DREAM_EQUATIONS_RE_PITCH_DISTRIBUTION_AVERAGED_BA_COEFF_HPP
#define _DREAM_EQUATIONS_RE_PITCH_DISTRIBUTION_AVERAGED_BA_COEFF_HPP

namespace DREAM { class REPitchDistributionAveragedBACoeff; }

#include "DREAM/DREAMException.hpp"
#include "DREAM/Equations/AnalyticDistributionRE.hpp"
#include "DREAM/Equations/PitchScatterFrequency.hpp"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

namespace DREAM {
    class REPitchDistributionAveragedBACoeff {
    private:
        FVM::RadialGrid *rGrid;
        AnalyticDistributionRE *distRE;

        /** 
         * The bounce averaged coefficient X is characterised by
         * (1) an elemental bounce average function 'BA_Func', being a simple function
         *     of xiOverXi0, BOverBmin, ROverR0 and NablaR2 parametrized by BA_Param 
         *     (or BA_Func_par, but assumed to be nullptr for now since we don't 
         *      know how to update it here - should otherwise pass instructions)
         * (2) a simple pitch-dependent prefactor function BA_PitchPrefactor(xi0)
         *     which will be integrated over the RE pitch distribution and factored
         *     into the 'REDistAverage' spline 
         * (3) a momentum and radially varying prefactor BA_MomentumPrefactor(ir,p)
         *     which multiplies the rest to produce the full coefficient
         */
        real_t(*BA_Func)(real_t,real_t,real_t,real_t,void*);
        void *BA_Func_par;
        int_t *BA_Param;
        std::function<real_t(real_t)> BA_PitchPrefactor;
        std::function<real_t(len_t,real_t)> BA_MomentumPrefactor;

        const gsl_interp_type *interp_mode; // interpolation method used by all splines
        gsl_integration_workspace *gsl_ad_w;
        static constexpr len_t N_BA_SPLINE = 50;
        gsl_spline **BA_Spline = nullptr;
        gsl_interp_accel **BA_Accel = nullptr;

        static constexpr len_t N_RE_DIST_SPLINE = 15;
        gsl_spline **REDistAverage_Spline = nullptr;
        gsl_interp_accel **REDistAverage_Accel = nullptr;

        len_t nr;

        void generateBASplines();

        void generateREDistAverageSplines();
        
        void Deallocate();

    public:
        REPitchDistributionAveragedBACoeff(
            FVM::RadialGrid*, AnalyticDistributionRE*,
            real_t(*BA_Func)(real_t,real_t,real_t,real_t,void*),
            void *BA_Func_par, int_t *BA_Param, 
            std::function<real_t(real_t)> BA_PitchPrefactor,
            std::function<real_t(len_t,real_t)> BA_MomentumPrefactor,
            const gsl_interp_type *t=gsl_interp_steffen
        );
        ~REPitchDistributionAveragedBACoeff();

        bool GridRebuilt();
        /** 
         * The RE pitch dist splines are stored and evaluated in the
         * variable X = A/(1+A) which maps the interval
         * [0,inf] in A to [0,1] in X. These methods
         * convert between A and X.
         */
        static real_t GetAFromX(real_t X)
            { return sqrt(X)/(1.0-sqrt(X)); }
        static real_t GetXFromA(real_t A)
            { real_t x = A/(1+A); return x*x; }

        /**
         * Main function of this class: evaluates 
         * the distribution-bounce averaged function, [[X]] 
         * at pitch distribution width parameter 'A'
         */
        real_t EvaluateREPitchDistAverageAtA(len_t ir, real_t p, real_t *dFdp=nullptr){
            real_t preFactor = BA_MomentumPrefactor(ir,p);
            real_t A = distRE->GetAatP(ir,p);
            real_t distIntegral = gsl_spline_eval(REDistAverage_Spline[ir], A, REDistAverage_Accel[ir]);

            if(dFdp != nullptr){
                /* TODO p derivative (if need be)
                *dFdp = preFactor * dAdp * gsl_spline_eval_deriv(REDistAverage_Spline[ir], A, REDistAverage_Accel[ir])
                    + DDPpreFactor * distIntegral;
                */
                throw DREAMException(
                    "REPitchDistributionAveragedCoeff: Evaluation of the p derivative" 
                    "of the averaged coefficients has not yet been implemented.");
            }
            return preFactor * distIntegral;
        }

        /**
         * Evaluate the splined bounce average {X} at xi0
         */
        /* Not used, and not really convenient?
        real_t EvaluateBounceAverageAtXi0(len_t ir, real_t xi0){
            return gsl_spline_eval(BA_Spline[ir], fabs(xi0), BA_Accel[ir]);
        }
        */


    };
}

#endif/*_DREAM_EQUATIONS_RE_PITCH_DISTRIBUTION_AVERAGED_BA_COEFF_HPP*/

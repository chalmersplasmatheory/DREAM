#ifndef _DREAM_MOLECULAR_RATE_INTERPOLATOR_HPP
  #define _DREAM_MOLECULAR_RATE_INTERPOLATOR_HPP

  #include <string>
  #include <gsl/gsl_interp.h>
  #include <gsl/gsl_interp2d.h>

  #include "FVM/config.h"
  #include "DREAM/MolecularRateData.hpp"

  namespace DREAM {
      class MolecularRateInterpolator {
      private:
          std::string name;

          len_t nn, nT;
          const real_t *logn;
          const real_t *logT;
          const real_t *coeff;
        
          //only need one pointer
          gsl_interp2d *interp;
          gsl_interp_accel *nacc;
          gsl_interp_accel *Tacc;

      public:
          MolecularRateInterpolator(const molecular_rate*);
          virtual ~MolecularRateInterpolator();

        real_t Eval(real_t n, real_t T) const;
        real_t Eval_deriv_n(const real_t n, const real_t T) const;
        real_t Eval_deriv_T(const real_t n, const real_t T) const;
      };
  }

  #endif

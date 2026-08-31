#include <cmath>

#include "DREAM/MolecularRateInterpolator.hpp"
#include <limits>

  using namespace DREAM;

  MolecularRateInterpolator::MolecularRateInterpolator(const molecular_rate *rate)
      : name(rate->name), nn(rate->nn), nT(rate->nT),
        logn(rate->logn), logT(rate->logT), coeff(rate->coeff) {
    

      this->interp = gsl_interp2d_alloc(gsl_interp2d_bilinear, nn, nT);
      gsl_interp2d_init(this->interp, logn, logT, coeff, nn, nT);

      this->nacc = gsl_interp_accel_alloc();
      this->Tacc = gsl_interp_accel_alloc();
  }

  MolecularRateInterpolator::~MolecularRateInterpolator() {
      gsl_interp2d_free(this->interp);
      gsl_interp_accel_free(this->nacc);
      gsl_interp_accel_free(this->Tacc);
  }

  real_t MolecularRateInterpolator::Eval(real_t n, real_t T) const {
      if (n <= 0 || T <= 0)
          return 0;

      const real_t ln = std::log10(n);
      const real_t lT = std::log10(T);

      real_t logValue = gsl_interp2d_eval_extrap(
          this->interp,
          this->logn,
          this->logT,
          this->coeff,
          ln,
          lT,
          this->nacc,
          this->Tacc
      );

      real_t value = std::pow(10.0, logValue);
      if (value < 0)
          return 0;

      return value;
  }

  real_t MolecularRateInterpolator::Eval_deriv_n(const real_t n, const real_t T) const {
    real_t eps = std::numeric_limits<real_t>::epsilon();
    real_t dn = sqrt(eps)*n + eps;
    return ( Eval(n+dn,T) - Eval(n,T) ) / dn;
}


real_t MolecularRateInterpolator::Eval_deriv_T(const real_t n, const real_t T) const {
      real_t eps = std::numeric_limits<real_t>::epsilon();
      real_t dT = sqrt(eps)*T + eps;
      return (Eval(n,T+dT) - Eval(n,T)) / dT;
  }


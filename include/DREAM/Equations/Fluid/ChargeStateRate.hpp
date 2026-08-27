//Goal is to remove this file later when mearching is easier
#ifndef _DREAM_CHARGE_STATE_RATE_HPP
#define _DREAM_CHARGE_STATE_RATE_HPP

#include <cmath>
#include <limits>
#include <string>
#include "FVM/config.h"
#include "DREAM/ADASRateInterpolator.hpp"
#include "DREAM/MolecularRateInterpolator.hpp"


namespace DREAM {

class ChargeStateRate {
    public:
        virtual ~ChargeStateRate() {}

        virtual real_t Eval(const len_t Z0,  real_t n,  real_t T) const = 0;

        virtual real_t Eval_deriv_n(const len_t Z0,  real_t n,  real_t T) const = 0;


        virtual real_t Eval_deriv_T(const len_t Z0,  real_t n,  real_t T) const = 0;

        virtual const std::string& GetName() const = 0;
    };

class ZeroChargeStateRate : public ChargeStateRate {
    private:
        std::string name;
            
    public:
        ZeroChargeStateRate(const std::string& name) : name(name) {}
        virtual real_t Eval(const len_t,  real_t,  real_t) const override {
                return 0;}
        virtual real_t Eval_deriv_n(const len_t,  real_t,  real_t) const override {
                return 0;}
        virtual real_t Eval_deriv_T(const len_t,  real_t,  real_t) const override {
                return 0;}
        virtual const std::string& GetName() const override {
                return name;}
};

class MolecularTableChargeStateRate : public ChargeStateRate {
  private:
      std::string name;
      len_t activeZ0;
      MolecularRateInterpolator *rate;

  public:
      MolecularTableChargeStateRate(
          const std::string& name,
          const len_t activeZ0,
          MolecularRateInterpolator *rate
      ) : name(name), activeZ0(activeZ0), rate(rate) {}

      virtual ~MolecularTableChargeStateRate() override {
          delete rate;
      }

      virtual real_t Eval(const len_t Z0, real_t n, real_t T) const override {
          if (Z0 != activeZ0)
              return 0;
          return rate->Eval(n, T);
      }

      virtual real_t Eval_deriv_n(const len_t Z0, real_t n, real_t T) const override {
          if (Z0 != activeZ0)
              return 0;
          return rate->Eval_deriv_n(n, T);
      }

      virtual real_t Eval_deriv_T(const len_t Z0, real_t n, real_t T) const override {
          if (Z0 != activeZ0)
              return 0;
          return rate->Eval_deriv_T(n, T);
      }

      virtual const std::string& GetName() const override {
          return name;
      }
  };




    class ADASChargeStateRate : public ChargeStateRate {
    private:
        std::string name;
        ADASRateInterpolator *adasrate;
    public:
        ADASChargeStateRate(const std::string& name, ADASRateInterpolator *adasrate) : name(name), adasrate(adasrate) {}
        virtual const std::string& GetName() const override { return name; }
        virtual real_t Eval(const len_t Z0,  real_t n,  real_t T) const override {
            return adasrate->Eval(Z0, n, T);
        }
        virtual real_t Eval_deriv_n(const len_t Z0,  real_t n,  real_t T) const override {
            return adasrate->Eval_deriv_n(Z0, n, T);
        }
        virtual real_t Eval_deriv_T(const len_t Z0,  real_t n,  real_t T) const override {
            return adasrate->Eval_deriv_T(Z0, n, T);
        }
    };
}
  #endif

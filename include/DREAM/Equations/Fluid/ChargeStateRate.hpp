
#ifndef _DREAM_CHARGE_STATE_RATE_HPP
#define _DREAM_CHARGE_STATE_RATE_HPP

#include <cmath>
#include <limits>
#include <string>
#include "FVM/config.h"
#include "DREAM/ADASRateInterpolator.hpp"

namespace DREAM {

class ChargeStateRate {
    public:
        virtual ~ChargeStateRate() {}

        virtual real_t Eval(const len_t Z0,  real_t n,  real_t T) const = 0;

        virtual real_t Eval_deriv_n(const len_t Z0,  real_t n,  real_t T) const {
            real_t eps = std::numeric_limits<real_t>::epsilon();
            real_t dn = sqrt(eps)*n + eps;
            return (Eval(Z0, n+dn, T) - Eval(Z0, n, T)) / dn;
        }

        virtual real_t Eval_deriv_T(const len_t Z0,  real_t n,  real_t T) const {
            real_t eps = std::numeric_limits<real_t>::epsilon();
            real_t dT = sqrt(eps)*T + eps;
            return (Eval(Z0, n, T+dT) - Eval(Z0, n, T)) / dT;
        }

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

class TabledChargeStateRate : public ChargeStateRate {
    private:
        std::string name;
        len_t activeZ0;
        len_t nT, nn;
        const real_t *coeff;
    public:
        TabledChargeStateRate(
            const std::string& name, const len_t activeZ0,
            const len_t nT, const len_t nn, const real_t *coeff
        ) : name(name), activeZ0(activeZ0), nT(nT), nn(nn), coeff(coeff) {}
        
        virtual const std::string& GetName() const override { return name; }
        virtual real_t Eval(const len_t Z0,  real_t n,  real_t T) const override {
                if (Z0 != activeZ0)
                    return 0;
                if(n>1e22)
                    n=1e22;
                else if(n<1e14)
                    n=1e14;

                const real_t n_cm3 = n * 1e-6 * 1e-8; //legacy conversoin
                const real_t logT = std::log(T);
                const real_t logn = std::log(n_cm3);

                real_t sum = 0;
                real_t pT = 1;

                for (len_t iT = 0; iT < nT; iT++) {
                    real_t pn = 1;
                    for (len_t in = 0; in < nn; in++) {
                        sum += coeff[iT*nn + in] * pT * pn;
                        pn *= logn;
                    }
                    pT *= logT;
                }
                return std::exp(sum) * 1e-6; // Factor 1e6 converts from cm^3 to m^3

        }
          virtual real_t Eval_deriv_n(const len_t Z0,  real_t n,  real_t T) const override {
            if (Z0 != activeZ0)
                return 0;

            if (n <= 0 || T <= 0)
                return 0;

            const real_t nScaled = n * 1e-14;
            const real_t lnT = std::log(T);
            const real_t lnn = std::log(nScaled);

            real_t sum = 0;
            real_t dsum_dn = 0;

            real_t pT = 1;
            for (len_t iT = 0; iT < nT; iT++) {
                real_t pn = 1;
                real_t pn_minus_1 = 1;

                for (len_t in = 0; in < nn; in++) {
                    const real_t c = coeff[iT*nn + in];

                    sum += c * pT * pn;

                    if (in > 0)
                        dsum_dn += c * pT * in * pn_minus_1 / n;

                    pn_minus_1 = pn;
                    pn *= lnn;
                }

                pT *= lnT;
            }

            return std::exp(sum) * dsum_dn * 1e-6;
        }


        virtual real_t Eval_deriv_T(const len_t Z0,  real_t n,  real_t T) const override {
            if (Z0 != activeZ0)
                return 0;

            if (n <= 0 || T <= 0)
                return 0;

            const real_t nScaled = n * 1e-14;
            const real_t lnT = std::log(T);
            const real_t lnn = std::log(nScaled);

            real_t sum = 0;
            real_t dsum_dT = 0;

            real_t pT = 1;
            real_t pT_minus_1 = 1;

            for (len_t iT = 0; iT < nT; iT++) {
                real_t pn = 1;

                for (len_t in = 0; in < nn; in++) {
                    const real_t c = coeff[iT*nn + in];

                    sum += c * pT * pn;

                    if (iT > 0)
                        dsum_dT += c * iT * pT_minus_1 * pn / T;

                    pn *= lnn;
                }

                pT_minus_1 = pT;
                pT *= lnT;
            }

            return std::exp(sum) * dsum_dT * 1e-6;
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

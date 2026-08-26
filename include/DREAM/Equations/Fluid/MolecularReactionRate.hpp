#ifndef _DREAM_MOLECULAR_REACTION_RATE_HPP
#define _DREAM_MOLECULAR_REACTION_RATE_HPP

  #include <string>
  #include "FVM/config.h"

  namespace DREAM {

  class MolecularReactionRate {
  public:
      virtual ~MolecularReactionRate() {}
      virtual real_t Eval(const len_t Z0, real_t n, real_t T) const = 0;
      virtual const std::string& GetName() const = 0;
  };

  class ZeroMolecularReactionRate : public MolecularReactionRate {
  private:
      std::string name;
  public:
      ZeroMolecularReactionRate(const std::string& name) : name(name) {}
      real_t Eval(const len_t Z0, real_t n, real_t T) const override { return 0; }
      const std::string& GetName() const override { return name; }
  };

  class ConstantMolecularReactionRate : public MolecularReactionRate {
  private:
      std::string name;
      real_t value;
  public:
      ConstantMolecularReactionRate(const std::string& name, real_t value)
          : name(name), value(value) {}
      real_t Eval(const len_t Z0, real_t n, real_t T) const override { return value; }
      const std::string& GetName() const override { return name; }
  };

  class EnergyInterpolatedMolecularReactionRate : public MolecularReactionRate {
  private:
      std::string name;
      const real_t *energyPoints;
      const real_t *crossSectionPoints;
    public:
      EnergyInterpolatedMolecularReactionRate(
          const std::string& name, const real_t *energyPoints, const real_t *crossSectionPoints
      ) : name(name), energyPoints(energyPoints), crossSectionPoints(crossSectionPoints) {}
      real_t Eval(const len_t Z0, real_t n, real_t T) const override {

            //just read table to get the rate for given density and temp 
          return 0; // Placeholder
      }
      const std::string& GetName() const override { return name; }

  };

  class AMJUELTable : public MolecularReactionRate {
    private:
        std::string name;
        len_t activeZ0;
        len_t nT, nn;
        const real_t *coeff;
    public:
        AMJUELTable(
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
        
    };
  
  } // namespace DREAM

  #endif

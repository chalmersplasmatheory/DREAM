#ifndef _DREAM_EQUATIONS_FLUID_RATEHANDLER_HPP
#define _DREAM_EQUATIONS_FLUID_RATEHANDLER_HPP


#include <string>
#include <unordered_map>
#include "FVM/config.h"
#include "DREAM/IonHandler.hpp"
#include "DREAM/ADAS.hpp"
 #include "DREAM/Equations/Fluid/ChargeStateRate.hpp"

namespace DREAM {

    struct ChargeStateRateSet {
      ChargeStateRate *acd;
      ChargeStateRate *scd;
  };

    class RateHandler {
    private:
        IonHandler *ions;
        ADAS *adas;

        //connecting the name of the molecule to the corresponding charge state rates
        std::unordered_map<std::string, ChargeStateRateSet> chargeStateRates;
        

        void AddMolecularChargeStateRates();
        //void AddMolecularRates()
        void AddAtomicChargeStateRates();

    public:
        RateHandler(IonHandler *ions, ADAS *adas);
        ~RateHandler();


        ChargeStateRate *GetACD(const std::string& name) const;
        ChargeStateRate *GetSCD(const std::string& name) const;
  };

  }

  #endif


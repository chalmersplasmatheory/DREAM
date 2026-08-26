#ifndef _DREAM_RATE_HANDLER_HPP
#define _DREAM_RATE_HANDLER_HPP


#include <string>
#include <unordered_map>
#include "FVM/config.h"
#include "DREAM/IonHandler.hpp"
#include "DREAM/ADAS.hpp"
#include "DREAM/Equations/Fluid/ChargeStateRate.hpp"
#include "DREAM/Equations/Fluid/RateData.hpp"
#include "DREAM/Equations/Fluid/MolecularReactionRate.hpp"
#include <vector>

namespace DREAM {

    //Struct to hold the charge state rates for a given species
    struct ChargeStateRateSet {
      ChargeStateRate *acd;
      ChargeStateRate *scd;

    };
  
    //Struct to hold the molecular reaction rates for a given pair of species
    struct MolecularRatePair {
    std::string reactant1Name;
    len_t reactant1Z0;

    std::string reactant2Name;
    len_t reactant2Z0;

    std::string product1Name;
    len_t product1Z0;
    std::string product2Name;
    len_t product2Z0;

    MolecularReactionRate *chargeExchange;
    MolecularReactionRate *dissociation;
    MolecularReactionRate *dissociativeRecombination;
  };


    class RateHandler {
    private:
        IonHandler *ions;
        ADAS *adas;

        std::unordered_map<std::string, ChargeStateRateSet> chargeStateRates;
        std::vector<MolecularRatePair> molecularRatePairs;

        MolecularReactionRate *CreateMolecularReactionRate(const std::string& name, const MolecularRateDefinition& def);
        
        void AddMolecularChargeStateRates();
        void AddAtomicChargeStateRates();
        void AddMolecularReactionRates();

    public:
        RateHandler(IonHandler *ions, ADAS *adas);
        ~RateHandler();


        const std::vector<MolecularRatePair>& GetMolecularRatePairs() const {
                return molecularRatePairs;
            }

        ChargeStateRate *GetACD(const std::string& name) const;
        ChargeStateRate *GetSCD(const std::string& name) const;

        
  };

  }

  #endif /*_DREAM_RATE_HANDLER_HPP*/


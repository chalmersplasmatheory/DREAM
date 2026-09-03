#ifndef _DREAM_RATE_HANDLER_HPP
#define _DREAM_RATE_HANDLER_HPP


#include <string>
#include <unordered_map>
#include "FVM/config.h"
#include "DREAM/IonHandler.hpp"
#include "DREAM/ADAS.hpp"
#include "DREAM/Equations/Fluid/ChargeStateRate.hpp"
#include "DREAM/Equations/Fluid/RateData.hpp"
#include "DREAM/MolecularRateInterpolator.hpp"
#include <vector>

namespace DREAM {

    //Struct to hold the charge state rates for a given species
    struct ChargeStateRateSet {
      ChargeStateRate *acd;
      ChargeStateRate *scd;

    };

   // enum class MolecularReactionType {
   //       CHARGE_EXCHANGE,
   //       IONIZATION,
   //       RECOMBINATION,
    //      DISSOCIATION,
     //     DISSOCIATIVE_IONIZATION,
      //    DISSOCIATIVE_RECOMBINATION,
      //    ION_MOLECULE_CONVERSION,
     //     REACTIVE_CHARGE_TRANSFER
     // };
  
    //Struct to hold the molecular reaction rates for a given pair of species
    struct MolecularReaction {
      const char *rateName;

          MolecularReactionProcess process;

          len_t nReactants;
          const MolecularReactionSpecies *reactants;

          len_t nProducts;
          const MolecularReactionSpecies *products;

        MolecularRateInterpolator *rate;
  };


    class RateHandler {
    private:
        IonHandler *ions;
        ADAS *adas;

        std::unordered_map<std::string, ChargeStateRateSet> chargeStateRates;
        std::vector<MolecularReaction> molecularReactions;

        MolecularRateInterpolator *GetMolecularRateByName(const char *rateName) const;
        
        void AddMolecularChargeStateRates();
        void AddAtomicChargeStateRates();
        void AddMolecularReactionRates();


    public:
        RateHandler(IonHandler *ions, ADAS *adas);
        ~RateHandler();



        //const std::vector<MolecularReaction>& GetMolecularRatePairs() const {
          //      return molecularReactions;
           // }

        ChargeStateRate *GetACD(const std::string& name) const;
        ChargeStateRate *GetSCD(const std::string& name) const;

        
  };

  }

  #endif /*_DREAM_RATE_HANDLER_HPP*/


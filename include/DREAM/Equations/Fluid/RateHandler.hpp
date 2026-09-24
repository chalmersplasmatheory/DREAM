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
#include "DREAM/Settings/Settings.hpp"
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
        MolecularInput temperatureInput;
        MolecularInput densityInput;

  };


    class RateHandler {
    private:
        FVM::Grid *grid;
        IonHandler *ions;
        ADAS *adas;
        bool separateNeutrals;
        len_t id_Wn;

        len_t id_ions, id_n_cold, id_Ni, id_T_cold, id_Wi;

        std::unordered_map<std::string, ChargeStateRateSet> chargeStateRates;
        std::vector<MolecularReaction> molecularReactions;

        MolecularRateInterpolator *GetMolecularRateByName(const char *rateName) const;
        
        void AddMolecularChargeStateRates();
        void AddAtomicChargeStateRates();
        void AddMolecularReactionRates( const std::vector<std::string>& enabledReactionNames);

        FVM::UnknownQuantityHandler *unknowns;

    public:
        RateHandler(FVM::Grid *grid, IonHandler *ions, ADAS *adas, FVM::UnknownQuantityHandler *unknowns, bool reactionsEnabled, const std::vector<std::string>& enabledReactionNames);
        ~RateHandler();



      const std::vector<MolecularReaction>& GetMolecularReactions() const {
      return molecularReactions;
      }

        ChargeStateRate *GetACD(const std::string& name) const;
        ChargeStateRate *GetSCD(const std::string& name) const;

      real_t ResolveDensity(
            const MolecularInput& input,
            const len_t ir
        ) ;

      real_t ResolveTemperature(
            const MolecularInput& input,
            const len_t ir
        ) ;
      real_t ResolveSpeciesTemperature(
            len_t species, int_t charge, len_t ir
        ) ;

        
  };

  }

  #endif /*_DREAM_RATE_HANDLER_HPP*/


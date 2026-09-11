#ifndef _DREAM_RATE_DATA_HPP
#define _DREAM_RATE_DATA_HPP

  #include "FVM/config.h"

  namespace DREAM {

      enum class MolecularReactionProcess {
          CHARGE_EXCHANGE,
          CHARGE_EXCHANGE_RESONANT,
          IONIZATION,
          IONIZATION_RUNAWAY,
          RECOMBINATION,
          DISSOCIATION,
          DISSOCIATION_MULTIPLE,
          DISSOCIATION_RUNAWAY,
      };



      struct MolecularReactionSpecies {
          const char *name;
          int_t Z0;
          len_t coefficient;
      };

      struct MolecularReactionDefinition {
          const char *rateName;

          MolecularReactionProcess process;

          len_t nReactants;
          const MolecularReactionSpecies *reactants;

          len_t nProducts;
          const MolecularReactionSpecies *products;
      };

      extern const len_t molecularReactionDefinitionCount;
      extern const MolecularReactionDefinition
      molecularReactionDefinitions[];

  }

  #endif /* _DREAM_RATE_DATA_HPP */

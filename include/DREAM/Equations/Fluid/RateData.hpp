#ifndef _DREAM_RATE_DATA_HPP
#define _DREAM_RATE_DATA_HPP

  #include "FVM/config.h"

  namespace DREAM {

      enum class MolecularReactionProcess {
          CHARGE_EXCHANGE,
          IONIZATION,
          RECOMBINATION,
          DISSOCIATION,
          DISSOCIATIVE_IONIZATION,
          DISSOCIATIVE_RECOMBINATION,
          ION_MOLECULE_CONVERSION,
          REACTIVE_CHARGE_TRANSFER
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

#ifndef _DREAM_RATE_DATA_HPP
#define _DREAM_RATE_DATA_HPP

#include "FVM/config.h"

namespace DREAM {

    //Type of molecular reaction process.
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
    
    //Type of molecular input for a reaction rate.
    enum class MolecularInputKind {
        NONE,
        ELECTRON,
        SPECIES,
        RELATIVE,
        NEUTRAL
    };

    //Structure describing a molecular input for a reaction rate.
    //It happens that the input depends on both the reactants
    struct MolecularInput {
        MolecularInputKind kind;

        const char *name1;
        int_t charge1;

        const char *name2;
        int_t charge2;
    };

    //Structure desribing a species involved in a molecular reaction.
    // Name, charge state and number nr of particles for the species. 
    struct MolecularReactionSpecies {
          const char *name;
          int_t Z0;
          len_t coefficient;
      };

    //Structure describing a molecular reaction.
    struct MolecularReactionDefinition {
        const char *rateName;

        MolecularReactionProcess process;

        len_t nReactants;
        const MolecularReactionSpecies *reactants;

        len_t nProducts;
        const MolecularReactionSpecies *products;

        MolecularInput temperatureInput;
        MolecularInput densityInput;

    };

    // Number of molecular reactions defined in the molecularReactionDefinitions array.
    extern const len_t molecularReactionDefinitionCount;
    // Array of molecular reaction definitions stored in RateData.cpp.
    extern const MolecularReactionDefinition molecularReactionDefinitions[];

  }

  #endif /* _DREAM_RATE_DATA_HPP */

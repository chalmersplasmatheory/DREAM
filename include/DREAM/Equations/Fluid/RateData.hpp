#ifndef _DREAM_RATE_DATA_HPP
#define _DREAM_RATE_DATA_HPP

#include "FVM/config.h"

namespace DREAM {



    //How to identify a molecule (should probably add iz)
    struct RateSpeciesState {
        const char *name;
        len_t Z0;
    };

    //Definition of a molecular reaction rate
    struct MolecularRateDefinition {
        const char *rateName;
    };

    //Definition of a pair of molecules and their reaction rates
    struct MolecularRatePairDefinition {
        RateSpeciesState reactant1;
        RateSpeciesState reactant2;

        RateSpeciesState product1;
        RateSpeciesState product2;

        MolecularRateDefinition chargeExchange;
        MolecularRateDefinition dissociation;
        MolecularRateDefinition dissociativeRecombination;
    };

    
    extern const len_t molecularRatePairDefinitionCount; // number of molecular rate pair definitions
    extern const MolecularRatePairDefinition molecularRatePairDefinitions[]; // array of molecular rate pair definitions
  }

  #endif/*_DREAM_RATE_DATA_HPP*/

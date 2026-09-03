#include "DREAM/Equations/Fluid/RateHandler.hpp"
#include "FVM/FVMException.hpp"
#include "DREAM/Equations/Fluid/RateData.hpp"
#include "DREAM/MoleculeHandler.hpp"
#include "DREAM/MolecularRateData.hpp"
#include "DREAM/MolecularRateInterpolator.hpp"



  using namespace DREAM;

/**
 * Constructor.
  */
    RateHandler::RateHandler(IonHandler *ions, ADAS *adas)
      : ions(ions), adas(adas) {
    
    //Set rates for IonRateEquation
    AddMolecularChargeStateRates();
    AddAtomicChargeStateRates();
      
    //Add and set rates for MolecularRateEquation (to be implemented)  
    AddMolecularReactionRates();

    //Add and set rates to runaway ionizationfluid equation (to be implemented)
  }

/**
 * Destructor.
 */
RateHandler::~RateHandler() {
      for (auto it = chargeStateRates.begin(); it != chargeStateRates.end(); it++) {
          delete it->second.acd;
          delete it->second.scd;
      }

      for (auto it = molecularReactions.begin(); it != molecularReactions.end(); it++) {
          delete it->rate;
      }
  }


  /**
   * Add charge-state rates for atomic species.
   */
  void RateHandler::AddAtomicChargeStateRates() {
    for (len_t iIon = 0; iIon <ions->GetNZ(); iIon++){
        const std::string& name = ions->GetName(iIon);
        //check if not a molecule
        if (MoleculeHandler().IsMolecule(name))
            continue;
        const len_t Z = ions->GetZ(iIon);

        if (!adas->HasElement(Z))
              continue;

        ChargeStateRateSet rates;
        rates.acd = new ADASChargeStateRate(name + "_ACD", adas->GetACD(Z));
        rates.scd = new ADASChargeStateRate(name + "_SCD", adas->GetSCD(Z));
        chargeStateRates[name] = rates;
        printf(
            "RateHandler: Added charge-state rates for atomic species '%s' (Z=%d).\n",
            name.c_str(), Z
        );
    } 
    printf("RateHandler: Added charge-state rates for %d species.\n", chargeStateRates.size());
  }

/**
 * Add charge-state rates for molecular species. 
 * Currently separate from the adas and only for D2
 */
void RateHandler::AddMolecularChargeStateRates() {
      MoleculeHandler molecules;

      for (len_t iIon = 0; iIon < ions->GetNZ(); iIon++) {
          const std::string& name = ions->GetName(iIon);
          if (!molecules.IsMolecule(name))
              continue;
          ChargeStateRateSet rates;
          if (name == "D2") {
              rates.acd = new ZeroChargeStateRate("D2_ACD_zero");
              rates.scd = new MolecularTableChargeStateRate(
                "D2_SCD_AMJUEL_2.2.9",
                0,
                GetMolecularRateByName("D_2_ionization")
            );
          } 
          chargeStateRates[name] = rates;
      }
      printf("RateHandler: Added charge-state rates for %d molecular species.\n", chargeStateRates.size());
  }


/**
 * 
 * Get the ACD (recombination) charge-state rate for a given species.
 */
  ChargeStateRate *RateHandler::GetACD(const std::string& name) const {
      auto it = chargeStateRates.find(name);
      if (it == chargeStateRates.end())
          throw FVM::FVMException(
              "RateHandler: No ACD-like charge-state rate available for '%s'.",
              name.c_str()
          );

      return it->second.acd;
  }


/**
 * Get the SCD (ionization) charge-state rate for a given species.
 */
  ChargeStateRate *RateHandler::GetSCD(const std::string& name) const {
      auto it = chargeStateRates.find(name);
      if (it == chargeStateRates.end())
          throw FVM::FVMException(
              "RateHandler: No SCD-like charge-state rate available for '%s'.",
              name.c_str()
          );

      return it->second.scd;
  }

/*
 * Add molecular reaction rates for all defined molecular rate pairs.
 */
void RateHandler::AddMolecularReactionRates() {
    //Loop though all the molecular reactions in Ratedata.cpp
    for (len_t i = 0; i < molecularReactionDefinitionCount; i++) { 
        
        ///Here we should add a check if the reactants and products are defined in the ionhandler and if not throw an exception
        const MolecularReactionDefinition& def = molecularReactionDefinitions[i];

        MolecularReaction reaction; //create a new molecularreaction
        reaction.rateName = def.rateName;
        reaction.process = def.process;
        reaction.nReactants = def.nReactants;
        reaction.reactants = def.reactants;
        reaction.nProducts = def.nProducts;
        reaction.products = def.products;
    
  
    printf(
            "RateHandler: Added %zu molecular reaction.\n",
            molecularReactions.size()
        );

        //Set the correct rate. No problem since every reaction only has 1 rate
        reaction.rate = GetMolecularRateByName(def.rateName);
        molecularReactions.push_back(reaction);
  }
  printf(
          "RateHandler: Added %zu molecular reaction pairs.\n",
          molecularReactions.size()
      );
}
MolecularRateInterpolator *RateHandler::GetMolecularRateByName(
      const char *rateName
  ) const {
      for (len_t i = 0; i < molecular_rate_n; i++) {
            //Pick up the correct rate in MolecularRateData.cpp by the correct name
          if (std::string(molecular_rate_table[i].name) == rateName) 
              return new MolecularRateInterpolator(&molecular_rate_table[i]);
      }

      throw FVM::FVMException(
          "RateHandler: No molecular rate table named '%s'.",
          rateName
      );
  }

#include "DREAM/Equations/Fluid/RateHandler.hpp"
#include "FVM/FVMException.hpp"
#include "DREAM/Equations/Fluid/RateData.hpp"
#include "DREAM/MoleculeHandler.hpp"
#include "DREAM/Equations/Fluid/AMJUELData.hpp"


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

      for (auto it = molecularRatePairs.begin(); it != molecularRatePairs.end(); it++) {
          delete it->chargeExchange;
          delete it->dissociation;
          delete it->dissociativeRecombination;
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
 * Currently separate from the adas. and only for D2
 */
void RateHandler::AddMolecularChargeStateRates() {
      MoleculeHandler molecules;

      for (len_t iIon = 0; iIon < ions->GetNZ(); iIon++) {
          const std::string& name = ions->GetName(iIon);
          if (!molecules.IsMolecule(name))
              continue;
          ChargeStateRateSet rates;
          if (name == "D2") {
              printf(
                  "RateHandler: Adding molecular charge-state rates for '%s'.\n",
                  name.c_str()
              );

              rates.acd = new ZeroChargeStateRate("D2_ACD_zero");
              rates.scd = new TabledChargeStateRate(
                  "D2_SCD_AMJUEL_2.2.9",
                  0,
                  D2_ioniz_229_nT,
                  D2_ioniz_229_nn,
                  D2_ioniz_229_coeff
              );
          } else {
              printf(
                  "WARNING: RateHandler: Molecule '%s' has no implemented charge-state rates. "
                  "Using zero ACD/SCD rates.\n",
                  name.c_str()
              );

              rates.acd = new ZeroChargeStateRate(name + "_ACD_zero");
              rates.scd = new ZeroChargeStateRate(name + "_SCD_zero");
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
      for (len_t i = 0; i < molecularRatePairDefinitionCount; i++) { //loop over all molecular rate pair definitions
          const MolecularRatePairDefinition& def = molecularRatePairDefinitions[i];

          MolecularRatePair pair; //create a new molecular rate pair object in ratehandler
          pair.reactant1Name = def.reactant1.name;
          pair.reactant1Z0   = def.reactant1.Z0;
          pair.reactant2Name = def.reactant2.name;
          pair.reactant2Z0   = def.reactant2.Z0;
        pair.product1Name  = def.product1.name;
            pair.product1Z0    = def.product1.Z0;
            pair.product2Name  = def.product2.name;
            pair.product2Z0    = def.product2.Z0;

          const std::string prefix = pair.reactant1Name + "_" + pair.reactant2Name;
            // create the molecular reaction rates for the pair and add them to the ratehandler
          pair.chargeExchange = CreateMolecularReactionRate(
              prefix + "_chargeExchange",
              def.chargeExchange
          );

          pair.dissociation = CreateMolecularReactionRate(
              prefix + "_dissociation",
              def.dissociation
          );

          pair.dissociativeRecombination = CreateMolecularReactionRate(
              prefix + "_dissociativeRecombination",
              def.dissociativeRecombination
          );

          molecularRatePairs.push_back(pair);
      }

      printf(
          "RateHandler: Added " LEN_T_PRINTF_FMT " molecular reaction pairs.\n",
          molecularRatePairs.size()
      );
  }
/* Create a molecular reaction rate object based on the given definition. */
MolecularReactionRate *RateHandler::CreateMolecularReactionRate(
      const std::string& name,
      const MolecularRateDefinition& def
  ) {
      switch (def.model) {
          case MolecularRateModel::ZERO:
              return new ZeroMolecularReactionRate(name + "_zero");

          case MolecularRateModel::CONSTANT:
              return new ConstantMolecularReactionRate(
                  name + "_constant",
                  def.constantValue
              );

          case MolecularRateModel::AMJUEL_POLYNOMIAL:
              return new AMJUELTable(
                  name + "_AMJUEL",
                  0, //activeZ0
                  def.nT,
                  def.nn,
                  def.coeff
              );
            case MolecularRateModel::INTERPOLATE_ENERGY:
                return new EnergyInterpolatedMolecularReactionRate(
                    name + "_energyInterp",
                    def.energyPoints,
                    def.crossSectionPoints
                );
      }

      throw FVM::FVMException("RateHandler: Unknown molecular rate model.");
  }
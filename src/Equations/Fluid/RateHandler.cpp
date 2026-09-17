#include "DREAM/Equations/Fluid/RateHandler.hpp"
#include "FVM/FVMException.hpp"
#include "DREAM/Equations/Fluid/RateData.hpp"
#include "DREAM/MoleculeHandler.hpp"
#include "DREAM/MolecularRateData.hpp"
#include "DREAM/MolecularRateInterpolator.hpp"
#include <algorithm>
#include "DREAM/Constants.hpp"



  using namespace DREAM;

/**
 * Constructor.
  */
    RateHandler::RateHandler(FVM::Grid *grid, IonHandler *ions, ADAS *adas, FVM::UnknownQuantityHandler *unknowns, bool reactionsEnabled, const std::vector<std::string>& enabledReactionNames)
      :  grid(grid), ions(ions), adas(adas), unknowns(unknowns) {

    this->unknowns  = unknowns;
    this->id_ions   = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
	this->id_n_cold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
	this->id_Ni  = unknowns->GetUnknownID(OptionConstants::UQTY_N_TOT);
	this->id_T_cold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->id_Wi     = unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER);
    
    //Set rates for IonRateEquation
    AddMolecularChargeStateRates();
    AddAtomicChargeStateRates();

      
    if (reactionsEnabled)
          AddMolecularReactionRates(enabledReactionNames);

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
          } else {
        rates.acd = new ZeroChargeStateRate("zero_ACD"); //this needs to be generalized
        rates.scd = new ZeroChargeStateRate("zero_SCD"); //this needs to be generalized
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
void RateHandler::AddMolecularReactionRates( const std::vector<std::string>& enabledReactionNames) {
      MoleculeHandler molecules;

      for (len_t i = 0; i < molecularReactionDefinitionCount; i++) {
          const MolecularReactionDefinition& def = molecularReactionDefinitions[i];

          // Skip definitions that Python did not select.
          if (std::find(
                  enabledReactionNames.begin(),
                  enabledReactionNames.end(),
                  std::string(def.rateName)
              ) == enabledReactionNames.end())
              continue;

          bool allReactantsExist = true;

          for (len_t j = 0; j < def.nReactants; j++) {
              const MolecularReactionSpecies& r = def.reactants[j];

              //we dont care if it has an electron
              if (std::string(r.name) == "e")
                  continue;

                //if not all reactants exist (remebers they are added as ions to the system before)
              if (!ions->HasIon(r.name)) {
                  allReactantsExist = false;
                  break;
              }
            
              //check if the charge state is valid for the reactant
              const len_t iIon = ions->GetIonIndex(r.name);
              if (r.Z0 < 0 || (len_t)r.Z0 > ions->GetZ(iIon)) {
                  allReactantsExist = false;
                  break;
              }
          }

          if (!allReactantsExist)
              continue;

          for (len_t j = 0; j < def.nProducts; j++) {
              const MolecularReactionSpecies& p = def.products[j];

              if (std::string(p.name) == "e")
                  continue;

              if (!ions->HasIon(p.name))
                  throw FVM::FVMException(
                      "RateHandler: Molecular reaction '%s' is enabled, but product species '%s' is not defined.",
                      def.rateName, p.name
                  );

              const len_t iIon = ions->GetIonIndex(p.name);
              if (p.Z0 < 0 || (len_t)p.Z0 > ions->GetZ(iIon))
                  throw FVM::FVMException(
                      "RateHandler: Molecular reaction '%s' references invalid charge state Z0=%d for product species '%s'.",
                      def.rateName, p.Z0, p.name
                  );

              if (molecules.IsMolecule(p.name))
                  molecules.GetMass(p.name); // throws if missing
          }


          MolecularReaction reaction;
          reaction.rateName = def.rateName;
          reaction.process = def.process;
          reaction.nReactants = def.nReactants;
          reaction.reactants = def.reactants;
          reaction.nProducts = def.nProducts;
          reaction.products = def.products;
          reaction.temperatureInput = def.temperatureInput;
          reaction.densityInput = def.densityInput;
          reaction.rate = GetMolecularRateByName(def.rateName);

          molecularReactions.push_back(reaction);
      }

      printf(
          "RateHandler: Added %zu molecular reactions.\n",
          molecularReactions.size()
      );

  printf(
          "RateHandler: Added %zu molecular reactions.\n",
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

      
real_t RateHandler::ResolveDensity(
      const MolecularInput& input,
      const len_t ir
  ) {
        if (input.kind == MolecularInputKind::NONE) {
              // Positive coordinate inside the current density grid.
              return 1e19;
        }
        if (input.kind == MolecularInputKind::ELECTRON) {
              return unknowns->GetUnknownData(id_n_cold)[ir];
        }
      else {
          throw std::runtime_error("Invalid density selector");
      }
  }
real_t RateHandler::ResolveTemperature(
      const MolecularInput& input,
      const len_t ir
  ) { ///TODO  add charge temperatures
      if (input.kind == MolecularInputKind::NONE) {
          return 1.0;
      }

      if (input.kind == MolecularInputKind::ELECTRON) {
          return unknowns->GetUnknownData(id_T_cold)[ir];
      }

      if (input.kind == MolecularInputKind::SPECIES) {
          const len_t i = ions->GetIonIndex(input.name1);
          const len_t index = i * grid->GetNr() + ir;

          const real_t W =
              unknowns->GetUnknownData(id_Wi)[index];
          const real_t N =
              unknowns->GetUnknownData(id_Ni)[index];

          return 2.0 * W / (3.0 * Constants::ec * N);
      }

      //if (input.kind == MolecularInputKind::NEUTRAL) {
      //return unknowns->GetUnknownData(id_T_neutral)[ir];
    //}

      if (input.kind == MolecularInputKind::RELATIVE) {
          const len_t i1 = ions->GetIonIndex(input.name1);
          const len_t i2 = ions->GetIonIndex(input.name2);

          const len_t index1 = i1 * grid->GetNr() + ir;
          const len_t index2 = i2 * grid->GetNr() + ir;

          const real_t W1 =
              unknowns->GetUnknownData(id_Wi)[index1];
          const real_t N1 =
              unknowns->GetUnknownData(id_Ni)[index1];

          const real_t W2 =
              unknowns->GetUnknownData(id_Wi)[index2];
          const real_t N2 =
              unknowns->GetUnknownData(id_Ni)[index2];

          const real_t T1 =
              2.0 * W1 / (3.0 * Constants::ec * N1);
          const real_t T2 =
              2.0 * W2 / (3.0 * Constants::ec * N2);

          return 0.5 * (T1 + T2);
      }

      throw std::runtime_error("Invalid temperature selector");
  }


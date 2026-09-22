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
RateHandler::RateHandler(FVM::Grid *grid, IonHandler *ions, 
    ADAS *adas, FVM::UnknownQuantityHandler *unknowns, bool reactionsEnabled, 
    const std::vector<std::string>& enabledReactionNames
) :  grid(grid), ions(ions), adas(adas), unknowns(unknowns) {

    this->unknowns  = unknowns;
    this->id_ions   = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
	this->id_n_cold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
	this->id_Ni  = unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS);
	this->id_T_cold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->id_Wi     = unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER);
    separateNeutrals = unknowns->HasUnknown(OptionConstants::UQTY_WN_ENER);

    if (separateNeutrals)
        id_Wn = unknowns->GetUnknownID(OptionConstants::UQTY_WN_ENER);
    
    //Set rates for IonRateEquation. TODO: Change so other parts of code use Ratehandler as well.
    AddMolecularChargeStateRates();
    AddAtomicChargeStateRates();
      
    if (reactionsEnabled)
          AddMolecularReactionRates(enabledReactionNames);

    //TODO: Add and set rates to runaway ionizationfluid equation 
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
    } 
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
                GetMolecularRateByName("D_2_ionization"));
          } else {
            rates.acd = new ZeroChargeStateRate("zero_ACD"); //this needs to be generalized
            rates.scd = new ZeroChargeStateRate("zero_SCD"); //this needs to be generalized
            }
            
        chargeStateRates[name] = rates;
      }
}


/**
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
            std::string(def.rateName)) == enabledReactionNames.end())
            continue;

        bool allReactantsExist = true;

        for (len_t j = 0; j < def.nReactants; j++) {
            const MolecularReactionSpecies& r = def.reactants[j];

            //We dont care if it has an electron
            if (std::string(r.name) == "e")
                continue;

            //If not all reactants exist (existance is handeled by python)
            if (!ions->HasIon(r.name)) {
                allReactantsExist = false;
                break;
            }
            
            //Check if the charge state is valid for the reactant
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
          }

          // Add the reaction to the list of molecular reactions.
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
          molecularReactions.size());

}

/**
 * Helper method to get a molecular rate by name.
 */
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

/**
 * Resolve the density for a given molecular input. Will help the 
 * molecular rate interpolator to get the correct density for the rate evaluation.
 */
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
/**
 * Resolve the temperature for a given molecular input. Will help the 
 * molecular rate interpolator to get the correct temperature for the rate evaluation.
 */
real_t RateHandler::ResolveTemperature(
      const MolecularInput& input,
      const len_t ir
  ) { 
    //Rate does not depend on temperature
    if (input.kind == MolecularInputKind::NONE) {
        return 1.0;
    }
    //Rate depends on electron temperature
    if (input.kind == MolecularInputKind::ELECTRON) {
        return unknowns->GetUnknownData(id_T_cold)[ir];
    }
    //Rate depends on a specific species temperature
    if (input.kind == MolecularInputKind::SPECIES) {
        const len_t i = ions->GetIonIndex(input.name1);
        return ResolveSpeciesTemperature(ions->GetIonIndex(input.name1), input.charge1, ir
    );}
    //Rate depends on the average temperature of two species
    if (input.kind == MolecularInputKind::RELATIVE) {
        const real_t T1 = ResolveSpeciesTemperature(ions->GetIonIndex(input.name1), input.charge1, ir);
        const real_t T2 = ResolveSpeciesTemperature(ions->GetIonIndex(input.name2), input.charge2, ir);

        return 0.5 * (T1 + T2);
      }

      throw std::runtime_error("Invalid temperature selector");
  }


/**
 * Calculate the temperature of a given species at a given radial position. 
 * This is used to resolve the temperature for molecular rates that depend on the 
 * temperature of a specific species.
 */
real_t RateHandler::ResolveSpeciesTemperature(
      len_t species, int_t charge, len_t ir
) {
    const len_t nr = grid->GetNr();
    const len_t offset = species * nr + ir;

    real_t n, W;

    //If we dont seperate the neutrals we just use the ion density and energy to calculate the temperature
    if (!separateNeutrals) {
        n = unknowns->GetUnknownData(id_Ni)[offset];
        W = unknowns->GetUnknownData(id_Wi)[offset];
      } else {
        //If we do separate the neutrals, we need to consider the neutral density and energy
        const real_t *densities = unknowns->GetUnknownData(id_ions);

        if (charge == 0) {
            n = densities[ions->GetIndex(species, 0) * nr + ir];
            W = unknowns->GetUnknownData(id_Wn)[offset];
        } else {
            n = 0;
            for (len_t z = 1; z <= ions->GetZ(species); z++)
                n += densities[ions->GetIndex(species, z) * nr + ir];

            W = unknowns->GetUnknownData(id_Wi)[offset];
          }
      }

      // Temperature is undefined for an absent population.
      return n > 0 ? W / (1.5 * Constants::ec * n) : 0;
  }



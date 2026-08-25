#include "DREAM/Equations/Fluid/RateHandler.hpp"
#include "FVM/FVMException.hpp"
#include "DREAM/Equations/Fluid/RateData.hpp"
#include "DREAM/MoleculeHandler.hpp"


  using namespace DREAM;

/**
 * Constructor.
  */
    RateHandler::RateHandler(IonHandler *ions, ADAS *adas)
      : ions(ions), adas(adas) {

      AddMolecularChargeStateRates();
      AddAtomicChargeStateRates();
  }

/**
 * Destructor.
 */
  RateHandler::~RateHandler() {
      for (auto it = chargeStateRates.begin(); it != chargeStateRates.end(); it++) {
          delete it->second.acd;
          delete it->second.scd;
      }
  }


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


//PLACEHOLDER ADAS IMPLEMENTATION OF RATES BETWEEN MOLECULAR AND ATOMIC SPECIES



/**
 * 
 * Get the ACD (recombination) charge-state rate for a given molecular species.
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
 * Get the SCD (ionization) charge-state rate for a given molecular species.
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

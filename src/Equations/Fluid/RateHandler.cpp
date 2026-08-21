#include "DREAM/Equations/Fluid/RateHandler.hpp"
#include "FVM/FVMException.hpp"
#include "DREAM/Equations/Fluid/RateData.hpp"


  using namespace DREAM;

/**
 * Constructor.
  */
    RateHandler::RateHandler(IonHandler *ions, ADAS *adas)
      : ions(ions), adas(adas) {

      AddMolecularChargeStateRates();
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

//PLACEHOLDER ADAS IMPLEMENTATION OF ADAS ATOM CHARGE STATE RATES

/**
 * Add charge-state rates for molecular species. 
 * Currently separate from the adas.
 */
  void RateHandler::AddMolecularChargeStateRates() {
    //maybe to this as a for loop over all possible molecules (which are hardcoded)
      if (ions->HasIon("D2")) { //this shouuld ideally not be checked fr name rather that iIon and loop thourgh
        printf("RateHandler::AddMolecularChargeStateRates: Adding molecular charge-state rates for D2.\n");
          ChargeStateRateSet rates;
          rates.acd = new ZeroChargeStateRate("D2_ACD_test");
          rates.scd = new TabledChargeStateRate(
              "D2_SCD_AMJUEL_2.2.9",
              0,  // active Z0: D2^0 -> D2^1
              D2_ioniz_229_nT,
              D2_ioniz_229_nn,
              D2_ioniz_229_coeff
          );

          chargeStateRates["D2"] = rates;
      }

  }

//PLACEHOLDER ADAS IMPLEMENTATION OF RATES BETWEEN MOLECULAR AND ATOMIC SPECIES



/**
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

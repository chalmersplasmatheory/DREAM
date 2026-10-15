/**
 * Implementation of a charge-exchange rate reaction,
 * of the form:
 *
 * A+ + B -> A + B+
 * 
 * This class in very similar to the IonRateEquation class, 
 * but is specifically designed to handle charge-exchange reactions 
 * involving molecular species. .
*/

#include "DREAM/Equations/Fluid/RateHandler.hpp"
#include "DREAM/Equations/Fluid/MoleculeRateReaction.hpp"
#include "DREAM/IonHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/MoleculeHandler.hpp"
#include "DREAM/Constants.hpp"


using namespace DREAM;

/**
 * Constructor
 */
MoleculeRateReaction::MoleculeRateReaction(
    FVM::Grid *g, IonHandler *ihdl,
    const len_t iIon, FVM::UnknownQuantityHandler *unknowns,
    RateHandler *ratehandler
) : IonEquationTerm<FVM::EquationTerm>(g, ihdl, iIon), adas(adas), ratehandler(ratehandler) {
    
    SetName("MoleculeRateReaction");

    this->unknowns  = unknowns;

    this->id_ions   = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
	this->id_n_cold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
	this->id_n_tot  = unknowns->GetUnknownID(OptionConstants::UQTY_N_TOT);
	this->id_T_cold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->id_Wi     = unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER);
    this->id_Ni = unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS);
    separateNeutrals = unknowns->HasUnknown(OptionConstants::UQTY_WN_ENER);
    if (separateNeutrals)
        id_Wn = unknowns->GetUnknownID(OptionConstants::UQTY_WN_ENER);

    AddUnknownForJacobian(unknowns, id_ions);
    AddUnknownForJacobian(unknowns, id_n_cold);
    AddUnknownForJacobian(unknowns, id_T_cold);
    AddUnknownForJacobian(unknowns, id_Ni);
    AddUnknownForJacobian(unknowns, id_Wi);
    if (separateNeutrals)
        AddUnknownForJacobian(unknowns, id_Wn);

    const len_t allocationSize = FindRelevantMolecularReactions();
    AllocateRateCoefficients(allocationSize);

    }
/*
 * Finds the relevant molecular reactions for this ion species and charge state.
 * Returns the number of relevant reactions found.
 */
len_t MoleculeRateReaction::FindRelevantMolecularReactions() {
    const auto& reactions = this->ratehandler->GetMolecularReactions();
    for (const MolecularReaction& reaction : reactions) {
        switch (reaction.process) {
            case MolecularReactionProcess::CHARGE_EXCHANGE:
            case MolecularReactionProcess::DISSOCIATION:
            case MolecularReactionProcess::DISSOCIATION_MULTIPLE:
            case MolecularReactionProcess::RECOMBINATION:
                moleculeReactions.push_back(reaction);
                break;
            default:
                // IONIZATION:               in IonRateEquation (D2 SCD)
                // CHARGE_EXCHANGE_RESONANT: no net particle change
                // *_RUNAWAY:                not supported yet
                break;
        }
    }
    return moleculeReactions.size();
}

/**
 * Allocates memory for the rate coefficients.
 */
void MoleculeRateReaction::AllocateRateCoefficients(const len_t allocationSize) {
    const len_t Nr  = this->grid->GetNr();

    if (allocationSize == 0) return;
    this->Rate = new real_t*[allocationSize];
    this->Rate_N = new real_t*[allocationSize];
    this->Rate_T = new real_t*[allocationSize];

    this->Rate[0] = new real_t[allocationSize * Nr];
    this->Rate_N[0] = new real_t[allocationSize * Nr];
    this->Rate_T[0] = new real_t[allocationSize * Nr];

    for (len_t i = 1; i < allocationSize; i++) {
        this->Rate[i] = this->Rate[i-1] + Nr;
        this->Rate_N[i] = this->Rate_N[i-1] + Nr;
        this->Rate_T[i] = this->Rate_T[i-1] + Nr;
    }
}

/**
 * Destructor.
 */
MoleculeRateReaction::~MoleculeRateReaction() {
    DeallocateRateCoefficients();
}

void MoleculeRateReaction::DeallocateRateCoefficients() {

    if (this->Rate != nullptr) {
      delete [] this->Rate[0];
      delete [] this->Rate;
    }
    if (this->Rate_N != nullptr) {
      delete [] this->Rate_N[0];
      delete [] this->Rate_N;
    }
    if (this->Rate_T != nullptr) {
      delete [] this->Rate_T[0];
        delete [] this->Rate_T;
    }
}

/*
* Rebuild the rate coefficients.
*/
void MoleculeRateReaction::Rebuild(
    const real_t, const real_t, FVM::UnknownQuantityHandler *unknowns
) {
    const len_t Nr = this->grid->GetNr();
     
    for (len_t rateIndex = 0; rateIndex < moleculeReactions.size(); rateIndex++) {
        const MolecularReaction& reaction = moleculeReactions[rateIndex];
      
        for (len_t ir = 0; ir < Nr; ir++) {
            real_t density, temperature;
            density = ratehandler->ResolveDensity(reaction.densityInput, ir);
            temperature = ratehandler->ResolveTemperature(reaction.temperatureInput, ir);
            // here we need to save the ID of density and temp and thats what the jacobian shold depend on
            Rate[rateIndex][ir]   = reaction.rate->Eval(density, temperature);
            Rate_N[rateIndex][ir] = reaction.rate->Eval_deriv_n(density, temperature);
            Rate_T[rateIndex][ir] = reaction.rate->Eval_deriv_T(density, temperature);
        }
    }
}



bool MoleculeRateReaction::SetCSJacobianBlock(
      const len_t uqtyId, const len_t derivId,
      FVM::Matrix *jac, const real_t *x,
      const len_t iIon, const len_t Z0, const len_t rOffset
  ) {
      // This term depends on ion densities, n_cold and T_cold.
    if (!HasJacobianContribution(derivId))
        return false;

    const len_t nr = this->grid->GetNr();
    const std::string name = ions->GetName(iIon);
    const auto& reactions = ratehandler->GetMolecularReactions();
    const real_t *ne = unknowns->GetUnknownData(id_n_cold);

    bool contributes = false;

    for (len_t rateIndex = 0; rateIndex < moleculeReactions.size(); rateIndex++) {
        const MolecularReaction& reaction = moleculeReactions[rateIndex];

        // Net production of the species and charge state
        // whose equation we are currently building.
        real_t net = 0;

        //check if this ion is part of the products, if so, the rate should be positive
        for (len_t j = 0; j < reaction.nProducts; j++) {
            const auto& p = reaction.products[j];

            if (name == p.name && p.Z0 == (int_t)Z0)
                net += p.coefficient;
          }

          //check if this ion is part of the reactants, if so, the rate should be negative
        for (len_t j = 0; j < reaction.nReactants; j++) {
            const auto& r = reaction.reactants[j];

            if (name == r.name && r.Z0 == (int_t)Z0)
                net -= r.coefficient; //coefient is how many of them so if there is 2 produced...
          }

        if (net == 0)
            continue;

        const auto& a = reaction.reactants[0];
        const auto& b = reaction.reactants[1];

        const bool electronA = std::string(a.name) == "e";
        const bool electronB = std::string(b.name) == "e";

        const len_t idA = electronA ? id_n_cold : id_ions;
        const len_t idB = electronB ? id_n_cold : id_ions;

        const len_t offsetA = electronA ? 0 : ions->GetIndex(ions->GetIonIndex(a.name), a.Z0) * nr;
        const len_t offsetB = electronB ? 0 : ions->GetIndex(ions->GetIonIndex(b.name), b.Z0) * nr;
        
        const MolecularInput& Tin = reaction.temperatureInput;
        for (len_t ir = 0; ir < nr; ir++) {
            const len_t row = rOffset + ir;

            const real_t nA = electronA ? ne[ir] : x[offsetA + ir];
            const real_t nB = electronB ? ne[ir] : x[offsetB + ir];
            //Dependance on the species densities, nA and nB.
            // dS/dnA  →  column n_cold (electron) or n_i (ion)
            if (derivId == idA)
                jac->SetElement(row, offsetA + ir, net * Rate[rateIndex][ir] * nB);

            // dS/dnB
            if (derivId == idB)
                jac->SetElement(row, offsetB + ir, net * Rate[rateIndex][ir] * nA);
            
                //If the input value depends on the electron density, then set the jacobian element for the electron density.
            if (derivId == id_n_cold && reaction.densityInput.kind == MolecularInputKind::ELECTRON) {
                jac->SetElement(row, ir, net * Rate_N[rateIndex][ir] * nA * nB);
            }
            //If the input value depends on the electron temperature, then set the jacobian element for the electron temperature.
            if (derivId == id_T_cold && reaction.temperatureInput.kind == MolecularInputKind::ELECTRON) {
                // The rate coefficient K depends on T_cold.
                jac->SetElement(row, ir, net * Rate_T[rateIndex][ir] * nA * nB);
            } 
            //If the input temperature depends on the species temperature, then set the jacobian element for the species temperature.
            if (Tin.kind == MolecularInputKind::SPECIES) {
                const len_t iz      = ions->GetIonIndex(Tin.name1);
                const bool neutral  = (Tin.charge1 == 0); //check if the input temperature is for a neutral species
                const len_t energyId = (separateNeutrals && neutral) ? id_Wn : id_Wi;

                // Same density as RateHandler::ResolveSpeciesTemperature
                const len_t firstZ = neutral ? 0 : 1;
                const len_t lastZ  = neutral ? 0 : ions->GetZ(iz);
                real_t n = 0;
                if (!separateNeutrals)
                    n = unknowns->GetUnknownData(id_Ni)[iz*nr + ir];
                else
                    for (len_t Z = firstZ; Z <= lastZ; Z++)
                        n += x[ions->GetIndex(iz, Z)*nr + ir];

                if (n > 0) {
                    const real_t T = ratehandler->ResolveTemperature(Tin, ir);

                    // dS/dW
                    if (derivId == energyId)
                        jac->SetElement(row, iz*nr + ir, net * Rate_T[rateIndex][ir] * nA * nB / (1.5*Constants::ec*n));

                    // dS/dn
                    if (!separateNeutrals && derivId == id_Ni)
                        jac->SetElement(row, iz*nr + ir, -(net * Rate_T[rateIndex][ir] * nA * nB) * T / n);
                    else if (separateNeutrals && derivId == id_ions)
                        for (len_t Z = firstZ; Z <= lastZ; Z++)
                            jac->SetElement(row, ions->GetIndex(iz, Z)*nr + ir, -(net * Rate_T[rateIndex][ir] * nA * nB) * T / n);
                }
            }


        }
          contributes = true;
    }
      return contributes;
}



void MoleculeRateReaction::SetCSMatrixElements(
      FVM::Matrix*, real_t*, const len_t, const len_t, const len_t
  ) { //TODO: implement this function if needed
  }

void MoleculeRateReaction::SetCSVectorElements(
      real_t *vec, const real_t *x, const len_t iIon, const len_t Z0, const len_t rOffset
  ) {
    const len_t nr = this->grid->GetNr();
    const std::string name = ions->GetName(iIon);
    const real_t *ne = unknowns->GetUnknownData(id_n_cold);

    for (len_t rateIndex = 0; rateIndex < moleculeReactions.size(); rateIndex++) { //go though all CS reactions and find the ones that are charge exchange
        const MolecularReaction& reaction = moleculeReactions[rateIndex];

        // Net number of particles produced in this equation. Check if we should add 
        // or remove density for this species and charge state.
        real_t net = 0;

        //check if this ion is part of the products, if so, the rate should be positive
        for (len_t j = 0; j < reaction.nProducts; j++) {
            const auto& p = reaction.products[j];

            if (name == p.name && p.Z0 == (int_t)Z0)
                net += p.coefficient;
        }

        //check if this ion is part of the reactants, if so, the rate should be negative
        for (len_t j = 0; j < reaction.nReactants; j++) {
            const auto& r = reaction.reactants[j];

            if (name == r.name && r.Z0 == (int_t)Z0)
                net -= r.coefficient;
        }

        if (net == 0)
            continue;
            
        const auto& a = reaction.reactants[0];
        const auto& b = reaction.reactants[1];
          
        const bool electronA = std::string(a.name) == "e";
        const bool electronB = std::string(b.name) == "e";

        const len_t offsetA = electronA ? 0 : ions->GetIndex(ions->GetIonIndex(a.name), a.Z0) * nr;
        const len_t offsetB = electronB ? 0 : ions->GetIndex(ions->GetIonIndex(b.name), b.Z0) * nr;
          
        for (len_t ir = 0; ir < nr; ir++) {
            const real_t nA = electronA ? ne[ir] : x[offsetA + ir];
            const real_t nB = electronB ? ne[ir] : x[offsetB + ir];
            const real_t R = Rate[rateIndex][ir] * nA * nB;

            vec[rOffset + ir] += net * R;
        }
      }
  }

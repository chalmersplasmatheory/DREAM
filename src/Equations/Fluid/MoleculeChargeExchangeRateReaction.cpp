



#include "DREAM/Equations/Fluid/RateHandler.hpp"
#include "DREAM/Equations/Fluid/MoleculeChargeExchangeRateReaction.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/Equations/Fluid/MoleculeChargeExchangeRateReaction.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/MoleculeHandler.hpp"


using namespace DREAM;

MoleculeChargeExchangeRateReaction::MoleculeChargeExchangeRateReaction(
    FVM::Grid *g, IonHandler *ihdl,
    const len_t iIon,ADAS *adas, FVM::UnknownQuantityHandler *unknowns,
    RateHandler *ratehandler, bool addFluidIonization, bool addFluidJacobian, bool isAbl 
) : IonEquationTerm<FVM::EquationTerm>(g, ihdl, iIon), adas(adas), 
    ratehandler(ratehandler), addFluidIonization(addFluidIonization), addFluidJacobian(addFluidJacobian) {
    
        SetName("MoleculeChargeExchangeRateReaction");
        printf("Species name called in MoleculeChargeExchangeRateReaction: %s\n", this->ions->GetName(iIon).c_str());

        this->unknowns  = unknowns;

        this->id_ions   = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
		this->id_n_cold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
		this->id_n_tot  = unknowns->GetUnknownID(OptionConstants::UQTY_N_TOT);
		this->id_T_cold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);

        const len_t allocationSize = FindRelevantMolecularReactions();
        AllocateRateCoefficients(allocationSize);

    }

len_t MoleculeChargeExchangeRateReaction::FindRelevantMolecularReactions() {
    const auto& reactions = this->ratehandler->GetMolecularReactions();
    len_t ReactionsCount = 0;
    for (const MolecularReaction& reaction : reactions) {
        if (reaction.process == MolecularReactionProcess::CHARGE_EXCHANGE)  {
        chargeExchangeReactions.push_back(reaction);
        ReactionsCount++;
        }
    }
    return ReactionsCount;
}

void MoleculeChargeExchangeRateReaction::AllocateRateCoefficients(const len_t allocationSize) {
    const len_t Nr  = this->grid->GetNr();
    //the size should be all charge exchange reaction, not Zion+1
    //FindRelevantMolecularReactions should have been called before this function to find the relevant reactions
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
MoleculeChargeExchangeRateReaction::~MoleculeChargeExchangeRateReaction() {
    DeallocateRateCoefficients();
}

void MoleculeChargeExchangeRateReaction::DeallocateRateCoefficients() {

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

void MoleculeChargeExchangeRateReaction::Rebuild(
    const real_t, const real_t, FVM::UnknownQuantityHandler *unknowns
) {
    const len_t Nr = this->grid->GetNr();

    real_t *T = unknowns->GetUnknownData(id_T_cold);
    real_t *n = unknowns->GetUnknownData(id_n_cold);
    
    for (len_t rateIndex = 0; rateIndex < chargeExchangeReactions.size(); rateIndex++) {
        const MolecularReaction& reaction = chargeExchangeReactions[rateIndex];
      
        for (len_t ir = 0; ir < Nr; ir++) {
            Rate[rateIndex][ir]   = reaction.rate->Eval(n[ir], T[ir]);
            Rate_N[rateIndex][ir] = reaction.rate->Eval_deriv_n(n[ir], T[ir]);
            Rate_T[rateIndex][ir] = reaction.rate->Eval_deriv_T(n[ir], T[ir]);
        }
    }
}

 bool MoleculeChargeExchangeRateReaction::SetCSJacobianBlock(
      const len_t uqtyId, const len_t derivId,
      FVM::Matrix *jac, const real_t *x,
      const len_t iIon, const len_t Z0, const len_t rOffset
  ) {
      // This term depends on ion densities, n_cold and T_cold.
    if (derivId != id_ions &&
          derivId != id_n_cold &&
          derivId != id_T_cold)
          return false;

    const len_t nr = this->grid->GetNr();
    const std::string name = ions->GetName(iIon);
    const auto& reactions = ratehandler->GetMolecularReactions();

    bool contributes = false;

    for (len_t rateIndex = 0; rateIndex < chargeExchangeReactions.size(); rateIndex++) {
        const MolecularReaction& reaction = chargeExchangeReactions[rateIndex];

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

        // Current charge-exchange reactions should only have 2 reactants
        const auto& a = reaction.reactants[0];
        const auto& b = reaction.reactants[1];

        const len_t ia = ions->GetIonIndex(a.name);
        const len_t ib = ions->GetIonIndex(b.name);

        const len_t offsetA = ions->GetIndex(ia, a.Z0) * nr;
        const len_t offsetB = ions->GetIndex(ib, b.Z0) * nr;

        for (len_t ir = 0; ir < nr; ir++) {
            const len_t row = rOffset + ir;

            const real_t nA = x[offsetA + ir];
            const real_t nB = x[offsetB + ir];

            if (derivId == id_ions) {
                // d(net * K * nA * nB) / dnA
                jac->SetElement(
                      row, offsetA + ir,
                      net * Rate[rateIndex][ir] * nB
                  );

                // d(net * K * nA * nB) / dnB
                jac->SetElement(
                      row, offsetB + ir,
                      net * Rate[rateIndex][ir] * nA
                  );
            }

            if (derivId == id_n_cold) {
                // The rate coefficient K depends on n_cold.
                jac->SetElement(
                      row, ir,
                      net * Rate_N[rateIndex][ir] * nA * nB
                  );
            }

            if (derivId == id_T_cold) {
                // The rate coefficient K depends on T_cold.
                jac->SetElement(
                      row, ir,
                      net * Rate_T[rateIndex][ir] * nA * nB
                  );
            } 
        }
          contributes = true;
    }
      return contributes;
}


  void MoleculeChargeExchangeRateReaction::SetCSMatrixElements(
      FVM::Matrix*, real_t*, const len_t, const len_t, const len_t
  ) {
  }

  void MoleculeChargeExchangeRateReaction::SetCSVectorElements(
      real_t *vec, const real_t *x, const len_t iIon, const len_t Z0, const len_t rOffset
  ) {
    const len_t nr = this->grid->GetNr();
    const std::string name = ions->GetName(iIon);
    for (len_t rateIndex = 0; rateIndex < chargeExchangeReactions.size(); rateIndex++) { //go though all CS reactions and find the ones that are charge exchange
        const MolecularReaction& reaction = chargeExchangeReactions[rateIndex];

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
            
        printf("Net production of %s (Z0=%d) in reaction %s: %g\n", name.c_str(), Z0, reaction.rateName, net);
        const auto& a = reaction.reactants[0];
        const auto& b = reaction.reactants[1];

        const len_t ia = ions->GetIonIndex(a.name);
        const len_t ib = ions->GetIonIndex(b.name);

        const len_t offsetA = ions->GetIndex(ia, a.Z0) * nr;
        const len_t offsetB = ions->GetIndex(ib, b.Z0) * nr;
          
        for (len_t ir = 0; ir < nr; ir++) {
            const real_t R =
                Rate[rateIndex][ir] *
                x[offsetA + ir] *
                x[offsetB + ir];
                printf("Added density for reaction %s at radial index %d: %g\n", reaction.rateName, ir, R);

            vec[rOffset + ir] += net * R;
          }
      }
  }

    
  



////then figure out how to set the matrix elements

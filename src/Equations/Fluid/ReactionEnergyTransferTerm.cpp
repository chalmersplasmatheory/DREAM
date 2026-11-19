/**
 * Implementation of the energy given to (or taken from) the heavy particles
 * of one ion species (W_i or W_n) in molecular reactions.
 *
 * (The thermal energy that atoms carry between W_n and W_i when they are
 * ionized or recombine is in NeutralIonizationEnergyTransferTerm.)
 */
#include "DREAM/Equations/Fluid/ReactionEnergyTransferTerm.hpp"
#include "DREAM/Constants.hpp"
#include <algorithm>
#include <string>

using namespace DREAM;

/**
 * Constructor.
 */
ReactionEnergyTransferTerm::ReactionEnergyTransferTerm(
    FVM::Grid *grid, len_t iz, FVM::UnknownQuantityHandler *unknowns,
    IonHandler *ions, RateHandler *ratehandler, bool neutralEquation,
    bool includeFluidREIonization
) : FVM::EquationTerm(grid), unknowns(unknowns), ions(ions),
    ratehandler(ratehandler), iz(iz), neutralEquation(neutralEquation) {

    SetName("ReactionEnergyTransferTerm");

    this->id_Te = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->separateNeutrals = unknowns->HasUnknown(OptionConstants::UQTY_WN_ENER);

    AddUnknownForJacobian(unknowns, unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES));
    AddUnknownForJacobian(unknowns, unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD));
    AddUnknownForJacobian(unknowns, unknowns->GetUnknownID(OptionConstants::UQTY_N_RE));
    AddUnknownForJacobian(unknowns, id_Te);
    AddUnknownForJacobian(unknowns, unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER));
    AddUnknownForJacobian(unknowns, unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS));
    if (separateNeutrals) {
        AddUnknownForJacobian(unknowns, unknowns->GetUnknownID(OptionConstants::UQTY_WN_ENER));
        AddUnknownForJacobian(unknowns, unknowns->GetUnknownID(OptionConstants::UQTY_NN_DENS));
    }

    this->iIonName = ions->GetName(iz);
    this->iIonMass = ions->GetIonSpeciesMass(iz);

    const std::vector<MolecularReaction>& allReactions = ratehandler->GetMolecularReactions();

    // Sort out relevant reactions for this species
    // For each reaction, determine whether this species is a product or reactant,
    // and if it is a product, what fraction of the energy released in the reaction
    // goes to this species (by momentum conservation).
    for (len_t k = 0; k < allReactions.size(); k++) {
        const MolecularReaction &reaction = allReactions[k];
        if (reaction.process == MolecularReactionProcess::DISSOCIATION_RUNAWAY)
            continue;
        if (reaction.process == MolecularReactionProcess::IONIZATION_RUNAWAY && !includeFluidREIonization)
            continue;

        real_t M = 0;
        len_t nHeavyProducts = 0;
        for (len_t j = 0; j < reaction.nProducts; j++) {
            const MolecularReactionSpecies &p = reaction.products[j];
            if (std::string(p.name) == "e" || std::string(p.name) == "RE")
                continue;
            M += p.coefficient * ions->GetIonSpeciesMass(ions->GetIonIndex(p.name));
            nHeavyProducts++;
        }

        real_t s = 0;
        for (len_t j = 0; j < reaction.nProducts; j++) {
            const MolecularReactionSpecies &p = reaction.products[j];
            const bool neutral = separateNeutrals && p.Z0 == 0;
            if (this->iIonName != p.name || neutral != neutralEquation)
                continue;

            if (nHeavyProducts == 1)
                s += 1;
            else
                s += (M - p.coefficient * this->iIonMass) / M;
        }

        bool reactant = false;
        for (len_t j = 0; j < reaction.nReactants; j++) {
            const MolecularReactionSpecies &r = reaction.reactants[j];
            const bool neutral = separateNeutrals && r.Z0 == 0;
            if (this->iIonName != r.name || neutral != neutralEquation)
                reactant = true;
        }

        if (s != 0 || reactant) {
            this->reactions.push_back(k);
            this->share.push_back(s);
        }
    }

    this->Energy = new real_t[reactions.size() * grid->GetNr()];
}

/**
 * Destructor.
 */
ReactionEnergyTransferTerm::~ReactionEnergyTransferTerm() {
    delete [] this->Energy;
}

/**
 * Rebuild the energy given to this species per reaction.
 */
void ReactionEnergyTransferTerm::Rebuild(
    const real_t, const real_t, FVM::UnknownQuantityHandler*
) {
    const len_t nr = grid->GetNr();
    const real_t ec = Constants::ec;
    const real_t *Te = unknowns->GetUnknownData(id_Te);

    for (len_t i = 0; i < reactions.size(); i++) {
        const MolecularReaction &reaction = ratehandler->GetMolecularReactions()[reactions[i]];

        // Energy released to the heavy particles [eV]. For electron impact the
        // electron pays an endothermic reaction.
        bool electronImpact = false;
        for (len_t j = 0; j < reaction.nReactants; j++)
            if (std::string(reaction.reactants[j].name) == "e" || std::string(reaction.reactants[j].name) == "RE")
                electronImpact = true;
        real_t Ereleased = reaction.energyReleased;
        if (electronImpact)
            Ereleased = std::max(Ereleased, 0.0) + reaction.franckCondonEnergy;

        for (len_t ir = 0; ir < nr; ir++) {
            real_t E = share[i] * ec * Ereleased;
            if (reaction.electronsConsumed > 0)
                E += share[i] * 1.5 * ec * reaction.electronsConsumed * Te[ir];

            // Thermal energy of the heavy reactants: given to the products (share),
            // and taken from the reactants of this species
            for (len_t j = 0; j < reaction.nReactants; j++) {
                const MolecularReactionSpecies &r = reaction.reactants[j];
                if (std::string(r.name) == "e" || std::string(r.name) == "RE")
                    continue;

                const real_t T = ratehandler->ResolveSpeciesTemperature(ions->GetIonIndex(r.name), r.Z0, ir);
                const bool neutral = separateNeutrals && r.Z0 == 0;
                E += share[i] * 1.5 * ec * r.coefficient * T;
                if (this->iIonName == r.name && neutral == neutralEquation)
                    E -= 1.5 * ec * r.coefficient * T;
            }
            Energy[i*nr + ir] = E;
        }
    }
}

/**
 * Build function vector.
 */
void ReactionEnergyTransferTerm::SetVectorElements(
    real_t *vec, const real_t*
) {
    const len_t nr = grid->GetNr();

    for (len_t i = 0; i < reactions.size(); i++) {
        const real_t *R = ratehandler->GetReactionRate(reactions[i]);
        for (len_t ir = 0; ir < nr; ir++)
            vec[iz*nr + ir] += R[ir] * Energy[i*nr + ir];
    }
}

/**
 * Build block of Jacobian matrix of Q = R(n, T_rate, nA, nB) * E(T_cold, T_j).
 *
 * derivId: ID of unknown quantity with respect to which differentiation
 *          should be carried out.
 * jac:     Jacobian matrix to build.
 */
bool ReactionEnergyTransferTerm::SetJacobianBlock(
    const len_t, const len_t derivId, FVM::Matrix *jac, const real_t*
) {
    if (!HasJacobianContribution(derivId))
        return false;

    const len_t nr = grid->GetNr();
    const real_t ec = Constants::ec;
    bool contributes = false;

    for (len_t i = 0; i < reactions.size(); i++) {
        const len_t k = reactions[i];
        const MolecularReaction &reaction = ratehandler->GetMolecularReactions()[k];
        const real_t *R = ratehandler->GetReactionRate(k);

        for (len_t ir = 0; ir < nr; ir++) {
            const len_t row = iz*nr + ir;

            bool added = ratehandler->SetReactionRateJacobian(k, derivId, jac, row, ir, Energy[i*nr + ir]);
            if (added)
                contributes = true;

            if (derivId == id_Te && reaction.electronsConsumed > 0) {
                jac->SetElement(row, ir, R[ir] * share[i] * 1.5 * ec * reaction.electronsConsumed);
                contributes = true;
            }
            for (len_t j = 0; j < reaction.nReactants; j++) {
                const MolecularReactionSpecies &r = reaction.reactants[j];
                if (std::string(r.name) == "e" || std::string(r.name) == "RE")
                    continue;

                const bool neutral = separateNeutrals && r.Z0 == 0;
                real_t dEdT = share[i] * 1.5 * ec * r.coefficient;
                if (this->iIonName == r.name && neutral == neutralEquation)
                    dEdT = dEdT - 1.5 * ec * r.coefficient;

                added = ratehandler->SetSpeciesTemperatureJacobian(
                    derivId, jac, row, ir, ions->GetIonIndex(r.name), r.Z0, R[ir] * dEdT
                );
                if (added)
                    contributes = true;
            }
        }
    }

    return contributes;
}
/**
 * Implementation of the molecular reactions in the ion rate equation
 */

#include "DREAM/Equations/Fluid/MoleculeRateReaction.hpp"
#include <string>


using namespace DREAM;

/**
 * Constructor.
 */
MoleculeRateReaction::MoleculeRateReaction(
    FVM::Grid *g, IonHandler *ihdl,
    const len_t iIon, FVM::UnknownQuantityHandler *unknowns,
    RateHandler *ratehandler, bool includeFluidREIonization
) : IonEquationTerm<FVM::EquationTerm>(g, ihdl, iIon), ratehandler(ratehandler) {

    SetName("MoleculeRateReaction");

    this->unknowns = unknowns;
    this->id_ions  = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);


    AddUnknownForJacobian(unknowns, id_ions);
    AddUnknownForJacobian(unknowns, unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD));
    AddUnknownForJacobian(unknowns, unknowns->GetUnknownID(OptionConstants::UQTY_N_RE));
    AddUnknownForJacobian(unknowns, unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD));
    AddUnknownForJacobian(unknowns, unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS));
    AddUnknownForJacobian(unknowns, unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER));
    if (unknowns->HasUnknown(OptionConstants::UQTY_WN_ENER)) {
        AddUnknownForJacobian(unknowns, unknowns->GetUnknownID(OptionConstants::UQTY_WN_ENER));
        AddUnknownForJacobian(unknowns, unknowns->GetUnknownID(OptionConstants::UQTY_NN_DENS));
    }

    // Sort ion relevant reactions for this iIon
    const std::vector<MolecularReaction>& allReactions = ratehandler->GetMolecularReactions();
    for (len_t k = 0; k < allReactions.size(); k++) {
        const MolecularReaction &reaction = allReactions[k];
        if (reaction.process == MolecularReactionProcess::CHARGE_EXCHANGE_RESONANT ||
            reaction.process == MolecularReactionProcess::DISSOCIATION_RUNAWAY)
            continue;
        if (reaction.process == MolecularReactionProcess::IONIZATION_RUNAWAY && !includeFluidREIonization)
            continue;

        bool involved = false;
        for (len_t j = 0; j < reaction.nReactants; j++)
            if (reaction.reactantIz[j] == iIon)
                involved = true;
        for (len_t j = 0; j < reaction.nProducts; j++)
            if (reaction.productIz[j] == iIon)
                involved = true;

        if (involved)
            this->reactions.push_back(k);
    }
}

bool MoleculeRateReaction::SetCSJacobianBlock(
    const len_t /*uqtyId*/, const len_t derivId, FVM::Matrix *jac, const real_t*,
    const len_t /*iIon*/, const len_t Z0, const len_t rOffset
) {
    if (!HasJacobianContribution(derivId))
        return false;

    const len_t nr = this->grid->GetNr();
    bool contributes = false;

    for (len_t k : reactions) {
        const MolecularReaction &reaction = ratehandler->GetMolecularReactions()[k];

        // Check if more than one reactant is the unknown with respect to which we are differentiating.
        real_t net = 0;
        for (len_t j = 0; j < reaction.nProducts; j++)
            if (reaction.productIz[j] == iIon && reaction.products[j].Z0 == (int_t)Z0)
                net += reaction.products[j].coefficient;
        for (len_t j = 0; j < reaction.nReactants; j++)
            if (reaction.reactantIz[j] == iIon && reaction.reactants[j].Z0 == (int_t)Z0)
                net -= reaction.reactants[j].coefficient;
        if (net == 0)
            continue;

        for (len_t ir = 0; ir < nr; ir++) {
            bool added = ratehandler->SetReactionRateJacobian(k, derivId, jac, rOffset + ir, ir, net);
            if (added)
                contributes = true;
        }
    }
    return contributes;
}


void MoleculeRateReaction::SetCSMatrixElements(
    FVM::Matrix *mat, real_t *rhs, const len_t /*iIon*/, const len_t Z0, const len_t rOffset
) {
    const len_t nr = this->grid->GetNr();

    for (len_t k : reactions) {
        const MolecularReaction &reaction = ratehandler->GetMolecularReactions()[k];

        real_t net = 0;
        for (len_t j = 0; j < reaction.nProducts; j++)
            if (reaction.productIz[j] == iIon && reaction.products[j].Z0 == (int_t)Z0)
                net += reaction.products[j].coefficient;
        for (len_t j = 0; j < reaction.nReactants; j++)
            if (reaction.reactantIz[j] == iIon && reaction.reactants[j].Z0 == (int_t)Z0)
                net -= reaction.reactants[j].coefficient;
        if (net == 0)
            continue;

        const real_t *K = ratehandler->GetRateCoefficient(k);
        const len_t idA = reaction.reactantId[0], offsetA = reaction.reactantOffset[0];
        const len_t idB = reaction.reactantId[1], offsetB = reaction.reactantOffset[1];
        const real_t *nA = unknowns->GetUnknownData(idA) + offsetA;
        const real_t *nB = unknowns->GetUnknownData(idB) + offsetB;

        for (len_t ir = 0; ir < nr; ir++) {
            const bool useB = (idB == id_ions) && (idA != id_ions || nB[ir] < nA[ir]);

            if (useB)
                mat->SetElement(rOffset + ir, offsetB + ir, net * K[ir] * nA[ir]);
            else if (idA == id_ions)
                mat->SetElement(rOffset + ir, offsetA + ir, net * K[ir] * nB[ir]);
            // No ion reactant: the whole term is a known source
            else if (rhs != nullptr)
                rhs[rOffset + ir] += net * K[ir] * nA[ir] * nB[ir];
        }
    }
}

void MoleculeRateReaction::SetCSVectorElements(
    real_t *vec, const real_t*, const len_t /*iIon*/, const len_t Z0, const len_t rOffset
) {
    const len_t nr = this->grid->GetNr();

    for (len_t k : reactions) {
        const MolecularReaction &reaction = ratehandler->GetMolecularReactions()[k];

        real_t net = 0;
        for (len_t j = 0; j < reaction.nProducts; j++)
            if (reaction.productIz[j] == iIon && reaction.products[j].Z0 == (int_t)Z0)
                net += reaction.products[j].coefficient;
        for (len_t j = 0; j < reaction.nReactants; j++)
            if (reaction.reactantIz[j] == iIon && reaction.reactants[j].Z0 == (int_t)Z0)
                net -= reaction.reactants[j].coefficient;
        if (net == 0)
            continue;

        const real_t *R = ratehandler->GetReactionRate(k);
        for (len_t ir = 0; ir < nr; ir++)
            vec[rOffset + ir] += net * R[ir];
    }
}

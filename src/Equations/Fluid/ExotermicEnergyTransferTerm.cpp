/** Exotermic heating using the birth-temperature prescription in HollmannModel. */
#include "DREAM/Equations/Fluid/ExotermicEnergyTransferTerm.hpp"
#include "DREAM/Constants.hpp"
#include <stdexcept>
#include <string>

using namespace DREAM;

/*
* Constructor
*/
ExotermicEnergyTransferTerm::ExotermicEnergyTransferTerm(
    FVM::Grid *grid, len_t iz, FVM::UnknownQuantityHandler *unknowns,
    IonHandler *ions, RateHandler *ratehandler
) : FVM::EquationTerm(grid), unknowns(unknowns), ions(ions),
    ratehandler(ratehandler), iz(iz) {

    SetName("ExotermicEnergyTransferTerm");
    mi = ions->GetIonSpeciesMass(iz);
    
    id_ions = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    id_Wn = unknowns->GetUnknownID(OptionConstants::UQTY_WN_ENER);
    id_ne = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    id_Te = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    
    AddUnknownForJacobian(unknowns, id_ions);
    AddUnknownForJacobian(unknowns, id_Wn);
    AddUnknownForJacobian(unknowns, id_ne);
    AddUnknownForJacobian(unknowns, id_Te);
}

/** Total released energy [eV], from Table 4 of HollmannModel. */
real_t ExotermicEnergyTransferTerm::getExotermicEnergy(
    const MolecularReaction& reaction
) const {
    //I think the easies thing is to have this here?
    struct ExothermicReaction { const char *name; real_t energy; };
    static const ExothermicReaction energies[] = {
        {"D_2_dissociation", 5.0},
        {"D_2_dissociation_2", 6.0},
        {"D_2+_dissociation", 7.0},
        {"D_2+_D_2_dissociation", 1.7},
        {"D_3+_dissociation", 9.2},
        {"D_3+_dissociation_2", 4.8},
        {"D_2+_Ar_dissociation", 1.8},
        {"D_3+_Ar_dissociation", 1.1},
        {"Ar+_D_2_dissociation", 1.5},
        {"Ar+_D_2_charge_exchange", 1.7}
    };
    const std::string name = reaction.rateName;
    for (const auto& entry : energies)
        if (name == entry.name)
            return entry.energy;
    return 0;
}

/** Total product mass [kg], including multiplicities and charged products. */
real_t ExotermicEnergyTransferTerm::getTotalMassReactionProducts(
    const MolecularReaction& reaction
) const {
    real_t mass = 0;
    for (len_t j = 0; j < reaction.nProducts; j++) {
        const auto& product = reaction.products[j];
        real_t productMass;
        if (std::string(product.name) == "e") {
            productMass = Constants::me;
        } else {
            const len_t speciesIndex = ions->GetIonIndex(product.name);
            productMass = ions->GetIonSpeciesMass(speciesIndex);
        }
        mass += product.coefficient * productMass;
    }
    return mass;
}

/** Add the reaction heating to the neutral energy equation. */
void ExotermicEnergyTransferTerm::SetVectorElements(
    real_t *vec, const real_t*
) {
    const len_t nr = grid->GetNr();
    const real_t *ni = unknowns->GetUnknownData(id_ions);
    const real_t *ne = unknowns->GetUnknownData(id_ne);
    const std::string name = ions->GetName(iz);

    for (const auto& reaction : ratehandler->GetMolecularReactions()) {
        const real_t deltaE = getExotermicEnergy(reaction);
        if (deltaE == 0) /// if to exotermix reaction exist
            continue;
        len_t multiplicity = 0; //in case of multiple same products
        for (len_t j = 0; j < reaction.nProducts; j++) {
            const auto& product = reaction.products[j];
            if (name == product.name && product.Z0 == 0)
                multiplicity += product.coefficient;
        }
        if (multiplicity == 0)
            continue;

        const real_t M = getTotalMassReactionProducts(reaction);
        const real_t Tbirth = (M - mi) / M * deltaE;
        const auto& a = reaction.reactants[0];
        const auto& b = reaction.reactants[1];
        // Locate each reactant in its density unknown.
        len_t densityIdA = id_ions;
        len_t densityIdB = id_ions;
        len_t offsetA = 0;
        len_t offsetB = 0;

        if (std::string(a.name) == "e") {
            densityIdA = id_ne;
        } else {
            const len_t speciesA = ions->GetIonIndex(a.name);
            offsetA = ions->GetIndex(speciesA, a.Z0) * nr;
        }

        if (std::string(b.name) == "e") {
            densityIdB = id_ne;
        } else {
            const len_t speciesB = ions->GetIonIndex(b.name);
            offsetB = ions->GetIndex(speciesB, b.Z0) * nr;
        }

        for (len_t ir = 0; ir < nr; ir++) {
            const len_t row = iz * nr + ir;
            real_t nA;
            real_t nB;
            if (densityIdA == id_ne) {
                nA = ne[ir];
            } else {
                nA = ni[offsetA + ir];
            }
            if (densityIdB == id_ne) {
                nB = ne[ir];
            } else {
                nB = ni[offsetB + ir];
            }

            const real_t density = ratehandler->ResolveDensity(reaction.densityInput, ir);
            const real_t T = ratehandler->ResolveTemperature(reaction.temperatureInput, ir);
            const real_t K = reaction.rate->Eval(density, T);

            // K*nA*nB is the number of reaction events per m^3 per second.
            const real_t reactionRate = K * nA * nB;
            const real_t productionRate = multiplicity * reactionRate;
            const real_t power = 1.5 * Constants::ec * productionRate * Tbirth;

            vec[row] += power;
        }
    }
}

/**
 * Differentiate Q = energyPerReaction*K*nA*nB.
 * The blocks below handle reactant densities, the rate's electron-density
 * input, its electron-temperature input, and its neutral-temperature input.
 */
bool ExotermicEnergyTransferTerm::SetJacobianBlock(
    const len_t, const len_t derivId, FVM::Matrix *jac, const real_t*
) {
    if (derivId != id_ions && derivId != id_ne &&
        derivId != id_Te && derivId != id_Wn) {
        return false;
    }

    const len_t nr = grid->GetNr();
    const real_t *ni = unknowns->GetUnknownData(id_ions);
    const real_t *ne = unknowns->GetUnknownData(id_ne);
    const std::string name = ions->GetName(iz);
    bool contributes = false;

    for (const auto& reaction : ratehandler->GetMolecularReactions()) {
        const real_t deltaE = getExotermicEnergy(reaction);
        if (deltaE == 0)
            continue;
        len_t multiplicity = 0;
        for (len_t j = 0; j < reaction.nProducts; j++) {
            const auto& product = reaction.products[j];
            if (name == product.name && product.Z0 == 0)
                multiplicity += product.coefficient;
        }
        if (multiplicity == 0)
            continue;

        // Supported by the current ten-reaction table.
        const auto& temperature = reaction.temperatureInput;
        if (temperature.kind != MolecularInputKind::NONE &&
            temperature.kind != MolecularInputKind::ELECTRON &&
            !(temperature.kind == MolecularInputKind::SPECIES && temperature.charge1 == 0))
            throw std::runtime_error("Exothermic term: unsupported temperature selector.");

        const real_t M = getTotalMassReactionProducts(reaction);
        const real_t Tbirth = (M - mi) / M * deltaE;
        const real_t energyPerReaction = 1.5 * Constants::ec * multiplicity * Tbirth;
        const auto& a = reaction.reactants[0];
        const auto& b = reaction.reactants[1];
        // Locate each reactant in its density unknown.
        len_t densityIdA = id_ions;
        len_t densityIdB = id_ions;
        len_t offsetA = 0;
        len_t offsetB = 0;

        if (std::string(a.name) == "e") {
            densityIdA = id_ne;
        } else {
            const len_t speciesA = ions->GetIonIndex(a.name);
            offsetA = ions->GetIndex(speciesA, a.Z0) * nr;
        }

        if (std::string(b.name) == "e") {
            densityIdB = id_ne;
        } else {
            const len_t speciesB = ions->GetIonIndex(b.name);
            offsetB = ions->GetIndex(speciesB, b.Z0) * nr;
        }

        for (len_t ir = 0; ir < nr; ir++) {
            const len_t row = iz * nr + ir;
            real_t nA;
            real_t nB;
            if (densityIdA == id_ne) {
                nA = ne[ir];
            } else {
                nA = ni[offsetA + ir];
            }
            if (densityIdB == id_ne) {
                nB = ne[ir];
            } else {
                nB = ni[offsetB + ir];
            }

            const real_t density = ratehandler->ResolveDensity(reaction.densityInput, ir);
            const real_t T = ratehandler->ResolveTemperature(temperature, ir);
            const real_t K = reaction.rate->Eval(density, T);

            // Explicit density factors. SetElement adds contributions when
            // reactants and interpolation inputs share the same column.
            if (derivId == densityIdA) {
                const real_t dQdnA = energyPerReaction * K * nB;
                jac->SetElement(row, offsetA + ir, dQdnA);
                contributes = true;
            }

            // Derivative with respect to the second reactant density.
            if (derivId == densityIdB) {
                const real_t dQdnB = energyPerReaction * K * nA;
                jac->SetElement(row, offsetB + ir, dQdnB);
                contributes = true;
            }

            // Additional dependence of K on electron density.
            if (derivId == id_ne && reaction.densityInput.kind == MolecularInputKind::ELECTRON) {
                // The interpolation coordinate is ne, in addition to any
                // explicit electron reactant already differentiated above.
                const real_t dKdne = reaction.rate->Eval_deriv_n(density, T);
                const real_t dQdne = energyPerReaction * nA * nB * dKdne;
                jac->SetElement(row, ir, dQdne);
                contributes = true;
            }

            // Dependence of K on electron temperature.
            if (temperature.kind == MolecularInputKind::ELECTRON && derivId == id_Te) {
                const real_t dKdTe = reaction.rate->Eval_deriv_T(density, T);
                const real_t dQdTe = energyPerReaction * nA * nB * dKdTe;
                jac->SetElement(row, ir, dQdTe);
                contributes = true;
            }

            // Dependence of K on a neutral temperature, through Wn and n.
            if (temperature.kind == MolecularInputKind::SPECIES &&
                (derivId == id_Wn || derivId == id_ions)) {
                const len_t species = ions->GetIonIndex(temperature.name1);
                const len_t col = ions->GetIndex(species, 0) * nr + ir;
                const real_t n = ni[col];
                // Match the zero-temperature branch in ResolveSpeciesTemperature.
                if (n <= 0)
                    continue;
                const real_t dKdT = reaction.rate->Eval_deriv_T(density, T);
                const real_t dQdT = energyPerReaction * nA * nB * dKdT;

                // Neutral temperature is Wn/(1.5*ec*n).
                if (derivId == id_Wn) {
                    const real_t dTdWn = 1 / (1.5 * Constants::ec * n);
                    const real_t dQdWn = dQdT * dTdWn;
                    jac->SetElement(row, species * nr + ir, dQdWn);
                } else {
                    const real_t dTdn = -T / n;
                    const real_t dQdn = dQdT * dTdn;
                    jac->SetElement(row, col, dQdn);
                }
                contributes = true;
            }
        }
    }
    return contributes;
}

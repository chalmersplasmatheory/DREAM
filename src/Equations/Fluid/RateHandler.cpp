/**
 * Implementation of the RateHandler class.
 * This class is responsible for managing and providing access to various rate coefficients
 * used in the fluid equations, including ADAS charge-state rates for atomic species,
 * as well as molecular reaction rates.
 */

#include "DREAM/Equations/Fluid/RateHandler.hpp"
#include "FVM/FVMException.hpp"
#include "DREAM/Equations/Fluid/RateData.hpp"
#include "DREAM/MolecularRateData.hpp"
#include "DREAM/MolecularRateInterpolator.hpp"
#include <algorithm>
#include <cmath>
#include "DREAM/Constants.hpp"


using namespace DREAM;

/**
* Constructor.
*/
RateHandler::RateHandler(FVM::Grid *grid, IonHandler *ions,
    ADAS *adas,
    real_t reionizScale, FVM::UnknownQuantityHandler *unknowns, bool reactionsEnabled,
    const std::vector<std::string>& enabledReactionNames,
    const std::unordered_map<std::string, UserMolecularRate>& userRates
) :  grid(grid), ions(ions), adas(adas), reionizScale(reionizScale),
     userRates(userRates), unknowns(unknowns) {



    this->unknowns  = unknowns;
    this->id_ions   = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    this->id_n_cold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    this->id_n_re   = unknowns->GetUnknownID(OptionConstants::UQTY_N_RE);
    this->id_T_cold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);

    // W_i and N_i only exist when the ion heat equation is included
    hasIonHeat = unknowns->HasUnknown(OptionConstants::UQTY_WI_ENER);
    if (hasIonHeat) {
        this->id_Ni = unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS);
        this->id_Wi = unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER);
    } else if (reactionsEnabled)
        throw FVM::FVMException(
            "RateHandler: Molecular reactions require the ion heat equation (T_i) to be included."
        );

    separateNeutrals = unknowns->HasUnknown(OptionConstants::UQTY_WN_ENER);
    if (separateNeutrals) {
        id_Wn = unknowns->GetUnknownID(OptionConstants::UQTY_WN_ENER);
        id_Nn = unknowns->GetUnknownID(OptionConstants::UQTY_NN_DENS);
    }


    if (reactionsEnabled){
        AddMolecularReactionRates(enabledReactionNames);
        AddMolecularLineRadiation();
    }

    const len_t n = molecularReactions.size() * grid->GetNr();
    this->Rate         = new real_t[n];
    this->Rate_N       = new real_t[n];
    this->Rate_T       = new real_t[n];
    this->ReactionRate = new real_t[n];
}

/**
 * Destructor.
 */
RateHandler::~RateHandler() {
    for (auto it = molecularReactions.begin(); it != molecularReactions.end(); it++) {
        delete it->rate;
    }
    for (auto& radiation : molecularLineRadiation)
        delete radiation.rate;

    delete [] this->Rate;
    delete [] this->Rate_N;
    delete [] this->Rate_T;
    delete [] this->ReactionRate;
}


/**
 * Get the atomic number used to look up the ADAS rates of ion species 'iz'.
 * Molecules are not in ADAS, they get their rates from MoleculeRateReaction.
 */
len_t RateHandler::GetADASZ(const len_t iz) const {
    if (ions->IsMolecule(iz))
        throw FVM::FVMException(
            "RateHandler: No ADAS rates available for molecule '%s'.",
            ions->GetName(iz).c_str()
        );
    return ions->GetZ(iz);
}

/**
 * Getters for the ADAS rates of ion species 'iz'.
 */
ADASRateInterpolator *RateHandler::GetACD(const len_t iz) const {return adas->GetACD(GetADASZ(iz));}

ADASRateInterpolator *RateHandler::GetSCD(const len_t iz) const {return adas->GetSCD(GetADASZ(iz));}

ADASRateInterpolator *RateHandler::GetPLT(const len_t iz) const {return adas->GetPLT(GetADASZ(iz));}

ADASRateInterpolator *RateHandler::GetPRB(const len_t iz) const {return adas->GetPRB(GetADASZ(iz));}

ADASRateInterpolator *RateHandler::GetCCD(const len_t iz) const {return adas->GetCCD(GetADASZ(iz));}

/**
 * Rebuild the rate coefficient K and the reaction rate R = K*nA*nB of every
 * molecular reaction. Called once per iteration by the solver, before the
 * equation terms are rebuilt, so that all terms (ion densities, ion heat,
 * neutral heat) use the same rates.
 */
void RateHandler::Rebuild() {
    const len_t nr = grid->GetNr();

    for (len_t k = 0; k < molecularReactions.size(); k++) {
        const MolecularReaction& reaction = molecularReactions[k];
        const real_t *nA = unknowns->GetUnknownData(reaction.reactantId[0]) + reaction.reactantOffset[0];
        const real_t *nB = unknowns->GetUnknownData(reaction.reactantId[1]) + reaction.reactantOffset[1];

        // Runaway ionization is rescaled by ions/reioniz_scale
        real_t scale = 1.0;
        if (reaction.process == MolecularReactionProcess::IONIZATION_RUNAWAY)
            scale = reionizScale;

        for (len_t ir = 0; ir < nr; ir++) {
            const real_t density = ResolveDensity(reaction.densityInput, ir);
            const real_t temperature = ResolveTemperature(reaction.temperatureInput, ir);

            Rate[k*nr + ir]   = scale * reaction.rate->Eval(density, temperature);
            Rate_N[k*nr + ir] = scale * reaction.rate->Eval_deriv_n(density, temperature);
            Rate_T[k*nr + ir] = scale * reaction.rate->Eval_deriv_T(density, temperature);
            ReactionRate[k*nr + ir] = Rate[k*nr + ir] * nA[ir] * nB[ir];
        }
    }
}

/**
 * Add factor * dR/dx to the Jacobian, where R = K(n, T)*nA*nB is the rate of
 * molecular reaction 'k' at radius 'ir', and x is the unknown 'derivId'.
 */
bool RateHandler::SetReactionRateJacobian(
    const len_t k, const len_t derivId, FVM::Matrix *jac,
    const len_t row, const len_t ir, const real_t factor
) {
    const len_t nr = grid->GetNr();
    const MolecularReaction& reaction = molecularReactions[k];
    const real_t nA = unknowns->GetUnknownData(reaction.reactantId[0])[reaction.reactantOffset[0] + ir];
    const real_t nB = unknowns->GetUnknownData(reaction.reactantId[1])[reaction.reactantOffset[1] + ir];
    bool contributes = false;
    bool added;

    if (derivId == reaction.reactantId[0]) {
        jac->SetElement(row, reaction.reactantOffset[0] + ir, factor * Rate[k*nr + ir] * nB);
        contributes = true;
    }
    if (derivId == reaction.reactantId[1]) {
        jac->SetElement(row, reaction.reactantOffset[1] + ir, factor * Rate[k*nr + ir] * nA);
        contributes = true;
    }

    if (derivId == id_n_cold && reaction.densityInput.kind == MolecularInputKind::ELECTRON) {
        jac->SetElement(row, ir, factor * Rate_N[k*nr + ir] * nA * nB);
        contributes = true;
    }

    const MolecularInput& Tin = reaction.temperatureInput;
    const real_t dRdT = factor * Rate_T[k*nr + ir] * nA * nB;
    if (Tin.kind == MolecularInputKind::ELECTRON && derivId == id_T_cold) {
        jac->SetElement(row, ir, dRdT);
        contributes = true;
    }
    else if (Tin.kind == MolecularInputKind::SPECIES) {
        added = SetSpeciesTemperatureJacobian(
            derivId, jac, row, ir, ions->GetIonIndex(Tin.name1), Tin.charge1, dRdT
        );
        if (added)
            contributes = true;
    }
    else if (Tin.kind == MolecularInputKind::RELATIVE) {
        added = SetSpeciesTemperatureJacobian(
            derivId, jac, row, ir, ions->GetIonIndex(Tin.name1), Tin.charge1, 0.5*dRdT
        );
        if (added)
            contributes = true;

        added = SetSpeciesTemperatureJacobian(
            derivId, jac, row, ir, ions->GetIonIndex(Tin.name2), Tin.charge2, 0.5*dRdT
        );
        if (added)
            contributes = true;
    }

    return contributes;
}

/**
 * Set the derivative of a quantity Q with respect to W and N of the
 * temperature T = W/(1.5*ec*N) of charge state Z0 of ion species.
 */
bool RateHandler::SetSpeciesTemperatureJacobian(
    const len_t derivId, FVM::Matrix *jac, const len_t row, const len_t ir,
    const len_t species, const int_t Z0, const real_t dQdT
) {
    const bool neutral = separateNeutrals && Z0 == 0;
    if (!neutral && !hasIonHeat)
        return false;

    const len_t id_W = neutral ? id_Wn : id_Wi;
    const len_t id_N = neutral ? id_Nn : id_Ni;
    if (derivId != id_W && derivId != id_N)
        return false;

    const len_t col = species * grid->GetNr() + ir;
    const real_t N = unknowns->GetUnknownData(id_N)[col];
    if (N <= 0)
        return false;

    if (derivId == id_W)
        jac->SetElement(row, col, dQdT / (1.5 * Constants::ec * N));
    else
        jac->SetElement(row, col, -dQdT * ResolveSpeciesTemperature(species, Z0, ir) / N);

    return true;
}

/*
 * Add molecular reaction rates for all defined molecular rate pairs.
 */
void RateHandler::AddMolecularReactionRates(const std::vector<std::string>& enabledReactionNames) {

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

            // We dont care if it has an electron or runaway electron
            if (std::string(r.name) == "e" || std::string(r.name) == "RE")
                continue;

            // If not all reactants exist (existence is handled by python)
            if (!ions->HasIon(r.name)) {
                allReactantsExist = false;
                break;
            }

            // Check if the charge state is valid for the reactant
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

            if (std::string(p.name) == "e" || std::string(p.name) == "RE")
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
        reaction.energyReleased = def.energyReleased;
        reaction.franckCondonEnergy = def.franckCondonEnergy;

        reaction.electronsConsumed = 0;
        for (len_t j = 0; j < def.nReactants; j++)
            if (std::string(def.reactants[j].name) == "e")
                reaction.electronsConsumed += (int_t)def.reactants[j].coefficient;
        for (len_t j = 0; j < def.nProducts; j++)
            if (std::string(def.products[j].name) == "e")
                reaction.electronsConsumed -= (int_t)def.products[j].coefficient;
        reaction.rate = GetMolecularRateByName(def.rateName);

        // Where the two reactant densities are stored
        if (def.nReactants != 2)
            throw FVM::FVMException(
                "RateHandler: Molecular reaction '%s' must have exactly two reactants.",
                def.rateName
            );
        for (len_t j = 0; j < 2; j++) {
            const MolecularReactionSpecies& r = def.reactants[j];
            reaction.reactantOffset[j] = 0;
            if (std::string(r.name) == "e")
                reaction.reactantId[j] = id_n_cold;
            else if (std::string(r.name) == "RE")
                reaction.reactantId[j] = id_n_re;
            else {
                reaction.reactantId[j] = id_ions;
                reaction.reactantOffset[j] = ions->GetIndex(ions->GetIonIndex(r.name), r.Z0) * grid->GetNr();
            }
        }

        molecularReactions.push_back(reaction);
    }

    printf(
        "RateHandler: Added %zu molecular reactions.\n",
        molecularReactions.size());

}

/**
 * Add the line radiation (molecular counterpart of PLT) of every molecular
 * species in the simulation that has a line-radiation rate in RateData.cpp.
 */
void RateHandler::AddMolecularLineRadiation() {
    for (len_t k = 0; k < molecularLineRadiationDefinitionCount; k++) {
        const auto& def = molecularLineRadiationDefinitions[k];
        if (!ions->HasIon(def.species))
            continue;

        const len_t iz = ions->GetIonIndex(def.species);
        if (def.Z0 < 0 || (len_t)def.Z0 > ions->GetZ(iz))
            throw FVM::FVMException(
                "RateHandler: Line radiation '%s' references invalid charge state Z0=%d for species '%s'.",
                def.rateName, def.Z0, def.species
            );

        molecularLineRadiation.push_back({iz, (len_t)def.Z0, GetMolecularRateByName(def.rateName)});
    }
}

/**
 * Helper method to get a molecular rate by name. A user-supplied
 * rate (from the settings) takes precedence over the pre-calculated
 * table in MolecularRateData.cpp.
 */
MolecularRateInterpolator *RateHandler::GetMolecularRateByName(
    const char *rateName
) const {
    auto it = userRates.find(rateName);
    if (it != userRates.end()) {
        const UserMolecularRate& r = it->second;
        molecular_rate rate = {
            rateName, r.logn.size(), r.logT.size(),
            r.logn.data(), r.logT.data(), r.coeff.data()
        };
        printf("RateHandler: Using user-supplied rate for '%s'.\n", rateName);
        return new MolecularRateInterpolator(&rate);
    }

    for (len_t i = 0; i < molecular_rate_n; i++) {
        // Pick up the correct rate in MolecularRateData.cpp by the correct name
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
    // Rate does not depend on temperature
    if (input.kind == MolecularInputKind::NONE) {
        return 1.0;
    }
    // Rate depends on electron temperature
    if (input.kind == MolecularInputKind::ELECTRON) {
        return unknowns->GetUnknownData(id_T_cold)[ir];
    }
    // Rate depends on a specific species temperature
    if (input.kind == MolecularInputKind::SPECIES) {
        return ResolveSpeciesTemperature(ions->GetIonIndex(input.name1), input.charge1, ir);
    }
    // Rate depends on the average temperature of two species
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

    if (separateNeutrals && charge == 0) {
        n = unknowns->GetUnknownData(id_Nn)[offset];
        W = unknowns->GetUnknownData(id_Wn)[offset];
    } else {
        n = unknowns->GetUnknownData(id_Ni)[offset];
        W = unknowns->GetUnknownData(id_Wi)[offset];
    }

    // Temperature is undefined for an absent population.
    return n > 0 ? W / (1.5 * Constants::ec * n) : 0;
}

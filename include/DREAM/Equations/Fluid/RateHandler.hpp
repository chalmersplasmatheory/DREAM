#ifndef _DREAM_RATE_HANDLER_HPP
#define _DREAM_RATE_HANDLER_HPP


#include <string>
#include <memory>
#include <unordered_map>
#include "FVM/config.h"
#include "DREAM/IonHandler.hpp"
#include "DREAM/ADAS.hpp"
#include "DREAM/Equations/Fluid/RateData.hpp"
#include "DREAM/MolecularRateInterpolator.hpp"
#include "DREAM/Settings/Settings.hpp"
#include <vector>

namespace DREAM {

//Struct to hold the molecular reaction rates for a given pair of species
//This is to connect the rate with the MolecularReaction struct.
struct MolecularReaction {
/// Name of the reaction
    const char *rateName;

    MolecularReactionProcess process;

    len_t nReactants;
    const MolecularReactionSpecies *reactants;

    len_t nProducts;
    const MolecularReactionSpecies *products;

    MolecularRateInterpolator *rate;
    /// The input for the temperature
    MolecularInput temperatureInput;
    /// The input for the density
    MolecularInput densityInput;
};

//Hold a user-supplied molecular rate table (given in the settings),
//which replaces the pre-calculated table with the same name.
struct UserMolecularRate {
    std::vector<real_t> logn;
    std::vector<real_t> logT;
    std::vector<real_t> coeff;
};


class RateHandler {
private:
    FVM::Grid *grid;
    IonHandler *ions;
    ADAS *adas;
    bool separateNeutrals;
    len_t id_Wn;

    len_t id_ions, id_n_cold, id_Ni, id_T_cold, id_Wi;

    std::vector<MolecularReaction> molecularReactions;
    std::unordered_map<std::string, UserMolecularRate> userRates;

    MolecularRateInterpolator *GetMolecularRateByName(const char *rateName) const;
    len_t GetADASZ(const len_t iz) const;

    void AddMolecularReactionRates( const std::vector<std::string>& enabledReactionNames);

    FVM::UnknownQuantityHandler *unknowns;

public:
    RateHandler(FVM::Grid *grid, IonHandler *ions, ADAS *adas, FVM::UnknownQuantityHandler *unknowns, bool reactionsEnabled, const std::vector<std::string>& enabledReactionNames, const std::unordered_map<std::string, UserMolecularRate>& userRates);
    ~RateHandler();

    const std::vector<MolecularReaction>& GetMolecularReactions() const {
    return molecularReactions;
    }

    ADASRateInterpolator *GetACD(const len_t iz) const;
    ADASRateInterpolator *GetSCD(const len_t iz) const;
    ADASRateInterpolator *GetPLT(const len_t iz) const;
    ADASRateInterpolator *GetPRB(const len_t iz) const;
    ADASRateInterpolator *GetCCD(const len_t iz) const;


    real_t ResolveDensity(const MolecularInput& input, const len_t ir);
    real_t ResolveTemperature(const MolecularInput& input,const len_t ir);
    real_t ResolveSpeciesTemperature(len_t species, int_t charge, len_t ir);
        
};

}

#endif /*_DREAM_RATE_HANDLER_HPP*/


#ifndef _DREAM_RATE_HANDLER_HPP
#define _DREAM_RATE_HANDLER_HPP


#include <string>
#include <memory>
#include <unordered_map>
#include "FVM/config.h"
#include "FVM/Matrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/ADAS.hpp"
#include "DREAM/AMJUEL.hpp"
#include "DREAM/Equations/Fluid/RateData.hpp"
#include "DREAM/MolecularRateInterpolator.hpp"
#include "DREAM/Settings/Settings.hpp"
#include <vector>

namespace DREAM {

    struct MolecularReaction {
        /// Name of the reaction
        const char *rateName;

        MolecularReactionProcess process;

        len_t nReactants;
        const MolecularReactionSpecies *reactants;

        len_t nProducts;
        const MolecularReactionSpecies *products;

        MolecularRateInterpolator *rate;
        MolecularInput temperatureInput;
        MolecularInput densityInput;
        real_t energyReleased;
        real_t franckCondonEnergy;
        int_t electronsConsumed;
    
        len_t reactantId[2];
        len_t reactantOffset[2];
        len_t reactantIz[2];
        len_t *productIz;     // nProducts elements
    };

    // Line radiation of one molecular charge state (molecular counterpart of PLT):
    struct MolecularLineRadiation {
        len_t iz;
        len_t Z0;
        MolecularRateInterpolator *rate;
    };

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
        bool hasIonHeat;        
        len_t id_Wn, id_Nn;
        len_t id_ions, id_n_cold, id_n_re, id_Ni, id_T_cold, id_Wi;
        real_t reionizScale;

        real_t *Rate = nullptr;
        real_t *Rate_N = nullptr;
        real_t *Rate_T = nullptr;
        real_t *ReactionRate = nullptr;

        std::vector<MolecularReaction> molecularReactions;
        std::unordered_map<std::string, UserMolecularRate> userRates;
        std::vector<MolecularLineRadiation> molecularLineRadiation;

        MolecularRateInterpolator *GetMolecularRateByName(const char *rateName) const;
        len_t GetADASZ(const len_t iz) const;

        void AddMolecularReactionRates(
            const std::vector<std::string>& enabledReactionNames,
            const int_t *reactantSpecies, const int_t *productSpecies,
            const int_t *temperatureSpecies
        );
        void AddMolecularLineRadiation(const int_t *lineRadiationSpecies);

        FVM::UnknownQuantityHandler *unknowns;

    public:
        RateHandler(
            FVM::Grid *grid, IonHandler *ions, ADAS *adas, real_t reionizScale,
            FVM::UnknownQuantityHandler *unknowns, bool reactionsEnabled,
            const std::vector<std::string>& enabledReactionNames,
            const std::unordered_map<std::string, UserMolecularRate>& userRates,
            const int_t *reactantSpecies = nullptr,
            const int_t *productSpecies = nullptr,
            const int_t *temperatureSpecies = nullptr,
            const int_t *lineRadiationSpecies = nullptr
        );
        ~RateHandler();

        void Rebuild();

        const std::vector<MolecularReaction>& GetMolecularReactions() const {return molecularReactions;}


        const real_t *GetRateCoefficient(const len_t k) const { return Rate + k*grid->GetNr(); }
        const real_t *GetReactionRate(const len_t k) const { return ReactionRate + k*grid->GetNr(); }


        bool SetReactionRateJacobian(
            const len_t k, const len_t derivId, FVM::Matrix *jac,
            const len_t row, const len_t ir, const real_t factor
        );
        bool SetSpeciesTemperatureJacobian(
            const len_t derivId, FVM::Matrix *jac, const len_t row, const len_t ir,
            const len_t species, const int_t Z0, const real_t dQdT
        );

        ADASRateInterpolator *GetACD(const len_t iz) const;
        ADASRateInterpolator *GetSCD(const len_t iz) const;
        ADASRateInterpolator *GetPLT(const len_t iz) const;
        ADASRateInterpolator *GetPRB(const len_t iz) const;
        ADASRateInterpolator *GetCCD(const len_t iz) const;

        real_t GetReionizScale() const { return reionizScale; }

        const std::vector<MolecularLineRadiation>& GetMolecularLineRadiation() const {
            return molecularLineRadiation;
        }

        // ADAS ionization rate by atomic number, for a species that is not part of
        // the plasma (e.g. the neutral beam)
        ADASRateInterpolator *GetSCDByZ(const len_t Z) const { return adas->GetSCD(Z); }


        real_t ResolveDensity(const MolecularInput& input, const len_t ir);
        real_t ResolveTemperature(const MolecularInput& input, const len_t ir);
        real_t ResolveSpeciesTemperature(len_t species, int_t charge, len_t ir);
    };
}

#endif/*_DREAM_RATE_HANDLER_HPP*/

#ifndef _DREAM_EQUATIONS_RUNAWAY_SOURCE_TERM_HANDLER_HPP
#define _DREAM_EQUATIONS_RUNAWAY_SOURCE_TERM_HANDLER_HPP

#include <string>
#include "DREAM/Equations/Fluid/AvalancheGrowthTerm.hpp"
#include "DREAM/Equations/Fluid/ComptonRateTerm.hpp"
#include "DREAM/Equations/Fluid/DreicerRateTerm.hpp"
#include "DREAM/Equations/Fluid/ExternalAvalancheTerm.hpp"
#include "DREAM/Equations/Fluid/HottailRateTerm.hpp"
#include "DREAM/Equations/Fluid/TritiumRateTerm.hpp"
#include "DREAM/Equations/Kinetic/AvalancheSourceRP.hpp"
#include "DREAM/Equations/Kinetic/TritiumSource.hpp"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Equation/Operator.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"

namespace DREAM {
    class RunawaySourceTermHandler {
    protected:
        FVM::EquationTerm *avalanche=nullptr;       // Avalanche source applied to n_re
        FVM::EquationTerm *avalanche_neg=nullptr;   // Avalanche source applied to n_re_neg (same p|| direction as, but subtracted from, n_re)
        FVM::EquationTerm *avalanche_negpos=nullptr;// Avalanche source applied to n_re_neg (oppositive p|| direction as n_re; positive contribution to RE generation)
        ComptonRateTerm *compton=nullptr;
        DreicerRateTerm *dreicer=nullptr;
        HottailRateTerm *hottail=nullptr;
        // There can be multiple tritium species in the simulation...
        std::vector<FVM::EquationTerm*> tritium;

        // Description of equation
        std::string description;        

        void applyToAll(std::function<void(FVM::EquationTerm*)>);
        void applyToAll(const std::function<void(FVM::EquationTerm*)>) const;

    public:
        RunawaySourceTermHandler();
        ~RunawaySourceTermHandler();

        void AddSourceTerm(const std::string& desc, AvalancheGrowthTerm *t) { this->description += desc; this->avalanche = t; }
        void AddSourceTerm(const std::string& desc, AvalancheSourceRP *t) { this->description += desc; this->avalanche = t; }
        void AddAvalancheNreNeg(AvalancheSourceRP *t) { this->avalanche_neg = t; }
        void AddAvalancheNreNegPos(AvalancheSourceRP *t) { this->avalanche_negpos = t; }
        void AddSourceTerm(const std::string& desc, ExternalAvalancheTerm *t) { this->description += desc; this->avalanche = t; }
        void AddSourceTerm(const std::string& desc, ComptonRateTerm *t) { this->description += desc; this->compton = t; }
        void AddSourceTerm(const std::string& desc, DreicerRateTerm *t) { this->description += desc; this->dreicer = t; }
        void AddSourceTerm(const std::string& desc, HottailRateTerm *t) { this->description += desc; this->hottail = t; }
        void AddSourceTerm(const std::string& desc, TritiumRateTerm *t) {
            if (this->tritium.empty())
                this->description += desc;

            this->tritium.push_back(t);
        }
        void AddSourceTerm(const std::string& desc, TritiumSource *t){
            if (this->tritium.empty())
                this->description += desc;

            this->tritium.push_back(t);
        }

        void AddToOperators(FVM::Operator*, FVM::Operator*, FVM::Operator*, FVM::Operator *nRE_neg=nullptr);

        const std::string& GetDescription() const { return description; }
    };
}

#endif/*_DREAM_EQUATIONS_RUNAWAY_SOURCE_TERM_HANDLER_HPP*/

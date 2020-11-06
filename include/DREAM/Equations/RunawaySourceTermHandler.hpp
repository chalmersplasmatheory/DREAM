#ifndef _DREAM_EQUATIONS_RUNAWAY_SOURCE_TERM_HANDLER_HPP
#define _DREAM_EQUATIONS_RUNAWAY_SOURCE_TERM_HANDLER_HPP

#include <string>
#include "DREAM/Equations/Fluid/AvalancheGrowthTerm.hpp"
#include "DREAM/Equations/Fluid/ComptonRateTerm.hpp"
#include "DREAM/Equations/Fluid/DreicerRateTerm.hpp"
#include "DREAM/Equations/Fluid/TritiumRateTerm.hpp"
#include "DREAM/Equations/Kinetic/AvalancheSourceRP.hpp"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Equation/Operator.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"

namespace DREAM {
    class RunawaySourceTermHandler {
    protected:
        FVM::EquationTerm *avalanche=nullptr;
        ComptonRateTerm *compton=nullptr;
        DreicerRateTerm *dreicer=nullptr;
        // There can be multiple tritium species in the simulation...
        std::vector<TritiumRateTerm*> tritium;

        // Description of equation
        std::string description;

        void applyToAll(std::function<void(FVM::EquationTerm*)>);
        void applyToAll(const std::function<void(FVM::EquationTerm*)>) const;

    public:
        RunawaySourceTermHandler();
        ~RunawaySourceTermHandler();

        void AddSourceTerm(const std::string& desc, AvalancheGrowthTerm *t) { this->description += desc; this->avalanche = t; }
        void AddSourceTerm(const std::string& desc, AvalancheSourceRP *t) { this->description += desc; this->avalanche = t; }
        void AddSourceTerm(const std::string& desc, ComptonRateTerm *t) { this->description += desc; this->compton = t; }
        void AddSourceTerm(const std::string& desc, DreicerRateTerm *t) { this->description += desc; this->dreicer = t; }
        void AddSourceTerm(const std::string& desc, TritiumRateTerm *t) {
            if (this->tritium.empty())
                this->description += desc;

            this->tritium.push_back(t);
        }

        void AddToOperators(FVM::Operator*, FVM::Operator*, FVM::Operator*);

        const std::string& GetDescription() const { return description; }
    };
}

#endif/*_DREAM_EQUATIONS_RUNAWAY_SOURCE_TERM_HANDLER_HPP*/

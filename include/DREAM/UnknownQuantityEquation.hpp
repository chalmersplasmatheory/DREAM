#ifndef _DREAM_UNKNOWN_QUANTITY_EQUATION_HPP
#define _DREAM_UNKNOWN_QUANTITY_EQUATION_HPP

#include <map>
#include <string>
#include "FVM/BlockMatrix.hpp"
#include "FVM/Equation/Equation.hpp"
#include "FVM/Equation/PredeterminedParameter.hpp"
#include "FVM/UnknownQuantity.hpp"


namespace DREAM {
    class UnknownQuantityEquation {
    private:
        // Pointer to associated data container (owned by 'UnknownQuantityHandler')
        FVM::UnknownQuantity *uqty;

        // List of equations associated with this quantity
        std::map<len_t, FVM::Equation*> equations;

        std::string description;

    public:
        UnknownQuantityEquation(FVM::UnknownQuantity *uqty, const std::string& desc="")
            : uqty(uqty), description(desc) { }
        ~UnknownQuantityEquation();

        void Evaluate(const len_t, real_t*, FVM::UnknownQuantityHandler*);

        const std::string& GetDescription() const { return this->description; }
        const std::map<len_t, FVM::Equation*>& GetEquations() const { return this->equations; }
        const FVM::Equation *GetEquation(const len_t i) const { return this->equations.at(i); }

        len_t NumberOfElements() const { return this->uqty->NumberOfElements(); }
        len_t NumberOfNonZeros();
        len_t NumberOfNonZeros_jac();

        bool IsEvaluable();
        FVM::PredeterminedParameter *GetPredetermined();
        bool IsPredetermined();
        void RebuildEquations(const real_t, const real_t, FVM::UnknownQuantityHandler*);

        void SetDescription(const std::string& desc) { this->description = desc; }

        void SetEquation(const len_t blockcol, FVM::Equation *eqn) {
            this->equations[blockcol] = eqn;
        }

        // Evaluation methods
        void SetVectorElements(
            real_t*, FVM::UnknownQuantityHandler*, FVM::BlockMatrix*,
            std::map<len_t, len_t>&
        );
    };
}

#endif/*_DREAM_UNKNOWN_QUANTITY_EQUATION_HPP*/

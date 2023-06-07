#ifndef _DREAM_UNKNOWN_QUANTITY_EQUATION_HPP
#define _DREAM_UNKNOWN_QUANTITY_EQUATION_HPP

#include <map>
#include <string>
#include "FVM/BlockMatrix.hpp"
#include "FVM/Equation/Operator.hpp"
#include "FVM/Equation/PredeterminedParameter.hpp"
#include "FVM/UnknownQuantity.hpp"


namespace DREAM {
    class UnknownQuantityEquation {
    private:
        // ID of associated unknown quantity
        len_t uqtyId;

        // Pointer to associated data container (owned by 'UnknownQuantityHandler')
        FVM::UnknownQuantity *uqty;

        // List of equations associated with this quantity
        std::map<len_t, FVM::Operator*> equations;

		// Flag indicating whether this quantity is supposed to be
		// evaluated using the external iterator
		bool solved_externally = false;

        std::string description;

    public:
        UnknownQuantityEquation(len_t uqtyId, FVM::UnknownQuantity *uqty, const std::string& desc="")
            : uqtyId(uqtyId), uqty(uqty), description(desc) { }
        ~UnknownQuantityEquation();

        void Evaluate(const len_t, real_t*, FVM::UnknownQuantityHandler*);

        const std::string& GetDescription() const { return this->description; }
        const std::map<len_t, FVM::Operator*>& GetOperators() const { return this->equations; }
        const FVM::Operator *GetOperator(const len_t i) const { return this->equations.at(i); }
        FVM::Operator *GetOperatorUnsafe(const len_t i) { return this->equations.at(i); }
        FVM::UnknownQuantity *GetUnknown() { return this->uqty; }

        len_t NumberOfElements() const { return this->uqty->NumberOfElements(); }
        len_t NumberOfNonZeros();
        len_t NumberOfNonZeros_jac();

        bool IsEvaluable();
        FVM::PredeterminedParameter *GetPredetermined();
        bool IsPredetermined();
		bool IsSolvedExternally() { return this->solved_externally; }
        void RebuildEquations(const real_t, const real_t, FVM::UnknownQuantityHandler*);

        void SetDescription(const std::string& desc) { this->description = desc; }
		void SetExternallySolved(bool s) { this->solved_externally = s; }

        void SetOperator(const len_t blockcol, FVM::Operator *eqn) {
            this->equations[blockcol] = eqn;
        }

        // Evaluation methods
        void SetVectorElements(
            real_t*, FVM::UnknownQuantityHandler*
        );
    };
}

#endif/*_DREAM_UNKNOWN_QUANTITY_EQUATION_HPP*/

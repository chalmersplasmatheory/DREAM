#ifndef _DREAM_UNKNOWN_QUANTITY_EQUATION_HPP
#define _DREAM_UNKNOWN_QUANTITY_EQUATION_HPP

#include <map>
#include <string>
#include "FVM/BlockMatrix.hpp"
#include "FVM/Equation/Operator.hpp"
#include "FVM/Equation/PredeterminedParameter.hpp"
#include "FVM/UnknownQuantity.hpp"
#include "Trigger/EquationTriggerCondition.hpp"


namespace DREAM {
    class UnknownQuantityEquation {
    private:
        // ID of associated unknown quantity
        len_t uqtyId;

        // Pointer to associated data container (owned by 'UnknownQuantityHandler')
        FVM::UnknownQuantity *uqty;

        // List of equations associated with this quantity
        std::map<len_t, FVM::Operator*> equations;
		std::map<len_t, FVM::Operator*> equations_alt;

		// Trigger condition
		EquationTriggerCondition *condition;
		real_t *eqn_cache = nullptr;
		len_t nElements = 0;

		// Flag indicating whether this quantity is supposed to be
		// evaluated using the external iterator
		bool solved_externally = false;

        std::string description;
		std::string description_alt;

    public:
        UnknownQuantityEquation(
			len_t uqtyId, FVM::UnknownQuantity *uqty, const std::string& desc="",
			EquationTriggerCondition *condition=nullptr
		);
        ~UnknownQuantityEquation();

        void Evaluate(const len_t, real_t*, FVM::UnknownQuantityHandler*);

		const bool *GetAlternativeEquationMask() { return this->condition->GetTriggerMask(); }
        std::string GetDescription() const { return this->description; }
        const std::string& GetDescriptionAlt() const { return this->description_alt; }
        const std::map<len_t, FVM::Operator*>& GetOperators(bool alt=false) const {
			if (alt) return this->equations_alt;
			else return this->equations;
		}
		std::string GetOperatorNames(const std::map<len_t, FVM::Operator*> *eqns=nullptr) const;
		std::string GetAltOperatorNames() const;
        const FVM::Operator *GetOperator(const len_t i) const { return this->equations.at(i); }
		const FVM::Operator *GetOperatorAlt(const len_t i) const { return this->equations_alt.at(i); }
        FVM::Operator *GetOperatorUnsafe(const len_t i) { return this->equations.at(i); }
        FVM::UnknownQuantity *GetUnknown() { return this->uqty; }
		EquationTriggerCondition *GetTriggerCondition() { return this->condition; }
		bool HasOperatorAt(const len_t i) const { return (this->equations.find(i) != this->equations.end()); }

		void CheckTriggerCondition(const real_t, FVM::UnknownQuantityHandler*);
		bool HasAlternativeEquation() const { return (this->condition != nullptr); }
		bool HasTransientTerm() const;
		void SaveTriggerState();

        len_t NumberOfElements() const { return this->uqty->NumberOfElements(); }
        len_t NumberOfNonZeros();
        len_t NumberOfNonZeros_jac();

        bool IsEvaluable();
        FVM::PredeterminedParameter *GetPredetermined();
        bool IsPredetermined();
		bool IsSolvedExternally() { return this->solved_externally; }
        void RebuildEquations(const real_t, const real_t, FVM::UnknownQuantityHandler*);

        void SetDescription(const std::string& desc) { this->description = desc; }
		void SetDescriptionAlt(const std::string& desc) { this->description_alt = desc; }
		void SetExternallySolved(bool s) { this->solved_externally = s; }
		void SetTriggerCondition(EquationTriggerCondition *c);

        void SetOperator(const len_t blockcol, FVM::Operator *eqn) {
            this->equations[blockcol] = eqn;
        }
		void SetOperatorAlt(const len_t blockcol, FVM::Operator *eqn) {
			this->equations_alt[blockcol] = eqn;
		}

        // Evaluation methods
        void SetVectorElements(
            real_t*, FVM::UnknownQuantityHandler*
        );

		// Write trigger condition diagnostic information
		void SaveTriggerConditionDiagnostics(SFile*, const std::string&);
    };
}

#endif/*_DREAM_UNKNOWN_QUANTITY_EQUATION_HPP*/

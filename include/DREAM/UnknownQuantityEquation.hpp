#ifndef _DREAM_UNKNOWN_QUANTITY_EQUATION_HPP
#define _DREAM_UNKNOWN_QUANTITY_EQUATION_HPP

#include <map>
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

    public:
        UnknownQuantityEquation(FVM::UnknownQuantity *uqty) { this->uqty = uqty; }
        ~UnknownQuantityEquation();

        const std::map<len_t, FVM::Equation*>& GetEquations() const { return this->equations; }

        len_t NumberOfElements() const { return this->uqty->NumberOfElements(); }
        len_t NumberOfNonZeros();
        len_t NumberOfNonZeros_jac();

        FVM::PredeterminedParameter *GetPredetermined();
        bool IsPredetermined();
        void RebuildEquations(const real_t, const real_t, FVM::UnknownQuantityHandler*);

        void SetEquation(const len_t blockcol, FVM::Equation *eqn) {
            this->equations[blockcol] = eqn;
        }
    };
}

#endif/*_DREAM_UNKNOWN_QUANTITY_EQUATION_HPP*/

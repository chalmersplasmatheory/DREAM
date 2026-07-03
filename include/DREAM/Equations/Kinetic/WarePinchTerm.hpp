#ifndef _DREAM_EQUATIONS_WAREPINCH_TERM_HPP
#define _DREAM_EQUATIONS_WAREPINCH_TERM_HPP

#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    class WarePinchTerm : public FVM::AdvectionTerm {
    private:
        enum OptionConstants::momentumgrid_type mgtype;
        len_t id_E_field;
        real_t *dFrMatrix = nullptr;

    public:
        WarePinchTerm(FVM::Grid*, enum OptionConstants::momentumgrid_type,FVM::UnknownQuantityHandler*);
        ~WarePinchTerm();
    
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
        virtual void SetPartialAdvectionTerm(const len_t, const len_t) override;
    };
}

#endif/*_DREAM_EQUATIONS_WAREPINCH_TERM_HPP*/

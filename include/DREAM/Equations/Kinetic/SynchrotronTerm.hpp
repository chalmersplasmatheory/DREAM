#ifndef _DREAM_EQUATIONS_KINETIC_SYNCHROTRON_TERM_HPP
#define _DREAM_EQUATIONS_KINETIC_SYNCHROTRON_TERM_HPP


namespace DREAM { class SynchrotronTerm; }

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include "FVM/config.h"
#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"


namespace DREAM {
    class SynchrotronTerm
        : public FVM::AdvectionTerm {
    private:
        enum OptionConstants::momentumgrid_type gridtype;
        const real_t constPrefactor = Constants::ec * Constants::ec * Constants::ec * Constants::ec / (
            6*M_PI * Constants::eps0 * Constants::me * Constants::me * Constants::me * 
            Constants::c * Constants::c * Constants::c);
    public:
        SynchrotronTerm(FVM::Grid*, enum OptionConstants::momentumgrid_type);
        
        
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_DREAM_EQUATIONS_KINETIC_SYNCHROTRON_TERM_HPP*/



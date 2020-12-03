#ifndef _DREAM_EQUATIONS_KINETIC_SLOWING_DOWN_TERM_HPP
#define _DREAM_EQUATIONS_KINETIC_SLOWING_DOWN_TERM_HPP


#include "FVM/config.h"
#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include "DREAM/Equations/SlowingDownFrequency.hpp"
#include "FVM/UnknownQuantityHandler.hpp"


namespace DREAM {
    class SlowingDownTerm
        : public FVM::AdvectionTerm {
    private:
        enum OptionConstants::momentumgrid_type gridtype;
        SlowingDownFrequency *nuS;
        virtual void SetPartialAdvectionTerm(len_t derivId, len_t nMultiples) override;
    public:

        SlowingDownTerm(FVM::Grid*,CollisionQuantityHandler*,enum OptionConstants::momentumgrid_type, 
                        FVM::UnknownQuantityHandler*, bool withKineticIonJacobian);
        
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_DREAM_EQUATIONS_KINETIC_SLOWING_DOWN_TERM_HPP*/



#ifndef _DREAM_EQUATIONS_KINETIC_BRAAMS_KARNEY_ADVECTION_HPP
#define _DREAM_EQUATIONS_KINETIC_BRAAMS_KARNEY_ADVECTION_HPP


#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include "FVM/config.h"
#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"


namespace DREAM {
    class BraamsKarneyAdvection
        : public FVM::AdvectionTerm {
    private:
        len_t id_f, id_pi_0, id_pi_1;
        virtual void SetPartialAdvectionTerm(len_t derivId, len_t nMultiples) override;

    public:
        BraamsKarneyAdvection(FVM::Grid*, FVM::UnknownQuantityHandler*, len_t id_f, len_t id_pi_0, len_t id_pi_1);
        
        
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_DREAM_EQUATIONS_KINETIC_BRAAMS_KARNEY_ADVECTION_HPP*/

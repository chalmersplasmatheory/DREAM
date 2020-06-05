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
        len_t id_ncold, id_ni, nzs;
        void GetDF(len_t derivId, real_t **&df1, real_t **&df2, len_t lenDeriv);
    public:
        SlowingDownTerm(FVM::Grid*,CollisionQuantityHandler*,enum OptionConstants::momentumgrid_type);
        
        virtual void SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;

        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_DREAM_EQUATIONS_KINETIC_SLOWING_DOWN_TERM_HPP*/



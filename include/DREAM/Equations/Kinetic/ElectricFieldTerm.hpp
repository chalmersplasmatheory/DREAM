#ifndef _DREAM_EQUATIONS_KINETIC_ELECTRIC_FIELD_TERM_HPP
#define _DREAM_EQUATIONS_KINETIC_ELECTRIC_FIELD_TERM_HPP


#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include "FVM/config.h"
#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"


namespace DREAM {
    class ElectricFieldTerm
        : public FVM::AdvectionTerm {
    private:
        len_t id_Eterm;
        bool gridtypePXI;
        bool gridtypePPARPPERP;
        virtual void SetPartialAdvectionTerm(len_t derivId, len_t nMultiples) override;

    public:
        ElectricFieldTerm(FVM::Grid*,FVM::UnknownQuantityHandler*, enum OptionConstants::momentumgrid_type);
        
        
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_DREAM_EQUATIONS_KINETIC_ELECTRIC_FIELD_TERM_HPP*/



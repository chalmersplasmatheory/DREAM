#ifndef _DREAM_EQUATIONS_KINETIC_ISOTROPIC_SOURCE_P_XI_HPP
#define _DREAM_EQUATIONS_KINETIC_ISOTROPIC_SOURCE_P_XI_HPP

#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include "DREAM/Equations/SlowingDownFrequency.hpp"
#include "FVM/Equation/BoundaryConditions/PInternalBoundaryCondition.hpp"

namespace DREAM {
    class BCIsotropicSourcePXi : public FVM::BC::PInternalBoundaryCondition {
    private:
        SlowingDownFrequency *slowingDownFreq;
        len_t id_f;

    public:
        BCIsotropicSourcePXi(FVM::Grid*, CollisionQuantityHandler*, len_t);

        bool Rebuild(const real_t, FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_DREAM_EQUATIONS_KINETIC_ISOTROPIC_SOURCE_P_XI_HPP*/

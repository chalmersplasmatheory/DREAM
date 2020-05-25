#ifndef _DREAM_EQUATIONS_KINETIC_ISOTROPIC_SOURCE_HPP
#define _DREAM_EQUATIONS_KINETIC_ISOTROPIC_SOURCE_HPP

#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include "DREAM/Equations/SlowingDownFrequency.hpp"
#include "FVM/Equation/BoundaryConditions/PInternalBoundaryCondition.hpp"

namespace DREAM {
    class BCIsotropicSource : public FVM::BC::PInternalBoundaryCondition {
    private:
        SlowingDownFrequency *slowingDownFreq;
    public:
        BCIsotropicSource(FVM::Grid*, CollisionQuantityHandler*);

        bool Rebuild(const real_t, FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_DREAM_EQUATIONS_KINETIC_ISOTROPIC_SOURCE_HPP*/

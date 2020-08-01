#ifndef _DREAM_EQUATIONS_PARTICLE_SOURCE_TERM_HPP
#define _DREAM_EQUATIONS_PARTICLE_SOURCE_TERM_HPP

#include "DREAM/Equations/Kinetic/FluidKineticSourceTerm.hpp"

namespace DREAM {
    class ParticleSourceTerm
        : public FluidKineticSourceTerm {

    public:
        enum ParticleSourceShape{
            PARTICLE_SOURCE_SHAPE_MAXWELLIAN = 1, // particle source takes the shape of the instantaneous Maxwellian
            PARTICLE_SOURCE_SHAPE_DELTA = 2       // particle source is zero at all points except p=0
        };

    private:
        ParticleSourceShape particleSourceShape;
        len_t id_Tcold;
    protected:
        virtual real_t GetSourceFunction(len_t ir, len_t i, len_t j) override;
        virtual real_t GetSourceFunctionJacobian(len_t ir, len_t i, len_t j, const len_t derivId) override;
    public:
        ParticleSourceTerm(FVM::Grid*, FVM::UnknownQuantityHandler*, ParticleSourceShape = PARTICLE_SOURCE_SHAPE_MAXWELLIAN);
    };
}

#endif/*_DREAM_EQUATIONS_PARTICLE_SOURCE_TERM_HPP*/



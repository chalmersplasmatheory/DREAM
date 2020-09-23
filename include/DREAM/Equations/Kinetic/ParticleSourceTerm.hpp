#ifndef _DREAM_EQUATIONS_PARTICLE_SOURCE_TERM_HPP
#define _DREAM_EQUATIONS_PARTICLE_SOURCE_TERM_HPP

#include "DREAM/Equations/FluidSourceTerm.hpp"

namespace DREAM {
    class ParticleSourceTerm
        : public FluidSourceTerm {

    public:
        enum ParticleSourceShape{
            PARTICLE_SOURCE_SHAPE_MAXWELLIAN = 1, // particle source takes the shape of a Maxwellian
            PARTICLE_SOURCE_SHAPE_DELTA = 2       // particle source is zero at all points except p=0
        };

    private:
        ParticleSourceShape particleSourceShape;
        len_t id_Tcold;
        real_t nRef; // reference density used in Maxwellian source
    protected:
        virtual real_t GetSourceFunction(len_t ir, len_t i, len_t j) override;
        virtual real_t GetSourceFunctionJacobian(len_t ir, len_t i, len_t j, const len_t derivId) override;
    public:
        ParticleSourceTerm(FVM::Grid*, FVM::UnknownQuantityHandler*, ParticleSourceShape = PARTICLE_SOURCE_SHAPE_MAXWELLIAN);

        // Rebuild: set normalization factor of Maxwellian source
        virtual void Rebuild(const real_t, const real_t dt, FVM::UnknownQuantityHandler*) override 
            {this->nRef = 1/dt;}

    };
}

#endif/*_DREAM_EQUATIONS_PARTICLE_SOURCE_TERM_HPP*/



#ifndef _DREAM_EQUATION_FLUID_ION_SPECIES_TRANSIENT_TERM_HPP
#define _DREAM_EQUATION_FLUID_ION_SPECIES_TRANSIENT_TERM_HPP

/**
 * Implementation of an equation term representing the
 * (backward euler) time derivative on an ion species 
 * with index iz. (which can be scaled using scalefactor)
 */

#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM{
    class IonSpeciesTransientTerm : public FVM::EquationTerm {
    private:
        const len_t unknownId;
        len_t iz;
        real_t scaleFactor;
        real_t dt;
        real_t *xPrev;
    public:
        IonSpeciesTransientTerm(FVM::Grid *g, len_t iz, const len_t id, real_t scaleFactor=1.0) 
            : FVM::EquationTerm(g), unknownId(id), iz(iz), scaleFactor(scaleFactor){}
        virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return 1; }
        virtual void Rebuild(const real_t, const real_t dt, FVM::UnknownQuantityHandler *u) {
            this->dt = dt;
            this->xPrev = u->GetUnknownDataPrevious(this->unknownId);
        }

        virtual void SetJacobianBlock(const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t*) {
            if(uqtyId==derivId)
                SetMatrixElements(jac, nullptr);
        }
        virtual void SetMatrixElements(FVM::Matrix *mat, real_t *rhs) override {
            for(len_t ir=0; ir<nr; ir++)
                mat->SetElement(iz*nr+ir,iz*nr+ir,scaleFactor/dt);
            if(rhs!=nullptr)
                for(len_t ir=0; ir<nr; ir++)
                    rhs[iz*nr+ir] -= scaleFactor/dt;
        }
        virtual void SetVectorElements(real_t *vec, const real_t *xi){
            for(len_t ir=0; ir<nr; ir++)
                vec[iz*nr+ir] += scaleFactor*(xi[iz*nr+ir]-xPrev[iz*nr+ir])/dt;
        }
    };
}

#endif/*_DREAM_EQUATION_FLUID_ION_SPECIES_TRANSIENT_TERM_HPP*/

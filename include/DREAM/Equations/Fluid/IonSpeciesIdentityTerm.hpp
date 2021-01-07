#ifndef _DREAM_EQUATION_FLUID_ION_SPECIES_IDENTITY_TERM_HPP
#define _DREAM_EQUATION_FLUID_ION_SPECIES_IDENTITY_TERM_HPP

#include "FVM/Equation/EvaluableEquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

/**
 * Implementation of an equation term representing an
 * identity operator acting on an ion species  with index iz
 * (which can be scaled using scalefactor) 
 */

namespace DREAM{
    class IonSpeciesIdentityTerm : public FVM::EvaluableEquationTerm {
    private:
        len_t iz;
        real_t scaleFactor;
    public:
        IonSpeciesIdentityTerm(FVM::Grid *g, len_t iz, real_t scaleFactor=1.0) 
            : FVM::EvaluableEquationTerm(g), iz(iz), scaleFactor(scaleFactor){}
        virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return 1; }
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override {}

        virtual void SetJacobianBlock(const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t*) override {
            if(uqtyId==derivId)
                SetMatrixElements(jac, nullptr);
        }
        virtual void SetMatrixElements(FVM::Matrix *mat, real_t*) override {
            for(len_t ir=0; ir<nr; ir++)
                mat->SetElement(iz*nr+ir,iz*nr+ir,scaleFactor);
        }
        virtual void SetVectorElements(real_t *vec, const real_t *xi) override {
            for(len_t ir=0; ir<nr; ir++)
                vec[iz*nr+ir] += scaleFactor*xi[iz*nr+ir];
        }
        virtual void EvaluableTransform(real_t *vec) override {
            for (len_t ir=0; ir<nr; ir++)
                vec[iz*nr+ir] = -vec[iz*nr+ir] / scaleFactor;
        }

    };
}

#endif/*_DREAM_EQUATION_FLUID_ION_SPECIES_IDENTITY_TERM_HPP*/

#ifndef _DREAM_EQUATION_FLUID_KINETIC_EQUATION_TERM_INTEGRATED_OVER_MOMENTUM_HPP
#define _DREAM_EQUATION_FLUID_KINETIC_EQUATION_TERM_INTEGRATED_OVER_MOMENTUM_HPP

#include "FVM/Equation/Operator.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    class KineticEquationTermIntegratedOverMomentum : public FVM::EquationTerm {
    private:
        FVM::Grid *kineticGrid;
        FVM::Operator *kineticOperator;
        const len_t id_f;
        FVM::UnknownQuantityHandler *unknowns;
        const PetscScalar scaleFactor;

        // Quantities for temporary storage of values used in the integrals:
        FVM::Matrix *kineticMatrix;      // matrix with size equal to the (id_f,id_f) block matrix
        FVM::Matrix *integrationMatrix;  // matrix with size equal to the (id_x,id_f) block matrix that contains integration weights
        PetscScalar *kineticVector;      // vector of size id_f
        PetscScalar *fluidVector;        // vector of size nr 

        PetscInt *idxFluid;
        len_t NCells;

        void zeroVecs();
        void allocateKineticStorage();
        void deallocateKineticStorage();
    public:
        KineticEquationTermIntegratedOverMomentum(FVM::Grid*, FVM::Grid*, FVM::Operator*, const len_t, FVM::UnknownQuantityHandler *u, real_t scaleFactor = 1.0);
        ~KineticEquationTermIntegratedOverMomentum();

        virtual len_t GetNumberOfNonZerosPerRow() const override {return 1;};
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override;

        virtual void Rebuild(const real_t t, const real_t dt, FVM::UnknownQuantityHandler *u) override 
            {kineticOperator->RebuildTerms(t,dt,u);}
        virtual bool GridRebuilt() override;
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
        virtual void SetJacobianBlock(const len_t uqtyId, const len_t derivId, FVM::Matrix*, const real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}
#endif/*_DREAM_EQUATION_FLUID_KINETIC_EQUATION_TERM_INTEGRATED_OVER_MOMENTUM_HPP*/
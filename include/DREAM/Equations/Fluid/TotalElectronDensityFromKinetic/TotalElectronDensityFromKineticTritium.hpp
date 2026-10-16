#ifndef _DREAM_EQUATIONS_TOTAL_ELECTRON_DENSITY_FROM_KINETIC_TRITIUM_HPP
#define _DREAM_EQUATIONS_TOTAL_ELECTRON_DENSITY_FROM_KINETIC_TRITIUM_HPP

#include "FVM/Equation/EquationTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/Equations/Kinetic/TritiumSource.hpp"

namespace DREAM {
    class TotalElectronDensityFromKineticTritium 
        : public FVM::EquationTerm {
    public:
        len_t id_nT;
        len_t indT;
        real_t pLower, pUpper, scaleFactor;
        real_t *weights = nullptr;
        bool hasBeenInitialized = false;
        
        
        TotalElectronDensityFromKineticTritium(FVM::Grid* g, real_t pLower, real_t pUpper, FVM::UnknownQuantityHandler *u, 
            IonHandler *ions, len_t iIon, real_t scaleFactor = 1.0);
        
        ~TotalElectronDensityFromKineticTritium();

        void AllocateWeights();
        void DeallocateWeights();
        void InitializeWeights();

        void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;

        bool SetJacobianBlock(const len_t , const len_t derivId, FVM::Matrix *jac, const real_t*);
        void SetMatrixElements(FVM::Matrix *mat, real_t* /*rhs*/) override;
        void SetVectorElements(real_t *vec, const real_t *x) override;
        void SetWeights();

        virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return GetNumberOfNonZerosPerRow(); }
        
    };
}

#endif/*_DREAM_EQUATIONS_TOTAL_ELECTRON_DENSITY_FROM_KINETIC_TRITIUM_HPP*/
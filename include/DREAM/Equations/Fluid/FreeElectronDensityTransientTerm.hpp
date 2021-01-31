#ifndef _DREAM_EQUATION_FLUID_FREE_ELECTRON_DENSITY_TRANSIENT_TERM_HPP
#define _DREAM_EQUATION_FLUID_FREE_ELECTRON_DENSITY_TRANSIENT_TERM_HPP

#include "DREAM/IonHandler.hpp"
#include "FVM/Equation/DiagonalLinearTerm.hpp"

namespace DREAM {
    class FreeElectronDensityTransientTerm : public FVM::DiagonalLinearTerm {
    private:
        IonHandler *ionHandler;
        const len_t id_ions;
        real_t scaleFactor;
        
        real_t dt;
        real_t *nions_prev;
    protected:
        virtual bool TermDependsOnUnknowns() override {return true;}
        virtual len_t GetNumberOfWeightsElements() override
            {return ionHandler->GetNzs()*grid->GetNCells();}        
        virtual void SetWeights() override;
    public:
        FreeElectronDensityTransientTerm(FVM::Grid*, IonHandler*, const len_t id_ions, real_t scaleFactor = 1.0);

        // This term sets 'nMultiples' diagonals in the matrix
        virtual len_t GetNumberOfNonZerosPerRow() const override { return this->ionHandler->GetNzs(); }
        
        virtual void Rebuild(const real_t t, const real_t dt, FVM::UnknownQuantityHandler*) override;        
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;

    };
}
#endif/*_DREAM_EQUATION_FLUID_FREE_ELECTRON_DENSITY_TRANSIENT_TERM_HPP*/

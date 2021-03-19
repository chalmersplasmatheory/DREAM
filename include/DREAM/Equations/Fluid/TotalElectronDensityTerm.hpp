#ifndef _DREAM_EQUATION_FLUID_TOTAL_ELECTRON_DENSITY_TERM_HPP
#define _DREAM_EQUATION_FLUID_TOTAL_ELECTRON_DENSITY_TERM_HPP

#include "DREAM/IonHandler.hpp"
#include "FVM/Equation/DiagonalLinearTerm.hpp"

namespace DREAM {
    class TotalElectronDensityTerm : public FVM::DiagonalLinearTerm {
    private:
        IonHandler *ionHandler;
        real_t scaleFactor;
    protected:

        virtual len_t GetNumberOfWeightsElements() override
            {return ionHandler->GetNzs()*grid->GetNCells();}
        
        virtual void SetWeights() override;
    public:
        TotalElectronDensityTerm(FVM::Grid*, IonHandler*, real_t scaleFactor = 1.0);

        virtual len_t GetNumberOfNonZerosPerRow() const override {
            // The matrix will consists of 'nMultiples' blocks with
            // one main diagonal each...
            return this->ionHandler->GetNzs();
        }
        
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;

    };
}
#endif/*_DREAM_EQUATION_FLUID_TOTAL_ELECTRON_DENSITY_TERM_HPP*/

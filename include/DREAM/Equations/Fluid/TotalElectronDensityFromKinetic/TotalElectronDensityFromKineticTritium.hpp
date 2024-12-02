#ifndef _DREAM_EQUATIONS_TOTAL_ELECTRON_DENSITY_FROM_KINETIC_TRITIUM_HPP
#define _DREAM_EQUATIONS_TOTAL_ELECTRON_DENSITY_FROM_KINETIC_TRITIUM_HPP

#include "FVM/Equation/DiagonalQuadraticTerm.hpp"

namespace DREAM {
    class TotalElectronDensityFromKineticTritium 
        : public FVM::DiagonalQuadraticTerm {
    public:
        real_t pLower, pUpper, scaleFactor;
        TotalElectronDensityFromKineticTritium(FVM::Grid* g, real_t pLower, real_t pUpper, FVM::UnknownQuantityHandler *u, real_t scaleFactor = 1.0) 
            : FVM::DiagonalQuadraticTerm(g,u->GetUnknownID(OptionConstants::UQTY_N_TOT),u), pLower(pLower), pUpper(pUpper), scaleFactor(scaleFactor) {}

        virtual void SetWeights() override {
            for(len_t i = 0; i<grid->GetNCells(); i++)
                weights[i] = scaleFactor * TritiumSource::EvaluateTotalTritiumNumber(pLower, pUpper);
        }
    };
}

#endif/*_DREAM_EQUATIONS_TOTAL_ELECTRON_DENSITY_FROM_KINETIC_TRITIUM_HPP*/
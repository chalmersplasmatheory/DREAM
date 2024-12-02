#ifndef _DREAM_EQUATIONS_TOTAL_ELECTRON_DENSITY_FROM_KINETIC_AVALANCHE_SOURCE_RP_HPP
#define _DREAM_EQUATIONS_TOTAL_ELECTRON_DENSITY_FROM_KINETIC_AVALANCHE_SOURCE_RP_HPP

#include "FVM/Equation/DiagonalQuadraticTerm.hpp"

/**
 * Implementation of an equation term which represents the total
 * number of electrons created by the kinetic Rosenbluth-Putvinski source
 */ 
namespace DREAM {
    class TotalElectronDensityFromKineticAvalancheRP : public FVM::DiagonalQuadraticTerm {
    public:
        real_t pLower, pUpper, scaleFactor;
        TotalElectronDensityFromKineticAvalancheRP(FVM::Grid* g, real_t pLower, real_t pUpper, FVM::UnknownQuantityHandler *u, real_t scaleFactor = 1.0) 
            : FVM::DiagonalQuadraticTerm(g,u->GetUnknownID(OptionConstants::UQTY_N_TOT),u), pLower(pLower), pUpper(pUpper), scaleFactor(scaleFactor) {}

        virtual void SetWeights() override {
            for(len_t i = 0; i<grid->GetNCells(); i++)
                weights[i] = scaleFactor * AvalancheSourceRP::EvaluateNormalizedTotalKnockOnNumber(pLower, pUpper);
        }
    };
} 

#endif/*_DREAM_EQUATIONS_TOTAL_ELECTRON_DENSITY_FROM_KINETIC_AVALANCHE_SOURCE_RP_HPP*/
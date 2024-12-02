#ifndef _DREAM_EQUATIONS_TOTAL_ELECTRON_DENSITY_FROM_KINETIC_COMPTON_HPP
#define _DREAM_EQUATIONS_TOTAL_ELECTRON_DENSITY_FROM_KINETIC_COMPTON_HPP

#include "FVM/Equation/DiagonalQuadraticTerm.hpp"
#include "FVM/Interpolator1D.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Equations/Kinetic/ComptonSource.hpp"

namespace DREAM {
    class TotalElectronDensityFromKineticCompton 
        : public FVM::DiagonalQuadraticTerm {
    private: 
        len_t limit;
        gsl_integration_workspace * wp;
        gsl_integration_workspace * wpOut;
    public:
        real_t pLower, pUpper;
        real_t integratedComptonSpectrum, C1, C2, C3;
        FVM::Interpolator1D *comptonPhotonFlux;
        real_t scaleFactor;
        real_t photonFlux;
        TotalElectronDensityFromKineticCompton(FVM::Grid* g, real_t pLower, real_t pUpper, FVM::UnknownQuantityHandler *u, FVM::Interpolator1D *comptonPhotonFlux, 
                real_t integratedComptonSpectrum, real_t C1, real_t C2, real_t C3, real_t scaleFactor = 1.0);
        
        virtual void Rebuild(const real_t t, const real_t, FVM::UnknownQuantityHandler*) override;

        virtual void SetWeights() override;
    };
}

#endif/*_DREAM_EQUATIONS_TOTAL_ELECTRON_DENSITY_FROM_KINETIC_COMPTON_HPP*/
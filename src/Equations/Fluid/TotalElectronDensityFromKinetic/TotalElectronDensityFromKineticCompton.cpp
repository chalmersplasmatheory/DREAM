#include "DREAM/Equations/Fluid/TotalElectronDensityFromKinetic/TotalElectronDensityFromKineticCompton.hpp"

/**
 * Implementation of an equation term which represents the total
 * number of electrons created by the kinetic Compton scattering 
 * source
 */ 
using namespace DREAM;

TotalElectronDensityFromKineticCompton::TotalElectronDensityFromKineticCompton(
    FVM::Grid* g, real_t pLower, real_t pUpper,
    FVM::Interpolator1D *comptonPhotonFlux, real_t integratedComptonSpectrum, 
    real_t C1, real_t C2, real_t C3, real_t scaleFactor
) : FVM::DiagonalLinearTerm(g), pLower(pLower), pUpper(pUpper), 
        integratedComptonSpectrum(integratedComptonSpectrum), C1(C1), C2(C2), C3(C3), comptonPhotonFlux(comptonPhotonFlux), scaleFactor(scaleFactor) 
{
    this->limit = 1000;
    this->wp = gsl_integration_workspace_alloc(limit);
    this->wpOut = gsl_integration_workspace_alloc(limit);
}

TotalElectronDensityFromKineticCompton::~TotalElectronDensityFromKineticCompton() {
    if (wp)
        gsl_integration_workspace_free(wp);
    if (wpOut) 
        gsl_integration_workspace_free(wpOut);
}

void TotalElectronDensityFromKineticCompton::Rebuild(const real_t t, const real_t, FVM::UnknownQuantityHandler*) {
    this->photonFlux = this->comptonPhotonFlux->Eval(t)[0];
    this->DiagonalLinearTerm::Rebuild(t,0,nullptr);
    this->SetWeights();
}
void TotalElectronDensityFromKineticCompton::SetWeights() {
    struct DREAM::ComptonSource::intparams params = {this->limit, this->wp, this->integratedComptonSpectrum, this->C1, this->C2, this->C3};
    struct DREAM::ComptonSource::intparams paramsOut = {this->limit, this->wpOut, this->integratedComptonSpectrum, this->C1, this->C2, this->C3};

    const real_t nCompton = ComptonSource::EvaluateTotalComptonNumber(pLower, &params, &paramsOut, pUpper);
    const real_t w = scaleFactor * this->photonFlux * nCompton;
    for(len_t i = 0; i<grid->GetNCells(); i++)
        weights[i] = w;
}

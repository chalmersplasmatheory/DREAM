/**
 * Implementation of equation term that sums ion densities over charge 
 * number to produce the free electron density at each radius.
 */

#include "DREAM/Equations/Fluid/TotalElectronDensityTerm.hpp"


using namespace DREAM;

/**
 * Constructor
 */
TotalElectronDensityTerm::TotalElectronDensityTerm(FVM::Grid *g, IonHandler *ih, real_t scaleFactor)
    : DiagonalLinearTerm(g), ionHandler(ih), scaleFactor(scaleFactor) {
    
    this->DiagonalTerm::SetName("TotalElectronDensityTerm");
}


/**
 * Set the linear operator matrix elements corresponding to this term.
 */
void TotalElectronDensityTerm::SetMatrixElements(FVM::Matrix *mat, real_t*) {
    len_t N = this->DiagonalTerm::grid->GetNCells();
    len_t nMultiples = ionHandler->GetNzs();

    for (len_t i = 0; i < N; i++)
        for(len_t n=0; n<nMultiples; n++)
            mat->SetElement(i, i+N*n, weights[N*n+i]);
}

/**
 * Set function vector for this term.
 */
void TotalElectronDensityTerm::SetVectorElements(real_t *vec, const real_t *x) {
    len_t N = this->DiagonalTerm::grid->GetNCells();
    len_t nMultiples = ionHandler->GetNzs();

    for (len_t i = 0; i < N; i++)
        for(len_t n=0; n<nMultiples; n++)
            vec[i] += weights[n*N+i] * x[n*N+i];
}

/**
 * Implementation of weights for this diagonal term
 */
void TotalElectronDensityTerm::SetWeights(){
    len_t N = grid->GetNCells();
    const len_t *Zs = ionHandler->GetZs();
    len_t NZ = ionHandler->GetNZ();
    for(len_t i=0; i<N; i++)
        for(len_t iz=0; iz<NZ; iz++)
            for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                len_t n = ionHandler->GetIndex(iz,Z0);
                weights[N*n+i] = scaleFactor*Zs[iz];
            }
}

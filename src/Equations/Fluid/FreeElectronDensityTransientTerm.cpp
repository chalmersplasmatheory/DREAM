/**
 * Implementation of equation term that sums ion densities over charge 
 * number to produce the change in the free electron density in the 
 * coming time step.
 */

#include "DREAM/Equations/Fluid/FreeElectronDensityTransientTerm.hpp"

using namespace DREAM;

/**
 * Constructor
 */
FreeElectronDensityTransientTerm::FreeElectronDensityTransientTerm(FVM::Grid *g, IonHandler *ih, const len_t id_ions, real_t scaleFactor)
    : DiagonalLinearTerm(g), ionHandler(ih), id_ions(id_ions), scaleFactor(scaleFactor) {
    
    this->DiagonalTerm::SetName("FreeElectronDensityTransientTerm");
}


/**
 * Rebuild the transient term.
 *
 * dt: Length of next time step to take.
 */
void FreeElectronDensityTransientTerm::Rebuild(const real_t t, const real_t dt, FVM::UnknownQuantityHandler *uqty) {
    this->dt = dt;
    this->nions_prev = uqty->GetUnknownDataPrevious(this->id_ions);
    this->DiagonalLinearTerm::Rebuild(t,dt,uqty);
}


/**
 * Set the linear operator matrix elements corresponding to this term.
 */
void FreeElectronDensityTransientTerm::SetMatrixElements(FVM::Matrix *mat, real_t *rhs) {
	// Transient term disabled? (mainly applicable for
	// equation system initializer, which seeks solutions for
	// which dX/dt -> 0)
	if (this->dt == 0)
		return;

    len_t N = this->DiagonalTerm::grid->GetNCells();
    len_t nMultiples = ionHandler->GetNzs();
    for (len_t i = 0; i < N; i++)
        for(len_t n=0; n<nMultiples; n++)
            mat->SetElement(i,n*N+i, weights[n*N+i]);
    if (rhs != nullptr)
        for (len_t i = 0; i < N; i++)
            for(len_t n=0; n<nMultiples; n++)
                rhs[i] -= weights[n*N+i]*nions_prev[n*N+i];
}


/**
 * Set function vector for this term.
 */
void FreeElectronDensityTransientTerm::SetVectorElements(real_t *vec, const real_t *nions) {
	// Transient term disabled? (mainly applicable for
	// equation system initializer, which seeks solutions for
	// which dX/dt -> 0)
	if (this->dt == 0)
		return;

    len_t N = this->DiagonalTerm::grid->GetNCells();
    len_t nMultiples = ionHandler->GetNzs();
    for (len_t i = 0; i < N; i++)
        for(len_t n=0; n<nMultiples; n++)
            vec[i] += weights[n*N+i] * (nions[n*N+i]-nions_prev[n*N+i]);
}


/**
 * Implementation of weights for this diagonal term
 */
void FreeElectronDensityTransientTerm::SetWeights(){
	// Transient term disabled? (mainly applicable for
	// equation system initializer, which seeks solutions for
	// which dX/dt -> 0)
	if (this->dt == 0)
		return;

    len_t N = grid->GetNCells();
    for(len_t i=0; i<N; i++)
        for(len_t iz=0; iz<ionHandler->GetNZ(); iz++)
            for(len_t Z0=0; Z0<=ionHandler->GetZ(iz); Z0++){
                len_t n = ionHandler->GetIndex(iz,Z0);
                weights[N*n+i] = scaleFactor*Z0 / dt;
            }
}

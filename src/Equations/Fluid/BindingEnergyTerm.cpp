#include "DREAM/Equations/Fluid/BindingEnergyTerm.hpp"
#include "DREAM/NIST.hpp"
#include "DREAM/Constants.hpp"

using namespace DREAM;

/**
 * Constructor.
 */
BindingEnergyTerm::BindingEnergyTerm(FVM::Grid* g, IonHandler *ionHandler, NIST *nist)
    : FVM::DiagonalLinearTerm(g), ionHandler(ionHandler), nist(nist) {
    
    SetName("BindingEnergyTerm");
}

/** 
 * Sets weights to the binding potential, which is _minus_
 * the binding energies since they are defined as positive.
 */
void BindingEnergyTerm::SetWeights(){
    len_t N = grid->GetNCells();
    len_t nZ = ionHandler->GetNZ();
    const len_t *Zs = ionHandler->GetZs();

    for(len_t iz = 0; iz<nZ; iz++){
        for(len_t Z0 = 0; Z0<=Zs[iz]; Z0++){
            len_t n  = ionHandler->GetIndex(iz,Z0);
            real_t w = nist->GetBindingEnergy(Zs[iz],Z0) * Constants::ec;

            for(len_t i = 0; i<N; i++)
                weights[n*N+i] = -w;
        }
    }
}


/**
 * Set matrix for all nnz ion species
 */
void BindingEnergyTerm::SetMatrixElements(FVM::Matrix* mat, real_t*){  
    len_t N   = grid->GetNCells();
    len_t nnz = ionHandler->GetNzs();

    for(len_t n=0; n<nnz;n++)
        for(len_t i=0; i<N; i++)
            mat->SetElement(i,n*N+i,weights[n*N+i]);
}

/**
 * Set vector for all nnz ion species
 */
void BindingEnergyTerm::SetVectorElements(real_t* vec, const real_t* x){ 
    len_t N   = grid->GetNCells();
    len_t nnz = ionHandler->GetNzs();

    for(len_t n=0; n < nnz; n++)
        for (len_t i = 0; i < N; i++)
            vec[i] += weights[n*N+i] * x[n*N+i];
}


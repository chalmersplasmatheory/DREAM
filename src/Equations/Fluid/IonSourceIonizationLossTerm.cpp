#include "DREAM/Equations/Fluid/IonSourceIonizationLossTerm.hpp"
#include "DREAM/NIST.hpp"
#include "DREAM/Constants.hpp"

using namespace DREAM;

/**
 * Constructor.
 */
IonSourceIonizationLossTerm::IonSourceIonizationLossTerm(FVM::Grid* g, IonHandler *ionHandler, NIST *nist, EquationTerm *ionSource)
    : FVM::EquationTerm(g), ionHandler(ionHandler), nist(nist), ionSource(ionSource) {

    len_t nZ = ionHandler->GetNZ();
    const len_t *Zs = ionHandler->GetZs();

    for(len_t iz = 0; iz<nZ; iz++){
        for(len_t Z0 = 0; Z0<=Zs[iz]; Z0++){
            real_t ETotIoniz[ionHandler->GetIndex(iz,Z0)] = (nist->GetBindingEnergy(0,Z0)-nist->GetBindingEnergy(Zs[iz],Z0)) * Constants::ec;
        }
    }
}

void IonSourceIonizationLossTerm::SetJacobianBlock(FVM::Matrix *jac, real_t *x){
    ionSource->SetJacobianBlock(jac,x);
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


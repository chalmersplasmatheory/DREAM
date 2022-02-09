#include "DREAM/Constants.hpp"
#include "DREAM/Equations/Fluid/IonNeutralAdvectionDiffusionTerm.hpp"
#include "DREAM/IonHandler.hpp"
//#include "DREAM/NotImplementedException.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM;

template<class T>
IonNeutralAdvectionDiffusionTerm<T>::IonNeutralAdvectionDiffusionTerm(FVM::Grid *grid, IonHandler *ihdl, const len_t iIon, bool allocCoefficients
	) : IonEquationTerm<T>(grid, ihdl, iIon, allocCoefficients){}

template<class T>
IonNeutralAdvectionDiffusionTerm<T>::~IonNeutralAdvectionDiffusionTerm(){}

template<class T>
bool IonNeutralAdvectionDiffusionTerm<T>::SetCSJacobianBlock(
    const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t *f,
    const len_t /*iIon*/, const len_t Z0, const len_t rOffset
){
    bool contributes = false;
    
    if(!this->HasJacobianContribution(derivId) || Z0>0)
        return contributes;
        
    len_t rowOffset0 = jac->GetRowOffset();
    len_t colOffset0 = jac->GetColOffset();
    jac->SetOffset(rowOffset0+rOffset,colOffset0);
    
    contributes = this->T::SetJacobianBlock(uqtyId, derivId, jac, f);
    
    jac->SetOffset(rowOffset0,colOffset0);

    return contributes;
}


/**
 * Sets the matrix elements of this equation term
 */
template<class T>
void IonNeutralAdvectionDiffusionTerm<T>::SetCSMatrixElements(
    FVM::Matrix *mat, real_t *rhs, const len_t /*iIon*/, const len_t Z0, const len_t rOffset
) {
    if(Z0>0)
        return;
    
    len_t rowOffset0 = mat->GetRowOffset();
    len_t colOffset0 = mat->GetColOffset();
    mat->SetOffset(rowOffset0+rOffset,colOffset0);
    this->T::SetMatrixElements(mat,rhs);
    mat->SetOffset(rowOffset0,colOffset0);
}


/**
 * Sets vector elements for this ion and charge state 
 */
template<class T>
void IonNeutralAdvectionDiffusionTerm<T>::SetCSVectorElements(
    real_t *vec, const real_t *f, const len_t /*iIon*/, const len_t Z0, const len_t rOffset
) {
    if(Z0>0)
        return;
    this->T::SetVectorElements(vec+rOffset, f+rOffset);
}

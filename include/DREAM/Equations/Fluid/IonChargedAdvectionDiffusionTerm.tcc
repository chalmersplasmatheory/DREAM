#include "DREAM/Constants.hpp"
#include "DREAM/Equations/Fluid/IonChargedAdvectionDiffusionTerm.hpp"
#include "DREAM/IonHandler.hpp"
//#include "DREAM/NotImplementedException.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM;

template<class T>
IonChargedAdvectionDiffusionTerm<T>::IonChargedAdvectionDiffusionTerm(FVM::Grid *grid, IonHandler *ihdl, const len_t iIon, bool allocCoefficients
	) : IonEquationTerm<T>(grid, ihdl, iIon, allocCoefficients)
{
	Allocate();
}

template<class T>
IonChargedAdvectionDiffusionTerm<T>::~IonChargedAdvectionDiffusionTerm(){
	Deallocate();
}

template<class T>
void IonChargedAdvectionDiffusionTerm<T>::Allocate(){
    Deallocate();

    const len_t nr = this->grid->GetNr();

	// Allocate memory for one diffusion coefficient of every charged state (excluding neutral) 
	// and every radial grid point
    this->CoeffsAllCS = new real_t*[this->Zion];
    for(len_t Z0=1; Z0<=this->Zion; Z0++)
        this->CoeffsAllCS[Z0-1] = new real_t[nr+1];
}

template<class T>
void IonChargedAdvectionDiffusionTerm<T>::Deallocate() {
    if(CoeffsAllCS == nullptr)
        return;
    for(len_t Z0=1; Z0<=this->Zion; Z0++)
        delete [] this->CoeffsAllCS[Z0-1];
    delete [] this->CoeffsAllCS;
}

template<class T>
void IonChargedAdvectionDiffusionTerm<T>::Rebuild(const real_t t, const real_t, FVM::UnknownQuantityHandler*){
	SetCoeffsAllCS(t);
	SetDiffCoeffsAllCS(t);
}

template<class T>
bool IonChargedAdvectionDiffusionTerm<T>::SetCSJacobianBlock(
    const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t *f,
    const len_t iIon, const len_t Z0, const len_t rOffset
){
    bool contributes = false;
    if(Z0<1)
        return contributes;
        
    if(uqtyId == derivId){
        contributes = true;
        this->SetCSMatrixElements(jac, nullptr, iIon, Z0, rOffset);
    }
    
    if(!this->HasJacobianContribution(derivId))
        return contributes;
        
    len_t rowOffset0 = jac->GetRowOffset();
    len_t colOffset0 = jac->GetColOffset();
    jac->SetOffset(rowOffset0+rOffset,colOffset0);
    
    this->Z0ForPartials = Z0;
    SetCoeffs(Z0);
    contributes = this->T::SetJacobianBlock(uqtyId, derivId, jac, f);
    
    jac->SetOffset(rowOffset0,colOffset0);

    return contributes;
}


/**
 * Sets the matrix elements of this equation term
 */
template<class T>
void IonChargedAdvectionDiffusionTerm<T>::SetCSMatrixElements(
    FVM::Matrix *mat, real_t *rhs, const len_t /*iIon*/, const len_t Z0, const len_t rOffset
) {
	if(Z0<1)
		return;
		
    SetCoeffs(Z0); 
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
void IonChargedAdvectionDiffusionTerm<T>::SetCSVectorElements(
    real_t *vec, const real_t *f, const len_t /*iIon*/, const len_t Z0, const len_t rOffset
) {
	if(Z0<1)
		return;
		
    SetCoeffs(Z0); 
    this->T::SetVectorElements(vec+rOffset, f+rOffset);
    // Indexing was originally copied from IonKineticIonizationTerm, 
    // but there the last rOffset was missing. Should it be like that?
}

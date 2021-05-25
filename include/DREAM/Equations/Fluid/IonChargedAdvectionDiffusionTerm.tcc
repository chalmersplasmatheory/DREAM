#include "DREAM/Constants.hpp"
#include "DREAM/Equations/Fluid/IonChargedAdvectionDiffusionTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/NotImplementedException.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM;

template<class T>
IonChargedAdvectionDiffusionTerm::IonChargedAdvectionDiffusionTerm(FVM::Grid *g, IonHandler *ihdl, const len_t iIon
	) : IonEquationTerm<T>(g, ihdl, iIon), T(g, true)
{
	Allocate();
}

template<class T>
IonChargedAdvectionDiffusionTerm::~IonChargedAdvectionDiffusionTerm(){
	Deallocate();
}

template<class T>
void IonChargedAdvectionDiffusionTerm::Allocate(){
    Deallocate();

    const len_t nr = this->g->GetNr();

	// Allocate memory for one diffusion coefficient of every charged state (excluding neutral) 
	// and every radial grid point
    this->CoeffsAllCS = new real_t*[Zion];
    for(len_t Z0=1; Z0<=Zion; Z0++)
        this->CoeffsAllCS[Z0-1] = new real_t[nr+1];
}

template<class T>
void IonChargedAdvectionDiffusionTerm::Deallocate() {
    if(CoeffsAllCS == nullptr)
        return;
    for(len_t Z0=1; Z0<=Zion; Z0++)
        delete [] this->CoeffsAllCS[Z0-1];
    delete [] this->CoeffsAllCS;
}

template<class T>
void IonChargedAdvectionDiffusionTerm::Rebuild(const real_t t, const real_t, FVM::UnknownQuantityHandler*){
	SetCoeffsAllCS(t);
	SetDiffCoeffsAllCS(t);
}

template<class T>
bool IonChargedAdvectionDiffusionTerm::SetCSJacobianBlock(
    const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t *f,
    const len_t iIon, const len_t Z0, const len_t rOffset
){
    bool contributes = false;
    
    if(!HasJacobianContribution(derivId) || Z0<1)
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
void IonChargedAdvectionDiffusionTerm::SetCSMatrixElements(
    FVM::Matrix *mat, real_t *rhs, const len_t /*iIon*/, const len_t Z0, const len_t rOffset
) {
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
void IonChargedAdvectionDiffusionTerm::SetCSVectorElements(
    real_t *vec, const real_t *f, const len_t /*iIon*/, const len_t Z0, const len_t rOffset
) {
    SetCoeffs(Z0); 
    this->T::SetVectorElements(vec+rOffset, f);
}

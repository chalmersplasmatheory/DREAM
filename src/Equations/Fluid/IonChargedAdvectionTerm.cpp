#include "DREAM/Equations/Fluid/IonChargedAdvectionTerm.hpp"
#include "FVM/Grid/Grid.hpp"

#include <iostream>

using namespace DREAM;

IonChargedAdvectionTerm::IonChargedAdvectionTerm(FVM::Grid *grid, IonHandler *ions, const len_t iIon, 
    bool allocCoefficients, FVM::AdvectionInterpolationCoefficient::adv_interpolation intp,
    OptionConstants::adv_jacobian_mode jac_mode, len_t id, real_t damping_factor) : IonChargedAdvectionDiffusionTerm<FVM::AdvectionTerm>(grid, ions, iIon, allocCoefficients){
	    AllocateInterpolationCoefficients();
	    SetAdvectionInterpolationMethod(intp, jac_mode, id, damping_factor);
	}
	
IonChargedAdvectionTerm::~IonChargedAdvectionTerm(){}

void IonChargedAdvectionTerm::AllocateInterpolationCoefficients(){
    if (!this->interpolationCoefficientsShared)
        DeallocateInterpolationCoefficients();
        
    const len_t nr = this->grid->GetNr();

    len_t idx = this->ions->GetIndex(iIon, 0);
    for(len_t Z0 = 1; Z0<=Zion; Z0++,idx++)
        this->deltars.push_back(new FVM::AdvectionInterpolationCoefficient(grid, FVM::FLUXGRIDTYPE_RADIAL, FVM::AdvectionInterpolationCoefficient::AD_BC_MIRRORED, FVM::AdvectionInterpolationCoefficient::AD_BC_DIRICHLET, idx*nr));
    
    this->interpolationCoefficientsShared = false;
}

void IonChargedAdvectionTerm::DeallocateInterpolationCoefficients(){
    for (auto it = deltars.begin(); it != deltars.end(); it++)
        if((*it)!=nullptr){
            delete (*it);
            (*it) = nullptr;
        }
}

void IonChargedAdvectionTerm::RebuildInterpolationCoefficients(FVM::UnknownQuantityHandler* uqty, real_t**drr, real_t** /*d11*/, real_t** /*d22*/){
    len_t Z0 = 1;
    for (auto it = deltars.begin(); it != deltars.end(); it++, Z0++){
        SetCoeffs(Z0);
        (*it)->SetCoefficient(this->fr, drr, uqty,  dampingWithIteration*fluxLimiterDampingFactor);
    }
}

void IonChargedAdvectionTerm::SetAdvectionInterpolationMethod(
                FVM::AdvectionInterpolationCoefficient::adv_interpolation intp,
                OptionConstants::adv_jacobian_mode jac_mode, 
                len_t id, real_t damping_factor
            ){
    this->fluxLimiterDampingFactor = damping_factor;
    for (auto it = deltars.begin(); it != deltars.end(); it++){
        (*it)->SetInterpolationMethod(intp); 
        (*it)->SetUnknownId(id);
        (*it)->SetJacobianMode(jac_mode);
    }
}            
            
void IonChargedAdvectionTerm::Rebuild(const real_t t, const real_t dt , FVM::UnknownQuantityHandler* uqty){
    this->IonChargedAdvectionDiffusionTerm<FVM::AdvectionTerm>::Rebuild(t, dt, uqty);
    
    this->FVM::AdvectionTerm::RebuildFluxLimiterDamping(t, dt);
    this->RebuildInterpolationCoefficients(uqty, nullptr, nullptr, nullptr);
}

bool IonChargedAdvectionTerm::SetCSJacobianBlock(
		    const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t *nions,
		    const len_t iIon, const len_t Z0, const len_t rOffset
		){
		if(Z0<1)
		    return false;
		    
		deltar = deltars[Z0-1];
		return this->IonChargedAdvectionDiffusionTerm<FVM::AdvectionTerm>::SetCSJacobianBlock(uqtyId, derivId, jac, nions, iIon, Z0, rOffset);
}

void IonChargedAdvectionTerm::SetCSMatrixElements(
		    FVM::Matrix *mat, real_t *rhs, const len_t iIon, const len_t Z0, const len_t rOffset
		){
		if(Z0<1)
		    return;
		    
		deltar = deltars[Z0-1];
		this->IonChargedAdvectionDiffusionTerm<FVM::AdvectionTerm>::SetCSMatrixElements(mat, rhs, iIon, Z0, rOffset);
}

void IonChargedAdvectionTerm::SetCSVectorElements(
		    real_t *vec, const real_t *nions, const len_t iIon, const len_t Z0, const len_t rOffset
		){
		if(Z0<1)
		    return;
		    
		deltar = deltars[Z0-1];
		this->IonChargedAdvectionDiffusionTerm<FVM::AdvectionTerm>::SetCSVectorElements(vec, nions, iIon, Z0, rOffset);
}
		
		
		
		

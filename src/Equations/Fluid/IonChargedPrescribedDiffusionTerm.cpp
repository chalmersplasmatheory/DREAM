#include "DREAM/Constants.hpp"
#include "DREAM/Equations/Fluid/IonChargedAdvectionDiffusionTerm.hpp"
#include "DREAM/Equations/Fluid/IonChargedPrescribedDiffusionTerm.hpp"
#include "DREAM/Equations/Fluid/IonEquationTerm.hpp"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "DREAM/IonHandler.hpp"
//#include "DREAM/NotImplementedException.hpp"
#include "FVM/Grid/Grid.hpp"

#include <iostream>


using namespace DREAM;

IonChargedPrescribedDiffusionTerm::IonChargedPrescribedDiffusionTerm(FVM::Grid *grid, IonHandler *ions,
	const len_t iIon, bool allocCoefficients, len_t offset, MultiInterpolator1D* DrrPrescribed
	) : IonChargedAdvectionDiffusionTerm<FVM::DiffusionTerm>(grid, 
	ions, iIon, allocCoefficients), offset(offset), DrrPrescribed(DrrPrescribed) {
	
    SetName("IonChargedPrescribedDiffusionTerm");

}

IonChargedPrescribedDiffusionTerm::~IonChargedPrescribedDiffusionTerm(){}


void IonChargedPrescribedDiffusionTerm::SetCoeffsAllCS(const real_t t){
    for(len_t Z0=1; Z0<=this->Zion; Z0++){
        const real_t *DrrZ0 = DrrPrescribed->Eval(offset+Z0-1,t); 
        for(len_t ir=0; ir<nr+1; ir++)
            CoeffsAllCS[Z0-1][ir] = DrrZ0[ir];
    }
}

void IonChargedPrescribedDiffusionTerm::SetDiffCoeffsAllCS(const real_t){}

void IonChargedPrescribedDiffusionTerm::SetCoeffs(const len_t Z0){
	if(Z0<1)
		return;
		
	const len_t nr = this->grid->GetNr();
	for(len_t ir=0; ir<nr+1; ir++){
	    // As the IonTransietTerm is defined with a different sign 
	    // than ordinary transient terms for some reason, 
	    // we have to include a minus sign here
		Drr(ir,0,0)=-CoeffsAllCS[Z0-1][ir];
	}
}


void IonChargedPrescribedDiffusionTerm::SetPartialDiffusionTerm(len_t, len_t){}





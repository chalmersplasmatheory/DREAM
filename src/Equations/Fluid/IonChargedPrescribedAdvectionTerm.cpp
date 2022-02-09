#include "DREAM/Constants.hpp"
#include "DREAM/Equations/Fluid/IonChargedAdvectionTerm.hpp"
#include "DREAM/Equations/Fluid/IonChargedAdvectionDiffusionTerm.hpp"
#include "DREAM/Equations/Fluid/IonChargedPrescribedAdvectionTerm.hpp"
#include "DREAM/Equations/Fluid/IonEquationTerm.hpp"
#include "FVM/Equation/AdvectionTerm.hpp"
#include "DREAM/IonHandler.hpp"
//#include "DREAM/NotImplementedException.hpp"
#include "FVM/Grid/Grid.hpp"

#include <iostream>


using namespace DREAM;

IonChargedPrescribedAdvectionTerm::IonChargedPrescribedAdvectionTerm(FVM::Grid *grid, IonHandler *ions, const len_t iIon,
    bool allocCoefficients, FVM::AdvectionInterpolationCoefficient::adv_interpolation intp,
    OptionConstants::adv_jacobian_mode jac_mode, len_t id, real_t damping_factor, len_t offset, MultiInterpolator1D* FrPrescribed
	) : IonChargedAdvectionTerm(grid, 
	ions, iIon, allocCoefficients, intp, jac_mode, id, damping_factor), offset(offset), FrPrescribed(FrPrescribed) {
	
    SetName("IonChargedPrescribedAdvectionTerm");

}

IonChargedPrescribedAdvectionTerm::~IonChargedPrescribedAdvectionTerm(){}


void IonChargedPrescribedAdvectionTerm::SetCoeffsAllCS(const real_t t){
    for(len_t Z0=1; Z0<=this->Zion; Z0++){
        const real_t *FrZ0 = FrPrescribed->Eval(offset+Z0-1,t); 
        for(len_t ir=0; ir<nr+1; ir++)
            CoeffsAllCS[Z0-1][ir] = FrZ0[ir];
    }
}

void IonChargedPrescribedAdvectionTerm::SetDiffCoeffsAllCS(const real_t){}

void IonChargedPrescribedAdvectionTerm::SetCoeffs(const len_t Z0){
	if(Z0<1)
		return;
		
	const len_t nr = this->grid->GetNr();
	for(len_t ir=0; ir<nr+1; ir++){
	    // As the IonTransietTerm is defined with a different sign 
	    // than ordinary transient terms for some reason, 
	    // we have to include a minus sign here
		Fr(ir,0,0)=-CoeffsAllCS[Z0-1][ir];
	}
}


void IonChargedPrescribedAdvectionTerm::SetPartialAdvectionTerm(len_t, len_t){}





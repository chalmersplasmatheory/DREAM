#include "DREAM/Constants.hpp"
#include "DREAM/Equations/Fluid/IonNeutralAdvectionDiffusionTerm.hpp"
#include "DREAM/Equations/Fluid/IonNeutralPrescribedDiffusionTerm.hpp"
#include "DREAM/Equations/Fluid/IonEquationTerm.hpp"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "DREAM/IonHandler.hpp"
//#include "DREAM/NotImplementedException.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM;

IonNeutralPrescribedDiffusionTerm::IonNeutralPrescribedDiffusionTerm(FVM::Grid *grid, IonHandler *ions, bool allocCoefficients,
	const len_t iIon, len_t offset, MultiInterpolator1D* DrrPrescribed
	) : IonNeutralAdvectionDiffusionTerm<FVM::DiffusionTerm>(grid, 
	ions, iIon, allocCoefficients), offset(offset), DrrPrescribed(DrrPrescribed) {
	
    SetName("IonNeutralPrescribedDiffusionTerm");

}

IonNeutralPrescribedDiffusionTerm::~IonNeutralPrescribedDiffusionTerm(){}


void IonNeutralPrescribedDiffusionTerm::Rebuild(const real_t t, const real_t, FVM::UnknownQuantityHandler*){
    const real_t *DrrEval = DrrPrescribed->Eval(offset,t); 
    for(len_t ir=0; ir<nr+1; ir++)
	    // As the IonTransietTerm is defined with a different sign 
	    // than ordinary transient terms for some reason, 
	    // we have to include a minus sign here
        Drr(ir,0,0) = -DrrEval[ir];

}

void IonNeutralPrescribedDiffusionTerm::SetPartialDiffusionTerm(len_t, len_t){}





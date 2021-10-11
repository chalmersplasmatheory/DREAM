#include "DREAM/Constants.hpp"
#include "DREAM/Equations/Fluid/IonNeutralAdvectionTerm.hpp"
#include "DREAM/Equations/Fluid/IonNeutralPrescribedAdvectionTerm.hpp"
#include "DREAM/Equations/Fluid/IonEquationTerm.hpp"
#include "FVM/Equation/AdvectionTerm.hpp"
#include "DREAM/IonHandler.hpp"
//#include "DREAM/NotImplementedException.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM;

IonNeutralPrescribedAdvectionTerm::IonNeutralPrescribedAdvectionTerm(FVM::Grid *grid, IonHandler *ions, const len_t iIon, 
    bool allocCoefficients, FVM::AdvectionInterpolationCoefficient::adv_interpolation intp,
    OptionConstants::adv_jacobian_mode jac_mode, len_t id, real_t damping_factor, len_t offset, MultiInterpolator1D* FrPrescribed
	) : IonNeutralAdvectionTerm(grid, 
	ions, iIon, allocCoefficients, intp, jac_mode, id, damping_factor), offset(offset), FrPrescribed(FrPrescribed) {
	
    SetName("IonNeutralPrescribedAdvectionTerm");

}

IonNeutralPrescribedAdvectionTerm::~IonNeutralPrescribedAdvectionTerm(){}


void IonNeutralPrescribedAdvectionTerm::SetCoeffs(const real_t t){
    const real_t *FrEval = FrPrescribed->Eval(offset,t); 
    for(len_t ir=0; ir<nr+1; ir++)
	    // As the IonTransietTerm is defined with a different sign 
	    // than ordinary transient terms for some reason, 
	    // we have to include a minus sign here
        Fr(ir,0,0) = -FrEval[ir];

}

void IonNeutralPrescribedAdvectionTerm::SetPartialAdvectionTerm(len_t, len_t){}





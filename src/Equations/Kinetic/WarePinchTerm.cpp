/**
 * This equation term implements the radial advection coefficient that 
results from the ware pinch effect:
 *
 *   Ar = <E*B>/(psi' * G *<1/R²>
 *
 */


#include "DREAM/Constants.hpp"
#include "DREAM/Equations/Kinetic/WarePinchTerm.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/UnknownQuantityHandler.hpp"


using namespace DREAM;


/**
 * Constructor.
 *
 * grid: Computational grid on which the equation term is defined.
 *   .
 */
WarePinchTerm::WarePinchTerm(
    FVM::Grid *grid, enum OptionConstants::momentumgrid_type mgtype, FVM::UnknownQuantityHandler * UQH
) : AdvectionTerm(grid), mgtype(mgtype)
 {

    SetName("WarePinchTerm");
    this->id_E_field = UQH->GetUnknownID(OptionConstants::UQTY_E_FIELD);
    AddUnknownForJacobian(UQH, id_E_field);
    this->dFrMatrix = new real_t[grid->GetNr() +1];
}

/**
 * Destructor.
 */
WarePinchTerm::~WarePinchTerm() {
    delete [] this->dFrMatrix;
}


/**
 * Rebuild this equation term.
 */
void WarePinchTerm::Rebuild(
    const real_t t, const real_t, FVM::UnknownQuantityHandler * UQH
) {

    // XXX Here we assume that all momentum grids are the same
    // at all radii...
    const len_t nr = this->grid->GetNr();
    auto mg = this->grid->GetMomentumGrid(0);

    real_t R0  = this->grid->GetRadialGrid()->GetR0();
    const len_t
        np1 = mg->GetNp1(),
        np2 = mg->GetNp2();
    const real_t *E_field = UQH->GetUnknownData(id_E_field); // <E*B>/<B^2>^0.5


    for (len_t ir = 0; ir < nr+1; ir++) {
        const real_t G = this->grid->GetRadialGrid()->GetBTorG_f(ir); //toroidal field function (I) G/R0
        const real_t ONE_OVER_R2 = this->grid->GetRadialGrid()->GetFSA_1OverR2_f(ir); // <1/R^2>
	    const real_t PSI = this->grid->GetRadialGrid()->GetPsiPrimeRef_f(ir); //d(psi)/dr
	    const real_t B2 = this->grid->GetRadialGrid()->GetFSA_B2_f(ir); // <B^2>  
	    real_t E_f = 0; // We need to evalutate E_field at the boundary instead of the centre
	    if (ir > 0 and ir < nr){
	        E_f = 0.5*(E_field[ir-1] + E_field[ir]);
	    }
	    else if(ir == 0){
	        E_f = E_field[0];
	    }
	    else if(ir == nr){
	        E_f = E_field[nr-1];
	    }    
	    for (len_t j = 0; j < np2; j++)  {              
            if (this->grid->IsTrapped_fr(ir, 0, j)){
                for (len_t i = 0; i < np1; i++) {
                    // Set Advection coefficient...
                    
                    real_t A =-std::sqrt(B2) * E_f/ (R0*G*PSI*ONE_OVER_R2);
                    Fr(ir, i, j) += A;  
                    
                }
            }   	
	    }
	    real_t dA =-std::sqrt(B2)/ (R0*G*PSI*ONE_OVER_R2);
        dFrMatrix[ir] = dA;
	}
}
// Set jacobian of the advection coefficients for this advection term
void WarePinchTerm::SetPartialAdvectionTerm(len_t /*derivId*/, len_t /*nMultiples*/){
    ResetDifferentiationCoefficients();
    const len_t nr = this->grid->GetNr();
    auto mg = this->grid->GetMomentumGrid(0);
    const len_t
        np1 = mg->GetNp1(),
        np2 = mg->GetNp2();
    for (len_t ir = 0; ir < nr+1; ir++) {
        for (len_t j = 0; j < np2; j++)  {              
            if (this->grid->IsTrapped_fr(ir, 0, j)){
                for (len_t i = 0; i < np1; i++) {
                dFr(ir, i, j, 0) = dFrMatrix[ir];
                }
            }
        }   
    }       
} 

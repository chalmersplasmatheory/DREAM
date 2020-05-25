/**
 * Implementation of the p*nu_s friction term in the kinetic equation.
 */

#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Equations/Kinetic/SlowingDownTerm.hpp"


using namespace DREAM;

/**
 * Constructor.
 */
SlowingDownTerm::SlowingDownTerm(
    FVM::Grid *g, CollisionQuantityHandler *cqh,
    enum OptionConstants::momentumgrid_type mgtype
) : FVM::AdvectionTerm(g) {

    this->gridtype = mgtype;
    this->nuS = cqh->GetNuS();
}

/**
 * Build the coefficients of this advection term.
 */
void SlowingDownTerm::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler *){
    const len_t nr = this->grid->GetNr();
 
    real_t *const* nu_s_f1 = nuS->GetValue_f1();
    real_t *const* nu_s_f2 = nuS->GetValue_f2();
  
    bool gridtypePXI, gridtypePPARPPERP;


    for (len_t ir = 0; ir < nr; ir++) {
        auto *mg = this->grid->GetMomentumGrid(ir);
        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();

        gridtypePXI        = (gridtype == OptionConstants::MOMENTUMGRID_TYPE_PXI);
        gridtypePPARPPERP  = (gridtype == OptionConstants::MOMENTUMGRID_TYPE_PPARPPERP);
        
        if (gridtypePXI || gridtypePPARPPERP) {
            for (len_t j = 0; j < np2; j++) {
                for (len_t i = 0; i < np1+1; i++) {
                    F1(ir, i, j) -= mg->GetP1_f(i) * nu_s_f1[ir][j*(np1+1)+i];
                }
            }
        }

        if (gridtypePPARPPERP) {
            for (len_t j = 0; j < np2+1; j++) {
                for (len_t i = 0; i < np1; i++) {
                    F2(ir, i, j) -= mg->GetP2_f(j) * nu_s_f2[ir][j*np1+i];
                }
            }
        }
    
    }
}



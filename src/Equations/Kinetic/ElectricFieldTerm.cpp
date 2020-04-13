/**
 * Implementation of the electric field advection term in the kinetic equation.
 */

#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include "DREAM/Equations/Kinetic/ElectricFieldTerm.hpp"
#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"


using namespace DREAM;

/**
 * Constructor.
 */
ElectricFieldTerm::ElectricFieldTerm(FVM::Grid *g, FVM::UnknownQuantityHandler *unknowns, enum OptionConstants::momentumgrid_type mgtype)
    : FVM::AdvectionTerm(g) {
        this->gridtype  = mgtype;
        this->grid      = g;
        this->id_Eterm  = unknowns->GetUnknownID( OptionConstants::UQTY_E_FIELD ); // E term should be <E*B>/sqrt(<B^2>)
}


/**
 * Build the coefficients of this advection (or diffusion) term.
 */
void ElectricFieldTerm::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler *x){
    const len_t nr = this->grid->GetNr();
    
    bool gridtypePXI, gridtypePPARPPERP;
    real_t xi0_f;
    real_t E_xi_bounceAvg_f1, E_xi_bounceAvg_f2;
    real_t sqrtB2OverB;
    real_t *E_term = x->GetUnknownData(id_Eterm);
    const real_t *xiAvgTerm_f1, *xiAvgTerm_f2;

    for (len_t ir = 0; ir < nr; ir++) {
        auto *mg = this->grid->GetMomentumGrid(ir);
        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();
        gridtypePXI         = (gridtype == OptionConstants::MOMENTUMGRID_TYPE_PXI);
        gridtypePPARPPERP   = (gridtype == OptionConstants::MOMENTUMGRID_TYPE_PPARPPERP);
        sqrtB2OverB =  sqrt(this->grid->GetRadialGrid()->GetFSA_B2(ir)) / this->grid->GetRadialGrid()->GetFSA_B(ir);
        xiAvgTerm_f1 = this->grid->GetRadialGrid()->GetBA_xi_f1(ir) ;
           ;
         
        //this->grid->GetRadialGrid()->GetBA_BOverBOverXi_f1(ir); 

        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1+1; i++) {
                E_xi_bounceAvg_f1 = Constants::ec * E_term[ir] * xiAvgTerm_f1[j*(np1+1)+i] * sqrtB2OverB;
                
                if (gridtypePXI) {
                    F1(ir, i, j)  += E_xi_bounceAvg_f1;
                } else if (gridtypePPARPPERP) {
                    xi0_f = mg->GetXi0_f1(i,j);
                    F1(ir, i, j) += E_xi_bounceAvg_f1/xi0_f;
                }
            }
        }

        xiAvgTerm_f2 = this->grid->GetRadialGrid()->GetBA_xi_f2(ir);
        for (len_t j = 0; j < np2+1; j++) {
            for (len_t i = 0; i < np1; i++) {
                E_xi_bounceAvg_f2 = Constants::ec * E_term[ir] * xiAvgTerm_f2[j*np1+i] * sqrtB2OverB;
                if (gridtypePXI) {
                        
                    // If hot tail grid, add to diffusion pp component 
                    if ( np2 != 1 ) {
                        xi0_f = mg->GetXi0_f2(i,j);
                        F2(ir, i, j)  += E_xi_bounceAvg_f2 * (1-xi0_f*xi0_f)/(xi0_f*mg->GetP1(i)) ;
                    }
                }
            }
        }


    }
}

/*void ElectricFieldTerm::SetJacobianBlock(const len_t uqtyId, const len_t derivId, Matrix *jac) {
    if (!this->coefficientsShared && id_f == derivId) {
        this->AdvectionTerm::SetJacobianBlock(uqtyId, derivId, jac);
    } else if (derivId == id_Efield) {
        RebuildElectricFieldTerm(t, dt, x, dF1, dF2, nullptr)
        this->AdvectionTerm::SetJacobianBlock(uqtyId, derivId, jac, dF1, dF2);
    }
}*/




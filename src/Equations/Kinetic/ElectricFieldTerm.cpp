/**
 * Implementation of the electric field advection term in the kinetic equation.
 * Note that we may view the below expression as the definition of E_term, where
 * A^p (i.e. F1 for P-Xi grid) is given by {e*E*xi}/m_e c,
 * so that E_term = {E*xi} / {xi} * <B>/sqrt(<B^2>) evaluated in the passing region
 * (since the RHS vanishes in the trapped region, and is independent of xi0,p for passing).
 * It is identical to E_term = <E*B>/sqrt(<B^2>).
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
    real_t xi0, p;
    real_t E_xi_bounceAvg_f1, E_xi_bounceAvg_f2;
    real_t sqrtB2OverB;
    real_t *E_term = x->GetUnknownData(id_Eterm);
    const real_t *xiAvgTerm_f1, *xiAvgTerm_f2;
    real_t E;
    for (len_t ir = 0; ir < nr; ir++) {
        auto *mg = this->grid->GetMomentumGrid(ir);
        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();
        gridtypePXI         = (gridtype == OptionConstants::MOMENTUMGRID_TYPE_PXI);
        gridtypePPARPPERP   = (gridtype == OptionConstants::MOMENTUMGRID_TYPE_PPARPPERP);
        sqrtB2OverB =  sqrt(this->grid->GetRadialGrid()->GetFSA_B2(ir)) / this->grid->GetRadialGrid()->GetFSA_B(ir);
        
        E =  Constants::ec * E_term[ir] /(Constants::me * Constants::c);
         

        xiAvgTerm_f1 = this->grid->GetRadialGrid()->GetBA_xi_f1(ir) ;
        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1+1; i++) {
                E_xi_bounceAvg_f1 = E * xiAvgTerm_f1[j*(np1+1)+i] * sqrtB2OverB; // (e/mc) {E xi}/xi0
                
                if (gridtypePXI) {
                    xi0 = mg->GetP2(j);
                    F1(ir, i, j) += xi0*E_xi_bounceAvg_f1;
                } else if (gridtypePPARPPERP) {
                    F1(ir, i, j) += E_xi_bounceAvg_f1;
                }
            }
        }

        if (gridtypePXI) {        
            xiAvgTerm_f2 = this->grid->GetRadialGrid()->GetBA_xi_f2(ir);
            for (len_t j = 0; j < np2+1; j++) {
                for (len_t i = 0; i < np1; i++) {
                    E_xi_bounceAvg_f2 = E * xiAvgTerm_f2[j*np1+i] * sqrtB2OverB;
                    xi0 = mg->GetP2_f(j);
                    p   = mg->GetP1(i);
                    F2(ir, i, j) += E_xi_bounceAvg_f2 * (1-xi0*xi0)/p ;
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




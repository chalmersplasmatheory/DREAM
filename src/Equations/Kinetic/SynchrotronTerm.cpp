/**
 * Implementation of the synchrotron radiation reaction advection term in the kinetic equation.
 */

#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include "DREAM/Equations/Kinetic/SynchrotronTerm.hpp"
#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"


using namespace DREAM;

/**
 * Constructor.
 */
SynchrotronTerm::SynchrotronTerm(FVM::Grid *g, enum OptionConstants::momentumgrid_type mgtype)
    : FVM::AdvectionTerm(g) {

    SetName("SynchrotronTerm");

    this->gridtype  = mgtype;
}


/**
 * Build the coefficients of this advection (or diffusion) term.
 */
void SynchrotronTerm::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*){
    const len_t nr = grid->GetNr();
    
    bool gridtypePXI, gridtypePPARPPERP;
    real_t xi0, gamma, p;
    real_t Bmin;
    FVM::RadialGrid *rGrid = grid->GetRadialGrid();
    const real_t *BA1_f1, *BA1_f2, *BA2_f1, *BA2_f2;
    for (len_t ir = 0; ir < nr; ir++) {
        auto *mg = grid->GetMomentumGrid(ir);
        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();
        gridtypePXI         = (gridtype == OptionConstants::MOMENTUMGRID_TYPE_PXI);
        gridtypePPARPPERP   = (gridtype == OptionConstants::MOMENTUMGRID_TYPE_PPARPPERP);

        BA1_f1 = grid->GetBA_B3_f1(ir);
        BA1_f2 = grid->GetBA_B3_f2(ir);
        
        BA2_f1 = grid->GetBA_xi2B2_f1(ir);
        BA2_f2 = grid->GetBA_xi2B2_f2(ir);
        
        Bmin = rGrid->GetBmin(ir);
        real_t preFactor = Bmin*Bmin*constPrefactor;
        if (gridtypePXI) {
            if(np2==1){ // pitch averaged p-component
                for(len_t i=0; i<np1+1; i++){
                    F1(ir,i,0) += getf1_PXI(i, 0, mg, BA1_f1, Bmin, ir, rGrid);
                }
                continue;
            }

            for (len_t j = 0; j < np2; j++)
                for (len_t i = 0; i < np1+1; i++) {
                    F1(ir, i, j)  +=  getf1_PXI(i, j, mg, BA1_f1, Bmin);
                }

            for (len_t j = 0; j < np2+1; j++)
                for (len_t i = 0; i < np1; i++) {
                    F2(ir, i, j)  += getf2_PXI(i, j, mg, BA2_f2, Bmin);
                }
        } else if (gridtypePPARPPERP) {
            for (len_t j = 0; j < np2; j++)
                for (len_t i = 0; i < np1+1; i++) {
                    xi0   = mg->GetXi0_f1(i,j);
                    p     = mg->GetP_f1(i,j);
                    gamma = sqrt(1+p*p);

                    F1(ir, i, j)  += -preFactor * (1-xi0*xi0) *( xi0*p*p*BA1_f1[j*(np1+1)+i] - p*xi0/gamma * BA2_f1[j*(np1+1)+i] ); 
                }

            for (len_t j = 0; j < np2+1; j++) 
                for (len_t i = 0; i < np1; i++) {
                    xi0   = mg->GetXi0_f2(i,j);
                    p     = mg->GetP_f2(i,j);
                    gamma = sqrt(1+p*p);

                    F2(ir, i, j)  += -preFactor * sqrt(1-xi0*xi0) *( (1-xi0*xi0)*p*p*BA1_f2[j*np1+i] + xi0*xi0*p/gamma*BA2_f2[j*np1+i] );
                }
        }
    }
}

real_t SynchrotronTerm::getf1_PXI(len_t i, len_t j, FVM::MomentumGrid *mg, const real_t *BA1_f1, real_t Bmin, len_t ir, FVM::RadialGrid *rGrid){
    real_t preFactor = Bmin*Bmin*this->constPrefactor;
    if (mg->GetNp2() == 1){
        real_t p = mg->GetP_f1(i,0);
        real_t gamma = mg->GetGamma_f1(i,0);
        return -2.0/3.0*preFactor*p*gamma*rGrid->GetFSA_B2(ir);
    }

    real_t xi0 = mg->GetP2(j);
    real_t p = mg->GetP1_f(i);

    return -preFactor * p*sqrt(1+p*p)*(1-xi0*xi0) * BA1_f1[j*(mg->GetNp1()+1)+i];
}

real_t SynchrotronTerm::getf2_PXI(len_t i, len_t j, FVM::MomentumGrid *mg, const real_t *BA2_f2, real_t Bmin){
    real_t preFactor = Bmin*Bmin*this->constPrefactor;
    
    real_t xi0 = mg->GetP2_f(j);
    real_t p = mg->GetP1(i);
    real_t gamma = sqrt(1+p*p);

    return +preFactor * (1-xi0*xi0)*xi0/gamma * BA2_f2[j*mg->GetNp1()+i];
}

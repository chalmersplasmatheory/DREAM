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
        
        Bmin = grid->GetRadialGrid()->GetBmin(ir);

        //this->grid->GetRadialGrid()->GetBA_BOverBOverXi_f1(ir); 


        if (gridtypePXI) {
            for (len_t j = 0; j < np2; j++) {
                for (len_t i = 0; i < np1+1; i++) {
                    xi0 = mg->GetP2(j);
                    p = mg->GetP1_f(i);

                    F1(ir, i, j)  += -constPrefactor * p*sqrt(1+p*p)*(1-xi0*xi0) * Bmin*Bmin * BA1_f1[j*(np1+1)+i] ;
                }
            }

            for (len_t j = 0; j < np2+1; j++) {
                for (len_t i = 0; i < np1; i++) {
                    xi0 = mg->GetP2_f(j);
                    p = mg->GetP1(i);
                    gamma = sqrt(1+p*p);

                    F2(ir, i, j)  += +constPrefactor * (1-xi0*xi0)*xi0/gamma * Bmin*Bmin * BA2_f2[j*np1+i] ;
                }
            }
        } else if (gridtypePPARPPERP) {
            for (len_t j = 0; j < np2; j++) {
                for (len_t i = 0; i < np1+1; i++) {
                    xi0   = mg->GetXi0_f1(i,j);
                    p     = mg->GetP_f1(i,j);

                    F1(ir, i, j)  += - constPrefactor * Bmin*Bmin * (1-xi0*xi0) *( xi0*p*p*BA1_f1[j*(np1+1)+i] - p*xi0/gamma * BA2_f1[j*(np1+1)+i] ); 
                }
            }

            for (len_t j = 0; j < np2+1; j++) {
                for (len_t i = 0; i < np1; i++) {
                    xi0   = mg->GetXi0_f2(i,j);
                    p     = mg->GetP_f2(i,j);

                    F2(ir, i, j)  += - constPrefactor * Bmin*Bmin* sqrt(1-xi0*xi0) *( (1-xi0*xi0)*p*p*BA1_f2[j*np1+i] + xi0*xi0*p/gamma*BA2_f2[j*np1+i] );
                }
            }
        }

    }
}

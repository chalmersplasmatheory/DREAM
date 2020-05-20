/**
 * Implementation of the pitch angle scattering term in the 
 * kinetic equation.
 */

#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Equations/Kinetic/PitchScatterTerm.hpp"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"


using namespace DREAM;

/**
 * Constructor.
 */
PitchScatterTerm::PitchScatterTerm(FVM::Grid *g, CollisionQuantityHandler *cqh, EquationSystem *es, enum OptionConstants::momentumgrid_type mgtype)
    : FVM::DiffusionTerm(g) {
        this->gridtype  = mgtype;
        this->nuD       = cqh->GetNuD();
        this->grid      = g;
        this->eqSys     = es;
}


/**
 * Build the coefficients of this advection (or diffusion) term.
 */
void PitchScatterTerm::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler *){
    const len_t nr = this->grid->GetNr();
    bool gridtypePXI, gridtypePPARPPERP;

    real_t xi0_f1, xi0_f2, p_f1, p_f2;
    const real_t *xiBAvg_f1, *xiBAvg_f2;
    real_t *const* nu_D_f1 = nuD->GetValue_f1();
    real_t *const* nu_D_f2 = nuD->GetValue_f2();

    for (len_t ir = 0; ir < nr; ir++) {
        auto *mg = this->grid->GetMomentumGrid(ir);
        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();
        real_t commonFactor_f1, commonFactor_f2;
        gridtypePXI         = (gridtype == OptionConstants::MOMENTUMGRID_TYPE_PXI);
        gridtypePPARPPERP   = (gridtype == OptionConstants::MOMENTUMGRID_TYPE_PPARPPERP);

        // No non-zero elements if np2<2.        
        if ( gridtypePXI && (np2 == 1) ){
            continue;
        }
        
        // Evaluates {xi^2(1-xi^2)Bmin^2/B^2} on flux grid 2
        //xiBAvg_f2 = this->grid->GetRadialGrid()->GetBA_xi21MinusXi2OverB2_f2(ir);
        // Retrieves the average {(Bmin/B)(xi^2/xi0^2)} on p2 flux grid. 
        xiBAvg_f2 = this->grid->GetRadialGrid()->GetBA_xi2OverB_f2(ir);
        
        for (len_t j = 0; j < np2+1; j++) {
            for (len_t i = 0; i < np1; i++) {
                
                commonFactor_f2 = 0.5 * xiBAvg_f2[j*np1+i]*nu_D_f2[ir][j*np1+i];
                xi0_f2 = mg->GetXi0_f2(i,j);
                if (gridtypePXI) {
                    D22(ir,i,j) +=  commonFactor_f2 * (1 - xi0_f2*xi0_f2) ;

                // If ppar-pperp grid
                } else if (gridtypePPARPPERP) {
                    p_f2 = mg->GetP_f2(i,j);
                    
                    D22(ir,i,j) +=  commonFactor_f2 * p_f2*p_f2 * xi0_f2*xi0_f2;
                    D21(ir,i,j) += -commonFactor_f2 * p_f2*p_f2 * xi0_f2*sqrt(1-xi0_f2*xi0_f2);
                    
                }
            }
        }

        
        if (gridtypePPARPPERP) {
            // Evaluates {xi^2(1-xi^2)Bmin^2/B^2} on flux grid 1
            //xiBAvg_f1 = this->grid->GetRadialGrid()->GetBA_xi21MinusXi2OverB2_f1(ir);
            // Retrieves the average {(Bmin/B)(xi^2/xi0^2)} on p2 flux grid. 
            xiBAvg_f1 = this->grid->GetRadialGrid()->GetBA_xi2OverB_f1(ir);
            for (len_t j = 0; j < np2; j++) {
                for (len_t i = 0; i < np1+1; i++) {
                    p_f1 = mg->GetP_f1(i,j);
                    commonFactor_f1 = 0.5 * nu_D_f1[ir][j*(np1+1)+i] *xiBAvg_f1[j*(np1+1)+i];
                    xi0_f1 = mg->GetXi0_f1(i,j);
                    D11(ir,i,j) +=  commonFactor_f1 * p_f1*p_f1 * (1 - xi0_f1*xi0_f1);
                    D12(ir,i,j) += -commonFactor_f1 * p_f1*p_f1 * xi0_f1*sqrt(1-xi0_f1*xi0_f1)  ;
                }
            }
        }

    }
}






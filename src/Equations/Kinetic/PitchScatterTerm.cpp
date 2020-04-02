

/**
 * Implementation of the electric field term in the kinetic equation.
 * If using a hot tail grid it is added as a diffusion operator,
 * otherwise as the regular advection term.
 * 
 * TODOs: Need a couple of new functions:
 * real_t GetB2Avg(ir): calculates <B^2>
 * real_t GetBOverXiAvg_f: calculates <B/xi>, where
 *      xi=xi(B,xi0) in p-xi or xi=xi(B,ppar0,pperp0) in ppar-pperp.
 *      I guess that this doesn't have to be stored, and could then exist 
 *      just as a function of xi0 (=ppar0/(ppar0^2+pperp0^2) in RadialGrid 
 *      (implemented by the radial grid generators).
 * bool IsTrapped_{f1/f2/fr}(ir,i,j): true if phase space 
 *                       coordinate denotes trapped orbit
 */

#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Equations/CollisionFrequencyCreator.hpp"
#include "DREAM/Equations/Kinetic/PitchScatterTerm.hpp"


using namespace DREAM;

/**
 * Constructor.
 */
PitchScatterTerm::PitchScatterTerm(FVM::Grid *g, CollisionFrequencyCreator *cfc, EquationSystem *es, enum SimulationGenerator::momentumgrid_type mgtype)
    : FVM::DiffusionTerm(g) {
        this->gridtype  = mgtype;
        this->collFreqs = cfc;
        this->eqSys     = es;
        this->grid      = g;
}


/**
 * Build the coefficients of this advection (or diffusion) term.
 */
void PitchScatterTerm::Rebuild(){
    const len_t nr = this->grid->GetNr();
    bool gridtypePXI, gridtypePPARPPERP;

    real_t xi0_f1, xi0_f2, p_f1, p_f2;
    const real_t *xiBAvg_f1, *xiBAvg_f2;

    for (len_t ir = 0; ir < nr; ir++) {
        auto *mg = this->grid->GetMomentumGrid(ir);
        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();
        real_t commonFactor_f1, commonFactor_f2;
        gridtypePXI         = (gridtype == SimulationGenerator::MOMENTUMGRID_TYPE_PXI);
        gridtypePPARPPERP   = (gridtype == SimulationGenerator::MOMENTUMGRID_TYPE_PPARPPERP);

        // Skip to next radial grid point if grid is hot-tail.        
        if ( !( gridtypePXI && (np2 == 1) ) ){
            continue;
        }
        
        // Evaluates {xi^2(1-xi^2)Bmin^2/B^2}
        xiBAvg_f1 = this->grid->GetRadialGrid()->GetXi21MinusXi2OverB2Avg_f1(ir);
        xiBAvg_f2 = this->grid->GetRadialGrid()->GetXi21MinusXi2OverB2Avg_f2(ir);
        
        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1; i++) {
                
                commonFactor_f1 = 0.5 * collFreqs->nu_D(ir,p_f1)*xiBAvg_f1[j*np1+i];
                commonFactor_f2 = 0.5 * collFreqs->nu_D(ir,p_f2)*xiBAvg_f2[j*np1+i];
                
                if (gridtypePXI) {
                    p_f1   = mg->GetP1_f(i);
                    p_f2   = mg->GetP1(i);
                    xi0_f1 = mg->GetP2(j);
                    xi0_f2 = mg->GetP2_f(j);

                    D22(ir,i,j) += commonFactor_f2 / (xi0_f2*xi0_f2) ;
                    

                // If ppar-pperp grid
                } else if (gridtypePPARPPERP) {
                    p_f1 = sqrt(mg->GetP1_f(i)*mg->GetP1_f(i) + mg->GetP2(j)*mg->GetP2(j));
                    p_f2 = sqrt(mg->GetP1(i)*mg->GetP1(i) + mg->GetP2_f(j)*mg->GetP2_f(j));
                    xi0_f1 = mg->GetXi0_f1(i,j);
                    xi0_f2 = mg->GetXi0_f2(i,j);
                    D11(ir,i,j) += commonFactor_f1 * p_f1*p_f1 / (xi0_f1*xi0_f1);
                    D22(ir,i,j) += commonFactor_f2 * p_f2*p_f2 / (1- xi0_f2*xi0_f2);
                    D12(ir,i,j) += -commonFactor_f1 * p_f1*p_f1 /( xi0_f1*sqrt(1-xi0_f1*xi0_f1) ) ;
                    D21(ir,i,j) += -commonFactor_f2 * p_f2*p_f2 /( xi0_f2*sqrt(1-xi0_f2*xi0_f2) );
                    
                }
            }
        }
    }
}






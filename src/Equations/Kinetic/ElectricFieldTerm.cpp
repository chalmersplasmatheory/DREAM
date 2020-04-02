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

#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Equations/CollisionFrequencyCreator.hpp"
#include "DREAM/Equations/Kinetic/ElectricFieldTerm.hpp"


using namespace DREAM;

/**
 * Constructor.
 */
ElectricFieldTerm::ElectricFieldTerm(FVM::Grid *g, CollisionFrequencyCreator *cfc, EquationSystem *es, enum SimulationGenerator::momentumgrid_type mgtype)
    : FVM::AdvectionTerm(g), FVM::DiffusionTerm(g) {
        this->gridtype  = mgtype;
        this->collFreqs = cfc;
        this->eqSys     = es;
        this->grid      = g;
}


/**
 * Build the coefficients of this advection (or diffusion) term.
 */
void ElectricFieldTerm::Rebuild(){
    const len_t id_Eterm = this->eqSys->GetUnknownID("E_field"); // E term should be <E*B>/sqrt(<B^2>)
    const len_t nr = this->grid->GetNr();
    bool gridtypePXI, gridtypePPARPPERP;
    real_t p, p_f, xi0_f, p_f1, p_f2;
    real_t E_xi_bounceAvg_f1, E_xi_bounceAvg_f2;
    real_t *E_term = this->eqSys->GetUnknownData(id_Eterm);
    const real_t *xiAvg_f1, *xiAvg_f2;

    for (len_t ir = 0; ir < nr; ir++) {
        auto *mg = this->grid->GetMomentumGrid(ir);
        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();
        gridtypePXI         = (gridtype == SimulationGenerator::MOMENTUMGRID_TYPE_PXI);
        gridtypePPARPPERP   = (gridtype == SimulationGenerator::MOMENTUMGRID_TYPE_PPARPPERP);
        
        xiAvg_f1 = this->grid->GetRadialGrid()->GetXiAvg_f1(ir);
        xiAvg_f2 = this->grid->GetRadialGrid()->GetXiAvg_f2(ir);
        

        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1; i++) {
                E_xi_bounceAvg_f1 = Constants::ec * E_term[ir] * xiAvg_f1[j*np1+i];
                E_xi_bounceAvg_f2 = Constants::ec * E_term[ir] * xiAvg_f2[j*np1+i];
                
                if (gridtypePXI) {
                        
                    // If hot tail grid, add to diffusion pp component 
                    if (np2 == 1) {
                        D11(ir, i, j) +=  1/3
                            * this->grid->GetRadialGrid()->GetEffPassFrac(ir) 
                            * Constants::ec * Constants::ec
                            * E_term[ir]*E_term[ir]
                            / this->collFreqs->nu_D(ir,p_f1);

                    // If runaway grid, add to advection
                    } else {
                        p     = mg->GetP1(i);
                        xi0_f = mg->GetP2_f(j);
                        F1(ir, i, j)  += E_xi_bounceAvg_f1;
                        F2(ir, i, j)  += E_xi_bounceAvg_f2 * (1-xi0_f*xi0_f)/(p*xi0_f) ;
                    }

                // If ppar-pperp grid
                } else if (gridtypePPARPPERP) {
                    xi0_f = mg->GetP1_f(i) / sqrt(mg->GetP1_f(i)*mg->GetP1_f(i) + mg->GetP2(i)*mg->GetP2(i));
                    
                    F1(ir, i, j) += E_xi_bounceAvg_f1/xi0_f;
                }
            }
        }
    }
}






/**
 * Implementation of the electric field term in the kinetic equation.
 * If using a hot tail grid it is added as a diffusion operator,
 * otherwise as the regular advection term.
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
    real_t xi0_f;
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
                            / this->collFreqs->nu_D(ir,mg->GetP1_f(i));

                    // If runaway p-xi grid, add to advection 
                    } else {
                        xi0_f = mg->GetP2_f(j);
                        F1(ir, i, j)  += E_xi_bounceAvg_f1;
                        F2(ir, i, j)  += E_xi_bounceAvg_f2 * (1-xi0_f*xi0_f)/(xi0_f*mg->GetP1(i)) ;
                    }

                // If ppar-pperp grid
                } else if (gridtypePPARPPERP) {
                    xi0_f = mg->GetXi0_f1(i,j);
                    F1(ir, i, j) += E_xi_bounceAvg_f1/xi0_f;
                }
            }
        }
    }
}






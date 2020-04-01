/**
 * Implementation of the p*nu_s friction term in the kinetic equation.
 */

#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Equations/CollisionFrequencyCreator.hpp"
#include "DREAM/Equations/Kinetic/SlowingDownTerm.hpp"


using namespace DREAM;

/**
 * Constructor.
 */
SlowingDownTerm::SlowingDownTerm(FVM::Grid *g, CollisionFrequencyCreator *cfc, enum SimulationGenerator::momentumgrid_type mgtype)
    : FVM::AdvectionTerm(g) {
        this->gridtype = mgtype;
        this->collFreqs = cfc;
}

/**
 * Build the coefficients of this advection term.
 */
void SlowingDownTerm::Rebuild(){
    const len_t nr = this->grid->GetNr();
    real_t p, p_f1, p_f2;
    for (len_t ir = 0; ir < nr; ir++) {
        auto *mg = this->grid->GetMomentumGrid(ir);
        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();


        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1; i++) {
                
                if (gridtype == SimulationGenerator::MOMENTUMGRID_TYPE_PXI) {
                    p = mg->GetP1_f(i);
                    F1(ir, i, j) += p * collFreqs->nu_s(ir,p);
                } else if (gridtype == SimulationGenerator::MOMENTUMGRID_TYPE_PPARPPERP) {
                    p_f1 = sqrt(mg->GetP1_f(i)*mg->GetP1_f(i) + mg->GetP2(j)*mg->GetP2(j));
                    p_f2 = sqrt(mg->GetP1(i)*mg->GetP1(i) + mg->GetP2_f(j)*mg->GetP2_f(j));
                    F1(ir, i, j) += mg->GetP1_f(i) * collFreqs->nu_s(ir,p_f1);
                    F2(ir, i, j) += mg->GetP2_f(j) * collFreqs->nu_s(ir,p_f2);
                }
            }
        }
    }
}



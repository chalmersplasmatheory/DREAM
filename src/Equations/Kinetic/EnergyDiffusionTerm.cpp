/**
 * Implementation of the p*nu_s friction term in the kinetic equation.
 */

#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Equations/CollisionFrequencyCreator.hpp"
#include "DREAM/Equations/Kinetic/EnergyDiffusionTerm.hpp"


using namespace DREAM;

/**
 * Constructor.
 */
EnergyDiffusionTerm::EnergyDiffusionTerm(FVM::Grid *g, CollisionQuantityHandler *cqh, EquationSystem* es, enum SimulationGenerator::momentumgrid_type mgtype)
    : FVM::DiffusionTerm(g) {
        this->gridtype = mgtype;
        this->collQty  = cqh;
        this->eqSys    = es;
}

/**
 * Build the coefficients of this advection term.
 */
void EnergyDiffusionTerm::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler *x){
    const len_t nr = this->grid->GetNr();
 
    real_t *const* nu_par_f1 = collQty->GetNuPar_f1();
    real_t *const* nu_par_f2 = collQty->GetNuPar_f2();
  
    bool gridtypePXI, gridtypePPARPPERP;


    for (len_t ir = 0; ir < nr; ir++) {
        auto *mg = this->grid->GetMomentumGrid(ir);
        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();

        real_t xi0_f1, xi0_f2;
        gridtypePXI        = (gridtype == SimulationGenerator::MOMENTUMGRID_TYPE_PXI);
        gridtypePPARPPERP  = (gridtype == SimulationGenerator::MOMENTUMGRID_TYPE_PPARPPERP);
        
        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1+1; i++) {
                if (gridtypePXI)
                    D11(ir, i, j) += nu_par_f1[ir][j*(np1+1)+i];
                else if (gridtypePPARPPERP){
                    xi0_f1 = mg->GetXi0_f1(i,j);
                    D11(ir, i, j) += xi0_f1*xi0_f1*nu_par_f1[ir][j*(np1+1)+i];
                    D12(ir, i,j ) += xi0_f1*sqrt(1-xi0_f1*xi0_f1)*nu_par_f1[ir][j*(np1+1)+i];
                } 
            }
        }

        if (gridtypePPARPPERP) {
            for (len_t j = 0; j < np2+1; j++) {
                for (len_t i = 0; i < np1; i++) {
                    xi0_f2 = mg->GetXi0_f2(i,j);
                    D22(ir, i, j) += (1-xi0_f2*xi0_f2) * nu_par_f2[ir][j*np1+i];
                    D21(ir, i, j) += xi0_f2*sqrt(1-xi0_f2*xi0_f2) * nu_par_f2[ir][j*np1+i];
                }
            }
        }
    
    }
}



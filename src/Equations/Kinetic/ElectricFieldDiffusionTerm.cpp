/**
 * Implementation of the electric field diffusion term in the kinetic equation,
 * which is used for hot-tail like grids.
 */

#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include "DREAM/Equations/Kinetic/ElectricFieldDiffusionTerm.hpp"


using namespace DREAM;

/**
 * Constructor.
 */
ElectricFieldDiffusionTerm::ElectricFieldDiffusionTerm(FVM::Grid *g, CollisionQuantityHandler *cqh, EquationSystem *es, enum SimulationGenerator::momentumgrid_type mgtype)
    : FVM::DiffusionTerm(g) {
        this->gridtype  = mgtype;
        this->collQty   = cqh;
        this->eqSys     = es;
        this->grid      = g;
        this->id_Eterm  = this->eqSys->GetUnknownID( SimulationGenerator::UQTY_E_FIELD ); // E term should be <E*B>/sqrt(<B^2>)
    
}


/**
 * Build the coefficients of this diffusion term. Realistically only used when np2 = 1, but let's keep it general.
 */
void ElectricFieldDiffusionTerm::Rebuild(
    const real_t, const real_t, FVM::UnknownQuantityHandler *x
){
    const len_t nr = this->grid->GetNr();
    real_t *E_term = x->GetUnknownData(id_Eterm);
    real_t *const *nu_D_f1 = collQty->GetNuD_f1();
    for (len_t ir = 0; ir < nr; ir++) {
        auto *mg = this->grid->GetMomentumGrid(ir);
        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();        

        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1+1; i++) {
                D11(ir, i, j) +=  1/3
                    * this->grid->GetRadialGrid()->GetEffPassFrac(ir) 
                    * Constants::ec * Constants::ec
                    * E_term[ir] * E_term[ir] / nu_D_f1[ir][j*(np1+1)+i];  
            }
        }
    }
}






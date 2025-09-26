#include "DREAM/Equations/Fluid/NBIElectronTerm.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"


using namespace DREAM;

/**
 * Constructor
 */
NBIElectronTerm::NBIElectronTerm(NBIHandler *h, FVM::Grid *grid, FVM::UnknownQuantityHandler *unknowns)
: EquationTerm(grid), handler(h), radialGrid(grid->GetRadialGrid()) {

    this->id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    this->id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->id_ion_density = unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS);
    this->id_ion_temperature = unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER);
}

/**
 * Rebuild the term (called once per time step). Call of the build function in NBIHandler
 */
void NBIElectronTerm::Rebuild(const real_t t, const real_t dt, FVM::UnknownQuantityHandler *unknowns){
    handler->Build(t, dt, unknowns);
    this->nr = radialGrid->GetNr();
}

/**
 * Set the vector elements corresponding to this term.
 */
void NBIElectronTerm::SetVectorElements(real_t *rhs, const real_t*){
    const real_t *Qe = handler->GetNBIHeatTerm_e();
    for (len_t ir=0; ir<nr; ++ir) {
        rhs[ir] -= Qe[ir];
    }
}

/**
 * Set the matrix elements corresponding to this term.
 */
void NBIElectronTerm::SetMatrixElements(FVM::Matrix*, real_t *rhs){ //TODO
    const real_t *Qe = handler->GetNBIHeatTerm_e();
    const real_t *nc = unknowns->GetUnknownData(id_ncold);

    for (len_t ir=0; ir<nr; ++ir){
        const real_t factor = 2.0/(3.0*nc[ir]);
        rhs[ir] -= factor * Qe[ir];
    }
}

/**
 * Set the Jacobian elements corresponding to this term.
 */
bool NBIElectronTerm::SetJacobianBlock(const len_t /*uqtyId*/, const len_t derivId,
                                       FVM::Matrix *jac, const real_t*){
    const real_t *H_r_dTe = handler->GetH_r_dTe();
    const real_t *H_r_dni = handler->GetH_r_dni();
    const real_t *H_r_dTi = handler->GetH_r_dTi();
    const real_t *H_r_dne = handler->GetH_r_dne();
    if (derivId != id_ncold && derivId != id_Tcold &&
        derivId != id_ion_density && derivId != id_ion_temperature)
        return false; 

    for (len_t ir = 0; ir < nr; ++ir){
        real_t dP = 0.0;

        if (derivId == id_ncold){
            dP = H_r_dne[ir];
        }
        else if (derivId == id_Tcold){
            dP = H_r_dTe[ir];
        }
        else if (derivId == id_ion_density){
            dP = H_r_dni[ir];
        }
        else if (derivId == id_ion_temperature){
            dP = H_r_dTi[ir];
        }

        jac->SetElement(ir, ir, dP);
    }

    return false;
}



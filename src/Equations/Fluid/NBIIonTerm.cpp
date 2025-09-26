#include "DREAM/Equations/Fluid/NBIIonTerm.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"


using namespace DREAM;

/**
 * Constructor
 */
NBIIonTerm::NBIIonTerm(NBIHandler *h, FVM::Grid *grid, IonHandler *ions, FVM::UnknownQuantityHandler *unknowns, len_t iz)
        : EquationTerm(grid), handler(h),
          radialGrid(grid->GetRadialGrid()), ions(ions) {

    this->id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    this->id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->id_ion_density = unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS);
    this->id_ion_temperature = unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER);
    this->iz = iz;
}


/**
 * Rebuild the term (called once per time step). Call of the build function in NBIHandler
 */
void NBIIonTerm::Rebuild(const real_t t, const real_t dt, FVM::UnknownQuantityHandler *unknowns){
    handler->Build(t, dt, unknowns);
    this->nr = radialGrid->GetNr();
}

/**
 * Set the vector elements corresponding to this term. Also calculates the weight factor for each species
 */
void NBIIonTerm::SetVectorElements(real_t *rhs, const real_t*){
    const real_t *Qi = handler->GetNBIHeatTerm_i();

    // Calculate the weight factor for each species
    len_t NZ = ions->GetNZ(); //Nr of ions species

    std::vector<real_t> ni_species(NZ, 0.0);   
    std::vector<real_t> Zi(NZ, 0.0);           
    std::vector<real_t> Mi(NZ, 0.0);  
    

    //Calculate the denominator, needed for all species and flux surfaces
    for (len_t ir=0; ir<nr; ++ir){
        real_t denom = 0;
        for (len_t iz_denom = 0; iz_denom < ions->GetNZ(); iz_denom++){
                len_t Zmax = ions->GetZ(iz_denom);
                Zi[iz_denom] = Zmax;
                Mi[iz_denom] = ions->GetIonSpeciesMass(iz_denom); //approx

                // Loop over charge states of this species and sum the density
                real_t sum_ni = 0.0;
                for (len_t Z0 = 0; Z0 <= Zmax; Z0++){
                    real_t niZ = ions->GetIonDensity(ir, iz_denom, Z0);
                    sum_ni += niZ; 
                }
                ni_species[iz_denom] = sum_ni;  
                denom += ni_species[iz_denom] * Zi[iz_denom] / Mi[iz_denom];
        }
        //Calculate the weight factor for this species and add to rhs
        len_t Zmax = ions->GetZ(iz);
        Zi[iz] = Zmax;
        Mi[iz] = ions->GetIonSpeciesMass(iz); 
        real_t sum_ni = 0.0;
        for (len_t Z0 = 0; Z0 <= Zmax; Z0++){
            real_t niZ = ions->GetIonDensity(ir, iz, Z0);
                sum_ni += niZ; 
        }
            ni_species[iz] = sum_ni;  
            real_t w = ni_species[iz]* Zi[iz]/Mi[iz]/denom;
            rhs[ir] -= w * Qi[ir]; 
    }
}
    


/**
 * Set the matrix elements corresponding to this term.
 */
void NBIIonTerm::SetMatrixElements(FVM::Matrix*, real_t *rhs){ //TODO
    
    
}
/**
 * Set the Jacobian elements corresponding to this term.
 */

bool NBIIonTerm::SetJacobianBlock(const len_t /*uqtyId*/, const len_t derivId, FVM::Matrix *jac, const real_t*){
    
    const real_t *H_r_dTe = handler->GetH_r_dTe();
    const real_t *H_r_dni = handler->GetH_r_dni();
    const real_t *H_r_dTi = handler->GetH_r_dTi();
    const real_t *H_r_dne = handler->GetH_r_dne();
    if (derivId != id_ncold && derivId != id_Tcold &&
        derivId != id_ion_density && derivId != id_ion_temperature)
        return false; 
        // Set to 0 for now
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

        jac->SetElement(ir, ir, 0);
    }
    return false;
}



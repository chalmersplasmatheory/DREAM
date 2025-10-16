#include "DREAM/Equations/Fluid/NBIElectronTerm.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"


using namespace DREAM;

/**
 * Constructor
 */
NBIElectronTerm::NBIElectronTerm(NBIHandler *h, FVM::Grid *grid, FVM::UnknownQuantityHandler *unknowns, IonHandler* ionHandler)
: EquationTerm(grid), handler(h), ions(ionHandler), radialGrid(grid->GetRadialGrid()) {
    this->unknowns = unknowns;
    this->nr = radialGrid->GetNr();
    this->id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    this->id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->id_ion_density = unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS);
    this->id_ion_temperature = unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER);
    this->id_ni = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    this-> Qe_1 = Qe_1;
    this-> Qe_2 = Qe_2;
    this-> Qe_3 = Qe_3;
}

/**
 * Rebuild the term (called once per time step). Call of the build function in NBIHandler
 */
void NBIElectronTerm::Rebuild(const real_t t, const real_t dt, FVM::UnknownQuantityHandler *unknowns){
   
    real_t beam_energy = handler->GetBeamEnergy();  
    //printf("Beam energy in NBIElectronTerm: %e J\n", beam_energy);
    handler->Build(t, dt, unknowns, beam_energy);
    Qe_1 = handler->GetNBIHeatTerm_e();
    //handler->Build(t, dt, unknowns, beam_energy/2);
    //printf("Beam energy/2 in NBIElectronTerm: %e J\n", beam_energy/2);
    //Qe_2 = handler->GetNBIHeatTerm_e();
    //handler->Build(t, dt, unknowns, beam_energy/3);
    //Qe_3 = handler->GetNBIHeatTerm_e();
}

/**
 * Set the vector elements corresponding to this term.
 */
void NBIElectronTerm::SetVectorElements(real_t *rhs, const real_t*){
    for (len_t ir=0; ir<nr; ++ir) {
        rhs[ir] -= this->Qe_1[ir]; // + this->Qe_2[ir] + this->Qe_3[ir];
        //printf("NBI heating term at ir = %d: %e W/m^3\n", ir, Qe[ir]);
    }
    
}

/**
 * Set the matrix elements corresponding to this term.
 */
void NBIElectronTerm::SetMatrixElements(FVM::Matrix*, real_t *rhs){ //TODO
    //const real_t *Qe = handler->GetNBIHeatTerm_e();
    //const real_t *nc = unknowns->GetUnknownData(id_ncold);

    //for (len_t ir=0; ir<nr; ++ir){
     //   const real_t factor = 2.0/(3.0*nc[ir]);
     //   rhs[ir] -= factor * Qe[ir];
    //} TODO
}

/**
 * Set the Jacobian elements corresponding to this term.
 */
bool NBIElectronTerm::SetJacobianBlock(const len_t uqtyId, const len_t derivId,
                                       FVM::Matrix *jac, const real_t*){
    // Check if this derivId is one we handle
    if (derivId != id_ncold && derivId != id_Tcold &&
        derivId != id_ni && derivId != id_ion_temperature)
        return false; 
    
    // Get derivative arrays from NBIHandler
    const real_t *d_NBIHeatTerm_e_d_Te = handler->Getd_NBIHeatTerm_e_d_Te();
    const real_t *d_NBIHeatTerm_e_d_ne = handler->Getd_NBIHeatTerm_e_d_ne();
    const std::vector<std::vector<std::vector<real_t>>>& d_NBIHeatTerm_e_d_n_ij = handler->Getd_NBIHeatTerm_e_d_n_ij();
    const std::vector<std::vector<std::vector<real_t>>>& d_NBIHeatTerm_e_d_T_ij = handler->Getd_NBIHeatTerm_e_d_T_ij();

    // Set Jacobian elements for each radial point
    if (derivId == id_ncold) {
        for (len_t ir = 0; ir < nr; ++ir) {
            jac->SetElement(ir, ir, 0);//d_NBIHeatTerm_e_d_ne[ir]);
        }
    }
    
    if (derivId == id_Tcold) {
        for (len_t ir = 0; ir < nr; ++ir) {
            jac->SetElement(ir, ir, 0);//d_NBIHeatTerm_e_d_Te[ir]);
        }
    }
    
    if (derivId == id_ni) {
        // here we set per species and chanrge
        for (len_t ir = 0; ir < nr; ++ir) {
            for (len_t iIon = 0; iIon < ions->GetNZ(); ++iIon) {
                const len_t Zmax = ions->GetZ(iIon);
                const len_t speciesStart = ions->GetIndex(iIon, 0) * nr;
                for (len_t Zp = 0; Zp <= Zmax; ++Zp) {
                    real_t dP = d_NBIHeatTerm_e_d_n_ij[ir][iIon][Zp];
                    len_t col = speciesStart + Zp * nr + ir;
                    jac->SetElement(ir, col, 0);//dP); //why not col, col
                }
            }
        }
    }
    
    if (derivId == id_ion_temperature) {
        // Derivative w.r.t. ion temperature: loop over all species
        for (len_t ir = 0; ir < nr; ++ir) {
            for (len_t iIon = 0; iIon < ions->GetNZ(); ++iIon) {
                const len_t Zmax = ions->GetZ(iIon);
                
                // Get total ion density for this species
                real_t ni_species = 0.0;
                for (len_t Zp = 0; Zp <= Zmax; ++Zp) {
                    const len_t ni_idx = ions->GetIndex(iIon, Zp) * nr + ir;
                    ni_species += unknowns->GetUnknownData(id_ion_density)[ni_idx];
                }
                
                // Sum derivatives over all charge states
                real_t dQe_dTi_species = 0.0;
                for (len_t Zp = 0; Zp <= Zmax; ++Zp) {  
                    dQe_dTi_species += d_NBIHeatTerm_e_d_T_ij[ir][iIon][Zp];
                }
                
                // Convert from dQ/dTi to dQ/dWi
                // Wi = (3/2) * ni * Ti, so dTi/dWi = 1/(3/2 * ni)
                const real_t eq = 1.602e-19; // eV to Joules
                real_t dQe_dWi = 0.0;
                if (ni_species > 0) {
                    dQe_dWi = dQe_dTi_species / (1.5 * ni_species * eq);
                }
                
                // WI_ENER has one value per species, not per charge state
                len_t col = iIon * nr + ir;
                jac->SetElement(ir, col, 0);//dQe_dWi); 
            }
        }
    }

    return true;
}



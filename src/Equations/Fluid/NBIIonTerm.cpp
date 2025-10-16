#include "DREAM/Equations/Fluid/NBIIonTerm.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"


using namespace DREAM;

/**
 * Constructor
 */
NBIIonTerm::NBIIonTerm(
    NBIHandler* h, FVM::Grid* grid, IonHandler* ionHandler,
    FVM::UnknownQuantityHandler* unknowns, const len_t iIon
)
: IonEquationTerm<FVM::EquationTerm>(grid, ionHandler, iIon)
, handler(h)
, ions(ionHandler)
, radialGrid(grid->GetRadialGrid()) {
    this->unknowns = unknowns;
    this->nr = radialGrid->GetNr();

    this->id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    this->id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->id_ion_density = unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS);
    this->id_ion_temperature = unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER);
    this->id_ni = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    this->iIon = iIon;
}


/**
 * Rebuild the term (called once per time step). Call of the build function in NBIHandler
 */
void NBIIonTerm::Rebuild(const real_t t, const real_t dt, FVM::UnknownQuantityHandler *unknowns){

    real_t beam_energy = handler->GetBeamEnergy();  
    //printf("Beam energy in NBIElectronTerm: %e J\n", beam_energy);
    handler->Build(t, dt, unknowns, beam_energy);
}

/**
 * Set the vector elements corresponding to this term. Also calculates the weight factor for each species. Called for every ion species iIon and charge state Z0.
 */
void NBIIonTerm::SetCSVectorElements(
            real_t *vec, const real_t *x,
        const len_t iIon, const len_t Z0, const len_t rOffset){

    
    const real_t *Qi = handler->GetNBIHeatTerm_i();

    // Calculate the weight factor for each species
    len_t NZ = ions->GetNZ(); //Nr of ions species

    std::vector<real_t> ni_species(NZ, 0.0);   
    std::vector<real_t> Zi(NZ, 0.0);           
    std::vector<real_t> Mi(NZ, 0.0);  
    

    for (len_t ir=0; ir<nr; ++ir){
        real_t denom = 0;
        //Calculate the weighing factor
        for (len_t iz_denom = 0; iz_denom < ions->GetNZ(); iz_denom++){
                len_t Zmax = ions->GetZ(iz_denom);
                Zi[iz_denom] = Zmax;
                Mi[iz_denom] = ions->GetIonSpeciesMass(iz_denom); //approx

                // Loop over charge states of this species and sum the density
                real_t sum_ni = 0.0;
                for (len_t Z0_inner = 0; Z0_inner <= Zmax; Z0_inner++){
                    real_t niZ = ions->GetIonDensity(ir, iz_denom, Z0_inner);
                    sum_ni += niZ; 
                }
                ni_species[iz_denom] = sum_ni;  
                denom += ni_species[iz_denom] * Zi[iz_denom] / Mi[iz_denom];
        }

        //Calculate the weight factor for this species
        // Only add once per species (when Z0 == 0)
        if (Z0 == 0) {
            real_t MiIon = ions->GetIonSpeciesMass(iIon);
            real_t zIon = ions->GetZ(iIon);
            real_t w = ni_species[iIon] * zIon / MiIon / denom;
            len_t idx = iIon * nr + ir;
            vec[idx] -= w * Qi[ir]; 
            //printf("NBI heating term at ir = %d, iIon = %d: w=%e, Qi=%e W/m^3\n", ir, iIon, w, Qi[ir]);
        }
    }
}
    


/**
 * Set the matrix elements corresponding to this term.
 */
void NBIIonTerm::SetCSMatrixElements(
            FVM::Matrix *mat, real_t *rhs,
        const len_t iIon, const len_t Z0, const len_t rOffset){ 
     const real_t *Qi = handler->GetNBIHeatTerm_i();
    //TODO
    
}
/**
 * Set the Jacobian elements corresponding to this term.
 */

bool NBIIonTerm::SetCSJacobianBlock(
           const len_t uqtyId, const len_t derivId,
                        FVM::Matrix *jac, const real_t *x,
                        const len_t iIon, const len_t Z0, const len_t rOffset){
    
    // Check if this derivId is one we handle
    if (derivId != id_ncold && derivId != id_Tcold &&
        derivId != id_ni && derivId != id_ion_temperature) {
        return false;
    }
    
    // For ncold, Tcold, and ion_temperature: only set once per species (not per charge state)
    // For ni: set for each charge state since it's the derivative w.r.t. that specific charge state
    if (derivId != id_ni && Z0 != 0) {
        return true;
    }
    
    // Get derivative arrays from NBIHandler
    const real_t *d_NBIHeatTerm_i_d_Te = handler->Getd_NBIHeatTerm_i_d_Te();
    const real_t *d_NBIHeatTerm_i_d_ne = handler->Getd_NBIHeatTerm_i_d_ne();
    const std::vector<std::vector<std::vector<real_t>>>& d_NBIHeatTerm_i_d_n_ij = handler->Getd_NBIHeatTerm_i_d_n_ij();
    const std::vector<std::vector<std::vector<real_t>>>& d_NBIHeatTerm_i_d_T_ij = handler->Getd_NBIHeatTerm_i_d_T_ij();
    
    // Set Jacobian elements for each radial point
    if (derivId == id_ncold) {
        for (len_t ir = 0; ir < nr; ++ir) {
            len_t row = iIon*nr+ir;
            len_t col = ir;
            jac->SetElement(row, col, 0);// d_NBIHeatTerm_i_d_ne[ir]); //why not col, col
        }
    }
    
    if (derivId == id_Tcold) {
        for (len_t ir = 0; ir < nr; ++ir) {
            len_t row = iIon*nr+ir;
            len_t col = ir;
            jac->SetElement(row, col, 0);//d_NBIHeatTerm_i_d_Te[ir]);
        }
    }
    
    if (derivId == id_ni) {
        for (len_t ir = 0; ir < nr; ++ir) {
            len_t row = iIon * nr + ir;
            len_t col = ions->GetIndex(iIon, Z0) * nr + ir;
            jac->SetElement(row, col, 0);//d_NBIHeatTerm_i_d_n_ij[ir][iIon][Z0]);
        }
    }
    
    if (derivId == id_ion_temperature) {
        for (len_t ir = 0; ir < nr; ++ir) {
            
            
            // Get total ion density for this species at this radial position
            const len_t Zmax = ions->GetZ(iIon);
            real_t ni_species = 0.0;
            real_t dQi_dTi_total = 0.0;
            for (len_t Zp = 0; Zp <= Zmax; ++Zp) {
                const len_t ni_idx = ions->GetIndex(iIon, Zp) * nr + ir;
                ni_species += unknowns->GetUnknownData(id_ion_density)[ni_idx];
                dQi_dTi_total += d_NBIHeatTerm_i_d_T_ij[ir][iIon][Zp];
            }

            
            
            // Convert from derivative w.r.t. Ti to derivative w.r.t. Wi
            // Wi = (3/2) * ni * Ti, so dTi/dWi = 1/(3/2 * ni)
            const real_t eq = 1.602e-19; // eV to Joules
            real_t dP = 0.0;
            if (ni_species > 0) {
                dP = dQi_dTi_total / (1.5 * ni_species * eq);
            }
            
            // WI_ENER has one value per species, not per charge state
            len_t col = iIon*nr+ir;
            len_t row = iIon*nr+ir;
            jac->SetElement(row, col, 0);//dP);
        }
    }
    
    return true;
}



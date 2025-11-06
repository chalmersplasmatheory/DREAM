#include "DREAM/Equations/Fluid/NBIIonTerm.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/DREAMException.hpp"



using namespace DREAM;

/**
 * Constructor
 */
NBIIonTerm::NBIIonTerm(
    NBIHandler* h, FVM::Grid* grid, IonHandler* ionHandler,
    FVM::UnknownQuantityHandler* unknowns, const len_t iIon
) : IonEquationTerm<FVM::EquationTerm>(grid, ionHandler, iIon), handler(h), ions(ionHandler), radialGrid(grid->GetRadialGrid()) {
    this->unknowns = unknowns;
    this->nr = radialGrid->GetNr();

    this->id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    this->id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->id_ion_density = unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS);
    this->id_ion_temperature = unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER);
    this->id_ni = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    this->iIon = iIon;
    AddUnknownForJacobian(unknowns, id_ncold);
    AddUnknownForJacobian(unknowns, id_Tcold);
    AddUnknownForJacobian(unknowns, id_ni);
    AddUnknownForJacobian(unknowns, id_ion_temperature);
    len_t nZ = ions->GetNZ();
    len_t Zmax_max = 0;
    for (len_t iz = 0; iz <nZ; ++iz)
        Zmax_max = std::max(Zmax_max, ions->GetZ(iz));   
    len_t nCharge = Zmax_max + 1;
    size_t deriv_size = static_cast<size_t>(nr) *
                    static_cast<size_t>(nZ) *
                    static_cast<size_t>(nCharge);
    

    Qe_1 = new real_t[nr];
    Qe_2 = new real_t[nr];
    Qe_3 = new real_t[nr];
    d_Qe1_d_Te = new real_t[nr];
    d_Qe1_d_ne = new real_t[nr];
    d_Qe1_d_n_ij = new real_t[deriv_size];
    d_Qe1_d_T_ij = new real_t[deriv_size];
    
    d_Qe2_d_Te = new real_t[nr];
    d_Qe2_d_ne = new real_t[nr];
    d_Qe2_d_n_ij = new real_t[deriv_size];
    d_Qe2_d_T_ij = new real_t[deriv_size];
    
    d_Qe3_d_Te = new real_t[nr];
    d_Qe3_d_ne = new real_t[nr];
    d_Qe3_d_n_ij = new real_t[deriv_size];
    d_Qe3_d_T_ij = new real_t[deriv_size];
}


/**
 * Rebuild the term (called once per time step). Call of the build function in NBIHandler
 */
void NBIIonTerm::Rebuild(const real_t t, const real_t dt, FVM::UnknownQuantityHandler *unknowns){
    real_t beam_energy = handler->GetBeamEnergy();
    

    auto copyData = [&](real_t energy, real_t* Qe_target,
                        real_t* d_Te_target, real_t* d_ne_target,
                        real_t* d_n_ij_target, real_t* d_T_ij_target) {
        handler->Build(t, dt, unknowns, energy);
        
        const real_t* Qe = handler->GetNBIHeatTerm_i();
        const real_t* d_Te = handler->Getd_NBIHeatTerm_i_d_Te();
        const real_t* d_ne = handler->Getd_NBIHeatTerm_i_d_ne();
        const real_t* d_n_ij = handler->Getd_NBIHeatTerm_i_d_n_ij();
        const real_t* d_T_ij = handler->Getd_NBIHeatTerm_i_d_T_ij();

        for (len_t ir = 0; ir < nr; ++ir) {
            Qe_target[ir] = Qe[ir];
            d_Te_target[ir] = d_Te[ir];
            d_ne_target[ir] = d_ne[ir];
        }
        
        for (len_t ir = 0; ir < nr; ++ir) {
            for (len_t iIon = 0; iIon < ions->GetNZ(); ++iIon) {
                for (len_t Zp = 0; Zp <= ions->GetZ(iIon); ++Zp) {
                    len_t idx = handler->idx(ir, iIon, Zp);
                    d_n_ij_target[idx] = d_n_ij[idx];
                    d_T_ij_target[idx] = d_T_ij[idx];
                }
            }
        }
    };
    
    copyData(beam_energy,     Qe_1, d_Qe1_d_Te, d_Qe1_d_ne, d_Qe1_d_n_ij, d_Qe1_d_T_ij);
    copyData(beam_energy / 2, Qe_2, d_Qe2_d_Te, d_Qe2_d_ne, d_Qe2_d_n_ij, d_Qe2_d_T_ij);
    copyData(beam_energy / 3, Qe_3, d_Qe3_d_Te, d_Qe3_d_ne, d_Qe3_d_n_ij, d_Qe3_d_T_ij);

}


/**
 * Set the vector elements corresponding to this term. Also calculates the weight factor for each species. Called for every ion species iIon and charge state Z0.
 */
void NBIIonTerm::SetCSVectorElements(
            real_t *vec, const real_t* /*x*/,
        const len_t iIon, const len_t Z0, const len_t /*rOffset*/
){
    const real_t* energy_fractions = handler->GetEnergyFractions();


    // Calculate the weight factor for each species
    len_t NZ = ions->GetNZ(); 
    if (NZ == 0) {
        throw DREAMException("No ion species found in NBIIonTerm");
    }
    std::vector<real_t> ni_species(NZ, 0.0);   
    std::vector<real_t> Zi(NZ, 0.0);           
    std::vector<real_t> Mi(NZ, 0.0);  
    
    
    for (len_t ir=0; ir<nr; ++ir){
        real_t denom = 0;
        for (len_t iz_denom = 0; iz_denom < ions->GetNZ(); iz_denom++){
                len_t Zmax = ions->GetZ(iz_denom);
                Zi[iz_denom] = Zmax;
                Mi[iz_denom] = ions->GetIonSpeciesMass(iz_denom);

                real_t sum_ni = 0.0;
                for (len_t Z0_inner = 0; Z0_inner <= Zmax; Z0_inner++){
                    real_t niZ = ions->GetIonDensity(ir, iz_denom, Z0_inner);
                    sum_ni += niZ; 
                }
                ni_species[iz_denom] = sum_ni;  
                denom += ni_species[iz_denom] * Zi[iz_denom] / Mi[iz_denom];
        }

        // Only add once per species (when Z0 == 0) since Wi is per species and not per state
        if (Z0 == 0) {
            real_t MiIon = ions->GetIonSpeciesMass(iIon);
            real_t zIon = ions->GetZ(iIon);
            const real_t denom_min = 1e10;
            if (denom < denom_min)
                denom = denom_min;
            real_t w = ni_species[iIon] * zIon / MiIon / denom;
            len_t idx = iIon * nr + ir;
            vec[idx] -= w * (energy_fractions[0] * this->Qe_1[ir] + energy_fractions[1] * this->Qe_2[ir] + energy_fractions[2] * this->Qe_3[ir]);
        }
    }
}
    


/**
 * Set the matrix elements corresponding to this term.
 */
void NBIIonTerm::SetCSMatrixElements(
            FVM::Matrix *mat, real_t *rhs,
        const len_t iIon, const len_t Z0, const len_t /*rOffset*/
){
            (void)mat;

    const real_t* energy_fractions = handler->GetEnergyFractions();


    len_t NZ = ions->GetNZ();
    if (NZ == 0) {
        throw DREAMException("No ion species found in NBIIonTerm");
    }
    std::vector<real_t> ni_species(NZ, 0.0);   
    std::vector<real_t> Zi(NZ, 0.0);           
    std::vector<real_t> Mi(NZ, 0.0);  
    
    
    for (len_t ir=0; ir<nr; ++ir){
        real_t denom = 0;
        for (len_t iz_denom = 0; iz_denom < ions->GetNZ(); iz_denom++){
                len_t Zmax = ions->GetZ(iz_denom);
                Zi[iz_denom] = Zmax;
                Mi[iz_denom] = ions->GetIonSpeciesMass(iz_denom);

                real_t sum_ni = 0.0;
                for (len_t Z0_inner = 0; Z0_inner <= Zmax; Z0_inner++){
                    real_t niZ = ions->GetIonDensity(ir, iz_denom, Z0_inner);
                    sum_ni += niZ; 
                }
                ni_species[iz_denom] = sum_ni;  
                denom += ni_species[iz_denom] * Zi[iz_denom] / Mi[iz_denom];
        }

        if (Z0 == 0) {
            real_t MiIon = ions->GetIonSpeciesMass(iIon);
            real_t zIon = ions->GetZ(iIon);
            const real_t denom_min = 1e10;
            if (denom < denom_min)
                denom = denom_min;
            real_t w = ni_species[iIon] * zIon / MiIon / denom;
            len_t idx = iIon * nr + ir;
            rhs[idx] -= w * (energy_fractions[0] * this->Qe_1[ir] + energy_fractions[1] * this->Qe_2[ir] + energy_fractions[2] * this->Qe_3[ir]);
        }
    }
    
}
/**
 * Set the Jacobian elements corresponding to this term.
 */

bool NBIIonTerm::SetCSJacobianBlock(
           const len_t /*uqtyId*/, const len_t derivId,
                        FVM::Matrix *jac, const real_t* /*x*/,
                        const len_t iIon, const len_t Z0, const len_t /*rOffset*/
){
    const real_t* energy_fractions = handler->GetEnergyFractions();

    if (derivId != id_ncold && derivId != id_Tcold &&
        derivId != id_ni && derivId != id_ion_temperature) {
        return false;
    }
    // For ion density derivatives, only set for Z0 == 0 (ion species), since Wi is per species
    if (derivId != id_ni && Z0 != 0) {
        return true;
    }

    // Set Jacobian elements for each radial point
    if (derivId == id_ncold) {
        for (len_t ir = 0; ir < nr; ++ir) {
            len_t row = iIon*nr+ir;
            len_t col = ir;
            real_t deriv = energy_fractions[0] * d_Qe1_d_ne[ir]
                         + energy_fractions[1] * d_Qe2_d_ne[ir]
                         + energy_fractions[2] * d_Qe3_d_ne[ir];
           jac->SetElement(row, col, deriv); 
        }
    }
    
    if (derivId == id_Tcold) {
        for (len_t ir = 0; ir < nr; ++ir) {
            len_t row = iIon*nr+ir;
            len_t col = ir;
            real_t deriv = energy_fractions[0] * d_Qe1_d_Te[ir]
                         + energy_fractions[1] * d_Qe2_d_Te[ir]
                         + energy_fractions[2] * d_Qe3_d_Te[ir];
           jac->SetElement(row, col, deriv); 
        }
    }
    
    if (derivId == id_ni) {
        for (len_t ir = 0; ir < nr; ++ir) {
            len_t row = iIon * nr + ir;
            len_t col = ions->GetIndex(iIon, Z0) * nr + ir;
            len_t idx = handler->idx(ir, iIon, Z0);
            real_t deriv = energy_fractions[0] * d_Qe1_d_n_ij[idx]
                         + energy_fractions[1] * d_Qe2_d_n_ij[idx]
                         + energy_fractions[2] * d_Qe3_d_n_ij[idx];
            jac->SetElement(row, col, deriv);
        }
    }
    
    if (derivId == id_ion_temperature) {
        for (len_t ir = 0; ir < nr; ++ir) {
            
            
            const len_t Zmax = ions->GetZ(iIon);
            real_t ni_species = 0.0;
            real_t dQi_dTi_total = 0.0;
            //sum over charges
            for (len_t Zp = 0; Zp <= Zmax; ++Zp) {
                ni_species += ions->GetIonDensity(ir, iIon, Zp);
                len_t idx = handler->idx(ir, iIon, Zp);
                dQi_dTi_total += energy_fractions[0] * d_Qe1_d_T_ij[idx]
                           + energy_fractions[1] * d_Qe2_d_T_ij[idx]
                           + energy_fractions[2] * d_Qe3_d_T_ij[idx];
            }

            
            const real_t eq = 1.602e-19; 
            real_t dP = 0.0;
            if (ni_species > 1e10) { 
                dP = dQi_dTi_total / (1.5 * ni_species * eq);
            }
            len_t col = iIon*nr+ir;
            len_t row = iIon*nr+ir;
            jac->SetElement(row, col, dP);
        }
    }
    return true;
}



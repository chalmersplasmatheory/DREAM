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
    len_t nZ = ions->GetNZ();
    len_t Zmax_max = 0;
    for (len_t iz = 0; iz <nZ; ++iz)
        Zmax_max = std::max(Zmax_max, ions->GetZ(iz));   
    len_t nCharge = Zmax_max + 1;
    size_t deriv_size = static_cast<size_t>(nr) *
                    static_cast<size_t>(nZ) *
                    static_cast<size_t>(nCharge);
    

    // Allocate temporary arrays for energy spectra and their derivatives
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

    // Register all unknowns that this term has Jacobian contributions for
    AddUnknownForJacobian(unknowns, id_ncold);
    AddUnknownForJacobian(unknowns, id_Tcold);
    AddUnknownForJacobian(unknowns, id_ni);
    AddUnknownForJacobian(unknowns, id_ion_temperature);
}


/**
 * Rebuild the term (called once per time step). Call of the build function in NBIHandler
 */
void NBIElectronTerm::Rebuild(const real_t t, const real_t dt, FVM::UnknownQuantityHandler *unknowns){
    real_t beam_energy = handler->GetBeamEnergy();
    
    // Helper to copy data
    auto copyData = [&](real_t energy, real_t* Qe_target,
                        real_t* d_Te_target, real_t* d_ne_target,
                        real_t* d_n_ij_target, real_t* d_T_ij_target) {
        handler->Build(t, dt, unknowns, energy);
        
        const real_t* Qe = handler->GetNBIHeatTerm_e();
        const real_t* d_Te = handler->Getd_NBIHeatTerm_e_d_Te();
        const real_t* d_ne = handler->Getd_NBIHeatTerm_e_d_ne();
        const real_t* d_n_ij = handler->Getd_NBIHeatTerm_e_d_n_ij();
        const real_t* d_T_ij = handler->Getd_NBIHeatTerm_e_d_T_ij();
        
        // Copy heating terms
        for (len_t ir = 0; ir < nr; ++ir) {
            Qe_target[ir] = Qe[ir];
            d_Te_target[ir] = d_Te[ir];
            d_ne_target[ir] = d_ne[ir];
        }
        
        // Copy ion derivatives (multi-dimensional)
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
    
    // Build for three energy components
    copyData(beam_energy,     Qe_1, d_Qe1_d_Te, d_Qe1_d_ne, d_Qe1_d_n_ij, d_Qe1_d_T_ij);
    copyData(beam_energy / 2, Qe_2, d_Qe2_d_Te, d_Qe2_d_ne, d_Qe2_d_n_ij, d_Qe2_d_T_ij);
    copyData(beam_energy / 3, Qe_3, d_Qe3_d_Te, d_Qe3_d_ne, d_Qe3_d_n_ij, d_Qe3_d_T_ij);

}

/**
 * Set the vector elements corresponding to this term.
 */
void NBIElectronTerm::SetVectorElements(real_t *rhs, const real_t*){
    const real_t* energy_fractions = handler->GetEnergyFractions();
    for (len_t ir=0; ir<nr; ++ir) {   
        rhs[ir] -= energy_fractions[0] * this->Qe_1[ir] + energy_fractions[1] * this->Qe_2[ir] + energy_fractions[2] * this->Qe_3[ir];
    }
}

/**
 * Set the matrix elements corresponding to this term.
 */
void NBIElectronTerm::SetMatrixElements(FVM::Matrix*, real_t *rhs){ 
    const real_t* energy_fractions = handler->GetEnergyFractions();
    for (len_t ir=0; ir<nr; ++ir) {
        rhs[ir] -= energy_fractions[0] * this->Qe_1[ir] + energy_fractions[1] * this->Qe_2[ir] + energy_fractions[2] * this->Qe_3[ir];
    }
}

/**
 * Set the Jacobian elements corresponding to this term.
 */
bool NBIElectronTerm::SetJacobianBlock(const len_t, const len_t derivId,
                                       FVM::Matrix *jac, const real_t*){
    if (derivId != id_ncold && derivId != id_Tcold &&
        derivId != id_ni && derivId != id_ion_temperature)
        return false;
    
    const real_t* energy_fractions = handler->GetEnergyFractions();    
    if (derivId == id_Tcold) {
        for (len_t ir = 0; ir < nr; ++ir) {
            real_t deriv = energy_fractions[0] * d_Qe1_d_Te[ir]
                         + energy_fractions[1] * d_Qe2_d_Te[ir]
                         + energy_fractions[2] * d_Qe3_d_Te[ir];
            jac->SetElement(ir, ir, deriv);
        }
    }
    
    if (derivId == id_ncold) {
        for (len_t ir = 0; ir < nr; ++ir) {
            real_t deriv = energy_fractions[0] * d_Qe1_d_ne[ir]
                         + energy_fractions[1] * d_Qe2_d_ne[ir]
                         + energy_fractions[2] * d_Qe3_d_ne[ir];
            jac->SetElement(ir, ir, deriv);
        }
    }
    
    if (derivId == id_ni) {
        for (len_t ir = 0; ir < nr; ++ir) {
            for (len_t iIon = 0; iIon < ions->GetNZ(); ++iIon) {
                const len_t Zmax = ions->GetZ(iIon);
                const len_t speciesStart = ions->GetIndex(iIon, 0) * nr;
                for (len_t Zp = 0; Zp <= Zmax; ++Zp) {
                    len_t idx = handler->idx(ir, iIon, Zp);
                    real_t deriv = energy_fractions[0] * d_Qe1_d_n_ij[idx]
                                 + energy_fractions[1] * d_Qe2_d_n_ij[idx]
                                 + energy_fractions[2] * d_Qe3_d_n_ij[idx];
                    len_t col = speciesStart + Zp * nr + ir;
                    jac->SetElement(ir, col, deriv);
                }
            }
        }
    }
    
    if (derivId == id_ion_temperature) {
        for (len_t ir = 0; ir < nr; ++ir) {
            for (len_t iIon = 0; iIon < ions->GetNZ(); ++iIon) {
                real_t deriv_sum = 0.0;
                for (len_t Zp = 0; Zp <= ions->GetZ(iIon); ++Zp) {
                    len_t idx = handler->idx(ir, iIon, Zp);
                    deriv_sum += energy_fractions[0] * d_Qe1_d_T_ij[idx]
                               + energy_fractions[1] * d_Qe2_d_T_ij[idx]
                               + energy_fractions[2] * d_Qe3_d_T_ij[idx];
                }
                len_t col = iIon * nr + ir;
                jac->SetElement(ir, col, deriv_sum);
            }
        }
    }
    
    return true;
}



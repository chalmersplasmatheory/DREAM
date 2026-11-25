#include "DREAM/Equations/Fluid/NBIIonTerm.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/DREAMException.hpp"
#include "DREAM/Constants.hpp"
#include "DREAM/IO.hpp"



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
    len_t deriv_size =(nr * nZ * nCharge);
    
    Qe_1 = new real_t[nr];
    Qe_2 = new real_t[nr];
    Qe_3 = new real_t[nr];
    // Derivatives are non-local: [ir*nr + k] and [ir*deriv_size + idx(k, iz, Z0)]
    this->deriv_size = deriv_size;
    d_Qe1_d_Te = new real_t[nr * nr];
    d_Qe1_d_ne = new real_t[nr * nr];
    d_Qe1_d_n_ij = new real_t[nr * deriv_size];
    d_Qe1_d_T_ij = new real_t[nr * deriv_size];
    
    d_Qe2_d_Te = new real_t[nr * nr];
    d_Qe2_d_ne = new real_t[nr * nr];
    d_Qe2_d_n_ij = new real_t[nr * deriv_size];
    d_Qe2_d_T_ij = new real_t[nr * deriv_size];
    
    d_Qe3_d_Te = new real_t[nr * nr];
    d_Qe3_d_ne = new real_t[nr * nr];
    d_Qe3_d_n_ij = new real_t[nr * deriv_size];
    d_Qe3_d_T_ij = new real_t[nr * deriv_size];
}
NBIIonTerm::~NBIIonTerm(){
        delete[] Qe_1;
        delete[] Qe_2;
        delete[] Qe_3;
        delete[] d_Qe1_d_Te;
        delete[] d_Qe1_d_ne;
        delete[] d_Qe1_d_n_ij;
        delete[] d_Qe1_d_T_ij;
        delete[] d_Qe2_d_Te;
        delete[] d_Qe2_d_ne;
        delete[] d_Qe2_d_n_ij;
        delete[] d_Qe2_d_T_ij;
        delete[] d_Qe3_d_Te;
        delete[] d_Qe3_d_ne;
        delete[] d_Qe3_d_n_ij;
        delete[] d_Qe3_d_T_ij;
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

        for (len_t ir = 0; ir < nr; ++ir)
            Qe_target[ir] = Qe[ir];
        for (len_t i = 0; i < nr * nr; ++i) {
            d_Te_target[i] = d_Te[i];
            d_ne_target[i] = d_ne[i];
        }
        for (len_t i = 0; i < nr * deriv_size; ++i) {
            d_n_ij_target[i] = d_n_ij[i];
            d_T_ij_target[i] = d_T_ij[i];
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
    if (Z0 != 0) return; 
    const real_t* energy_fractions = handler->GetEnergyFractions();
    
    len_t NZ = ions->GetNZ(); 
    if (NZ == 0) {
        throw DREAMException("No ion species found in NBIIonTerm");
    }
    std::vector<real_t> ni_species(NZ, 0.0);   
    std::vector<real_t> Zi(NZ, 0.0);           
    std::vector<real_t> Mi(NZ, 0.0);  
    
    for (len_t ir=0; ir<nr; ++ir){
        real_t w = ComputeWeightFactor(ir, iIon);
        len_t idx = iIon * nr + ir;
        vec[idx] += w * (energy_fractions[0] * this->Qe_1[ir] + energy_fractions[1] * this->Qe_2[ir] + energy_fractions[2] * this->Qe_3[ir]); 
    }
}
    

/**
 * Set the matrix elements corresponding to this term.
 */
void NBIIonTerm::SetCSMatrixElements(
            FVM::Matrix* /*mat*/, real_t *rhs,
        const len_t iIon, const len_t Z0, const len_t /*rOffset*/
){
    if (Z0 != 0) return; 
    const real_t* energy_fractions = handler->GetEnergyFractions();
    len_t NZ = ions->GetNZ();
    if (NZ == 0) {
        throw DREAMException("No ion species found in NBIIonTerm");
    }
    std::vector<real_t> ni_species(NZ, 0.0);   
    std::vector<real_t> Zi(NZ, 0.0);           
    std::vector<real_t> Mi(NZ, 0.0);  
    
    for (len_t ir=0; ir<nr; ++ir){
        real_t w = ComputeWeightFactor(ir, iIon);
        len_t idx = iIon * nr + ir;
        rhs[idx] += w * (energy_fractions[0] * this->Qe_1[ir] + energy_fractions[1] * this->Qe_2[ir] + energy_fractions[2] * this->Qe_3[ir]);
    }
    
}
/**
 * Set the Jacobian elements corresponding to this term.
 * The ion heating of species iIon in cell ir is w_iIon(ir) * Q_i(ir). Q_i depends on the
 * plasma in every cell k that the beam passes before ir (non-local), while the weight
 * w depends on the local charge-state densities of all species.
 */
bool NBIIonTerm::SetCSJacobianBlock(
           const len_t /*uqtyId*/, const len_t derivId,
                        FVM::Matrix *jac, const real_t* /*x*/,
                        const len_t iIon, const len_t Z0, const len_t /*rOffset*/
){
    // The ion energy equation has one row per species, so only act once (Z0 = 0)
    if (Z0 != 0)
        return false;

    const real_t* energy_fractions = handler->GetEnergyFractions(); 

    if (derivId == id_ncold || derivId == id_Tcold) {
        const real_t *d1 = (derivId == id_ncold ? d_Qe1_d_ne : d_Qe1_d_Te);
        const real_t *d2 = (derivId == id_ncold ? d_Qe2_d_ne : d_Qe2_d_Te);
        const real_t *d3 = (derivId == id_ncold ? d_Qe3_d_ne : d_Qe3_d_Te);
        for (len_t ir = 0; ir < nr; ++ir) {
            len_t row = iIon*nr+ir;
            real_t w = ComputeWeightFactor(ir, iIon);
            for (len_t k = 0; k < nr; ++k) {
                real_t deriv = energy_fractions[0] * d1[ir * nr + k]
                             + energy_fractions[1] * d2[ir * nr + k]
                             + energy_fractions[2] * d3[ir * nr + k];
                if (deriv != 0)
                    jac->SetElement(row, k, w*deriv);
            }
        }
    }
    
    else if (derivId == id_ni) {
        for (len_t ir = 0; ir < nr; ++ir) {
            len_t row = iIon * nr + ir;
            real_t w = ComputeWeightFactor(ir, iIon);
            real_t Qi = energy_fractions[0] * Qe_1[ir] + energy_fractions[1] * Qe_2[ir] + energy_fractions[2] * Qe_3[ir];
            for (len_t k = 0; k < nr; ++k) {
                for (len_t jIon = 0; jIon < ions->GetNZ(); ++jIon) {
                    for (len_t Zp = 0; Zp <= ions->GetZ(jIon); ++Zp) {
                        len_t col = ions->GetIndex(jIon, Zp) * nr + k;
                        len_t idx = ir * deriv_size + handler->idx(k, jIon, Zp);
                        real_t deriv = energy_fractions[0] * d_Qe1_d_n_ij[idx]
                                     + energy_fractions[1] * d_Qe2_d_n_ij[idx]
                                     + energy_fractions[2] * d_Qe3_d_n_ij[idx];
                        real_t val = w * deriv;
                        // The weight factor only depends on the local densities
                        if (k == ir)
                            val += ComputeWeightFactorDerivative(ir, iIon, jIon, Zp) * Qi;
                        if (val != 0)
                            jac->SetElement(row, col, val);
                    }
                }
            }
        }
    }
    
    else if (derivId == id_ion_temperature) {
        // The unknown is W_j = 1.5 e n_j T_j, and T_j enters through the CX rates
        const real_t eq = Constants::ec;
        for (len_t ir = 0; ir < nr; ++ir) {
            len_t row = iIon*nr+ir;
            real_t w = ComputeWeightFactor(ir, iIon);
            for (len_t k = 0; k < nr; ++k) {
                for (len_t jIon = 0; jIon < ions->GetNZ(); ++jIon) {
                    real_t nj = 0.0;
                    real_t dQi_dTj = 0.0;
                    //Sum over charges
                    for (len_t Zp = 0; Zp <= ions->GetZ(jIon); ++Zp) {
                        nj += ions->GetIonDensity(k, jIon, Zp);
                        len_t idx = ir * deriv_size + handler->idx(k, jIon, Zp);
                        dQi_dTj += energy_fractions[0] * d_Qe1_d_T_ij[idx]
                                 + energy_fractions[1] * d_Qe2_d_T_ij[idx]
                                 + energy_fractions[2] * d_Qe3_d_T_ij[idx];
                    }
                    len_t col = jIon*nr+k;
                    if (dQi_dTj != 0)
                        jac->SetElement(row, col, w * dQi_dTj / (1.5 * nj * eq));
                }
            }
        }
    }
    else {
        return false;
    }
    return true;
}

/**
 * Compute the weight factor for ion species iIon at radial index ir.
 * The ion share of the beam power is split between species in proportion to
 * sum_{Z0>=1} n_{j,Z0} Z0^2 / M_j (Coulomb collisions; neutrals receive nothing).
 */
real_t NBIIonTerm::ComputeWeightFactor(len_t ir, len_t iIon) {
    len_t NZ = ions->GetNZ();
    if (NZ == 0) {
        throw DREAMException("No ion species found in NBIIonTerm");
    }

    real_t denom = 0;
    real_t numer = 0;
    for (len_t iz_denom = 0; iz_denom < NZ; iz_denom++){
            real_t Mi = ions->GetIonSpeciesMass(iz_denom);
            real_t sum_nZ2 = 0.0;
            for (len_t Z0_inner = 1; Z0_inner <= ions->GetZ(iz_denom); Z0_inner++){
                real_t niZ = ions->GetIonDensity(ir, iz_denom, Z0_inner);
                sum_nZ2 += niZ * Z0_inner * Z0_inner;
            }
            denom += sum_nZ2 / Mi;
            if (iz_denom == iIon)
                numer = sum_nZ2 / Mi;
    }

    const real_t denom_min = 1e5;
    if (denom < denom_min){
       DREAM::IO::PrintWarning(
				DREAM::IO::WARNING_ION_DENSITY_TOO_SMALL,
				"Warning: NBIIonTerm weight factor denom too small."
			);
        denom = denom_min;
    }
    real_t w = numer / denom;
    return w;
}

/**
 * Derivative of the weight factor of species iIon with respect to the charge state
 * density n_{jIon,Z0} (same cell):
 *   dw_i/dn_{j,Z0} = (delta_ij - w_i) * (Z0^2/M_j) / denom
 */
real_t NBIIonTerm::ComputeWeightFactorDerivative(len_t ir, len_t iIon, len_t jIon, len_t Z0) {
    len_t NZ = ions->GetNZ();
    if (NZ == 0) {
        throw DREAMException("No ion species found in NBIIonTerm");
    }

    real_t denom = 0;
    for (len_t iz_denom = 0; iz_denom < NZ; iz_denom++){
            real_t Mi = ions->GetIonSpeciesMass(iz_denom);
            for (len_t Z0_inner = 1; Z0_inner <= ions->GetZ(iz_denom); Z0_inner++){
                real_t niZ = ions->GetIonDensity(ir, iz_denom, Z0_inner);
                denom += niZ * Z0_inner * Z0_inner / Mi;
            }
    }

    const real_t denom_min = 1e5;
    if (denom < denom_min){
        DREAM::IO::PrintWarning(
				DREAM::IO::WARNING_ION_DENSITY_TOO_SMALL,
				"Warning: NBIIonTerm weight factor denom too small."
			);
        denom = denom_min;
    }
    real_t Mj = ions->GetIonSpeciesMass(jIon);
    real_t w = ComputeWeightFactor(ir, iIon);
    real_t dwdni = ((iIon == jIon ? 1.0 : 0.0) - w) * (real_t)(Z0 * Z0) / (Mj * denom);
    return dwdni;
}

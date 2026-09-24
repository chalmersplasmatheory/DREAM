/**
 * Implementation of a neutral collisional energy transfer term. 
 * Calculates the energy transfer between two neutral species due to collisions.
 */

#include "DREAM/Equations/Fluid/NeutralCollisionalEnergyTransferTerm.hpp"
#include "DREAM/Constants.hpp"
#include <cmath>

using namespace DREAM;

/**
 * Constructor
 */
NeutralCollisionalEnergyTransferTerm::NeutralCollisionalEnergyTransferTerm(
    FVM::Grid *grid, len_t iz, len_t jz, FVM::UnknownQuantityHandler *unknowns,IonHandler *ions
): FVM::EquationTerm(grid), unknowns(unknowns), ions(ions), iz(iz), jz(jz){
    
    SetName("NeutralCollisionalEnergyTransferTerm");

    this->mi = ions->GetIonSpeciesMass(iz);
    this->mj = ions->GetIonSpeciesMass(jz);
    this->id_ions =unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    this->id_Wn =unknowns->GetUnknownID(OptionConstants::UQTY_WN_ENER);

    AddUnknownForJacobian(unknowns,id_ions);
    AddUnknownForJacobian(unknowns,id_Wn);
}

/**
 * Gets the parameters of a neutral species at a given radial position.
 */
 void NeutralCollisionalEnergyTransferTerm::GetNeutralParameters(
      len_t ir,
      real_t& n,
      real_t& W,
      real_t& T
  ) {
    const real_t *ni =unknowns->GetUnknownData(id_ions);
    const real_t *Wn =unknowns->GetUnknownData(id_Wn);
    const len_t neutralIndex =ions->GetIndex(iz, 0);

    n = ni[neutralIndex*nr + ir];
    W = Wn[iz*nr + ir];
    if (n > 0)
        T = W / (1.5 * Constants::ec * n);
    else
        T = 0;
  }

/**
 * Calculates the momentum scattering radius for a given neutral species. In units of Ångström.
 */
real_t NeutralCollisionalEnergyTransferTerm::getMomentumScatteringRadius(len_t IZ){
    //TODO: Look this up
    const std::string name = ions->GetName(IZ);
    if (name == "D")
        return 1.35;

    else if (name == "D2")
          return 2.01;

    else if (name == "D3")
        return 2.01;   // must be confirmed

    else if (name == "He")
        return 1.47;

    else if (name == "Ne")
        return 1.91;

    else if (name == "Ar")
          return 2.45;
    else{
        return 0;
    }
}

/**
 * Calculates the neutral scattering rate between two neutral species at a given radial position.
 * The scattering rate is calculated using the momentum scattering radius of the two species, their masses,
 * and their temperatures.
 */
real_t NeutralCollisionalEnergyTransferTerm::CalculateNeutralScatteringRate(
      len_t iz, len_t jz, len_t ir, real_t mi, real_t mj, real_t Ti, real_t Tj, real_t nj
) {

    real_t ri = getMomentumScatteringRadius(iz) * 1e-10; // convert to meters
    real_t rj = getMomentumScatteringRadius(jz) * 1e-10; // convert to meters
    const real_t miAMU = ions->GetIonSpeciesMass(iz)/ Constants::mu;
    const real_t mjAMU =ions->GetIonSpeciesMass(jz)/ Constants::mu;

    const real_t reducedTemperature =(mjAMU*Ti + miAMU*Tj) / (miAMU + mjAMU);
    const real_t reducedMass =(miAMU*mjAMU)/ (miAMU + mjAMU);

    const real_t v_ij =std::sqrt(2 * Constants::ec * reducedTemperature/ (reducedMass * Constants::mu));
    //const real_t r_ij = 0.5 *(ri+rj);
    const real_t r_ij = (ri+rj);
    //const real_t vNN_inj = 2*std::sqrt(M_PI)/3 * nj *v_ij*(mj/(mi + mj))* r_ij*r_ij;
    const real_t vNN_inj = 2*2*std::sqrt(M_PI)/3 * nj *v_ij* r_ij*r_ij;

    return vNN_inj;
  }

/**
 * Sets the elements of the vector based on the neutral collisional energy transfer.
 */
void NeutralCollisionalEnergyTransferTerm::SetVectorElements(
    real_t *vec,
    const real_t*
) {
    const real_t *densities = unknowns->GetUnknownData(id_ions);
    const real_t *energies = unknowns->GetUnknownData(id_Wn);

    const len_t neutralIndexI = ions->GetIndex(iz, 0);
    const len_t neutralIndexJ = ions->GetIndex(jz, 0);

    const real_t mi = ions->GetIonSpeciesMass(iz);
    const real_t mj = ions->GetIonSpeciesMass(jz);
    
    const real_t massFactor = 2 * mi * mj / ((mi + mj) * (mi + mj));

    for (len_t ir = 0; ir < nr; ir++) {
        const real_t ni = densities[neutralIndexI * nr + ir];
        const real_t nj = densities[neutralIndexJ * nr + ir];

        if (ni <= 0 || nj <= 0)
            continue;

        const real_t Wi = energies[iz * nr + ir];
        const real_t Wj = energies[jz * nr + ir];

        const real_t Ti = Wi / (1.5 * Constants::ec * ni);
        const real_t Tj = Wj / (1.5 * Constants::ec * nj);

        const real_t nuNN_ij = CalculateNeutralScatteringRate(iz, jz, ir, mi, mj, Ti, Tj, nj);
        const real_t nuT_ij = massFactor * nuNN_ij;

        // Power density transferred from species j to species i [W/m^3].
        const real_t Q =1.5 * Constants::ec * ni * nuT_ij * (Tj - Ti);

        // Update the energy densities of species i and j.
        vec[iz * nr + ir] += Q;
        vec[jz * nr + ir] -= Q;
    } 
  }

/**
 * Sets the elements of the Jacobian matrix based on the neutral 
 * collisional energy transfer. TODO: Check this math
 */
bool NeutralCollisionalEnergyTransferTerm::SetJacobianBlock(
    const len_t,const len_t derivId,
    FVM::Matrix *jac,
    const real_t*
  ) {
      if (derivId != id_Wn && derivId != id_ions)
          return false;

    const real_t *densities = unknowns->GetUnknownData(id_ions);
    const real_t *energies = unknowns->GetUnknownData(id_Wn);

    const len_t neutralIndexI = ions->GetIndex(iz, 0);
    const len_t neutralIndexJ = ions->GetIndex(jz, 0);

    const real_t A = 1.5 * Constants::ec;
    const real_t massFactor = 2 * mi * mj / ((mi + mj) * (mi + mj));

    bool contributes = false;

    for (len_t ir = 0; ir < nr; ir++) {
        const len_t rowI = iz * nr + ir;
        const len_t rowJ = jz * nr + ir;
        const len_t densityColI = neutralIndexI * nr + ir;
        const len_t densityColJ = neutralIndexJ * nr + ir;

        const real_t ni = densities[densityColI];
        const real_t nj = densities[densityColJ];

        if (ni <= 0 || nj <= 0)
            continue;

        const real_t Ti = energies[rowI] / (A * ni);
        const real_t Tj = energies[rowJ] / (A * nj);

        const real_t nuT = massFactor * CalculateNeutralScatteringRate(
            iz, jz, ir, mi, mj, Ti, Tj, nj);

        // At Ti = Tj = 0, the first derivatives vanish.
        if (Ti == 0 && Tj == 0)
            continue;

        const real_t s = Ti / mi + Tj / mj;
        const real_t dT = Tj - Ti;
        const real_t prefactor = A * ni * nuT;
        const real_t Q = prefactor * dT;

        // Temperature derivatives at fixed densities.
        const real_t dQdTi =
            prefactor * (dT / (2 * mi * s) - 1);
        const real_t dQdTj =
            prefactor * (dT / (2 * mj * s) + 1);

        real_t dQdxI, dQdxJ;
        len_t colI, colJ;

        if (derivId == id_Wn) {
            colI = rowI;
            colJ = rowJ;

            dQdxI = dQdTi / (A * ni);
            dQdxJ = dQdTj / (A * nj);
        } else {
            colI = densityColI;
            colJ = densityColJ;

            // Density derivatives at fixed energy densities.
            dQdxI = (Q - Ti * dQdTi) / ni;
            dQdxJ = (Q - Tj * dQdTj) / nj;
        }

        // Species i receives Q; species j loses exactly Q.
        jac->SetElement(rowI, colI,  dQdxI);
        jac->SetElement(rowI, colJ,  dQdxJ);
        jac->SetElement(rowJ, colI, -dQdxI);
        jac->SetElement(rowJ, colJ, -dQdxJ);

        contributes = true;
    }

    return contributes;
  }


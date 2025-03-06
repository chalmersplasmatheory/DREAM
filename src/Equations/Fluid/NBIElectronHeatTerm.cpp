/**
 * Implementation of equation term to the heating of cold electrons from an NBI beam 
 */
#include <cmath>
#include <array>
#include "DREAM/Equations/Fluid/NBIElectronHeatTerm.hpp"
#include "DREAM/Constants.hpp"
#include "DREAM/ADAS.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM{
/**
 * Constructor
 */
NBIElectronHeatTerm::NBIElectronHeatTerm(FVM::Grid* grid, FVM::UnknownQuantityHandler *unknowns):EquationTerm(grid){
//Id for unknown quantaties, allows to query actual numerical data during a simulation step.
id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD); //density
id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD); // Plasma Temp
//add more dependependenciesdencaces


//Add dependencies for jacobian
AddUnknownForJacobian(unknowns,id_ncold);
AddUnknownForJacobian(unknowns, id_Tcold);
//add more

//get the grid number
this->nr = grid->GetNr();
NBIHeatTerm = new real_t[nr];

}

/**
 * Destructor, deallocate memory in some way
 */
NBIElectronHeatTerm::~NBIElectronHeatTerm() {
    delete[] NBIHeatTerm;
}

/**
 * Rebuild: Called at the start of every timestep.
 * Compute added heating term for each timestep and flux surface radius
 */
void NBIElectronHeatTerm::Rebuild(const real_t t, const real_t dt, UnknownQuantityHandler* unknowns){
    this-> dt = dt;

    for (len_t ir = 0; ir < nr;ir ++){
        real_t ncold = unknowns ->GetUnknownData(id_ncold)[ir];
        real_t Tcold = unknowns ->GetUnknownData(id_Tcold)[ir];
        //add more dependancies
        //put more dependancies into compute depositions profile
        NBIHeatTerm[ir] = beamCurrentDensity/beamIntensity *ComputeDepositionProfile(lambda_s, survivalProb); //compute the NBI heating term * times the current density and beam intensisty
    }
}

/**
 * Converts a flux surface index ir into cylindrical coordinates along the beamline. To map plasma profiles (like density) into the beam coordinate system
 */
void NBIElectronHeatTerm::ComputeCylindricalCoordinates(len_t ir, real_t &r_B, real_t &theta_B, real_t &s_B){
//Part of it is already inbuild

//takes in ir
//turns it to x,y,z 
//turns it to cylindrical

}

/**
 *  Computes the total  cross-section at a flux surface.
 */
real_t NBIElectronHeatTerm::GetTotalCrossSection(len_t ir, real_t ncold, real_t Tcold){
//Gets and sums the total cross section from other files, ADAS, need more dependancies
}

/**
 * Computes the Mean Free Path
 */
real_t NBIElectronHeatTerm::ComputeMeanFreePath(real_t ncold, real_t sigma){
//Computes the mean free path depending on the total cross section and plasma density, called at every flux surface
    return 1.0 / (ncold * sigma);
}

/**
 * Computes the exponent corresponding to the survival probability
 */
real_t NBIElectronHeatTerm::ComputeSurvivalProbability(real_t s, real_t lambda_s){
//Should be some integration along the beam path up to that point. Depends on s and the mean free path
}


/**
 * Computes the  deposition profile
 */
real_t NBIElectronHeatTerm::ComputeDepositionProfile(real_t lambda_s, real_t survivalProb){
// depends on lambda and the survival probability, but should also be inegrated over the angles. Gives H(r). Is Jb and Ib constant look up
        r_B,theta_B,s_B = ComputeCylindricalCoordinates(ir, r_B, theta_B, s_B);

        real_t sigma = GetTotalCrossSection(ir,ncold,Tcold); //compute the total cross section
        real_t lambda_s = ComputeMeanFreePath(ncold, sigma); // compute the mean free path
        real_t survivalProb = ComputeSurvivalProbability(s_B, lambda_s); //compute the exponent survavalprobability

        //add the actual computing/integrating for that flux surface
}




void NBIElectronHeatTerm::SetJacobianBlock(const len_t uqtyId, const len_t derivId, Matrix*, const real_t*) {}

void NBIElectronHeatTerm::SetMatrixElements(FVM::Matrix *mat, real_t *rhs) {}

void NBIElectronHeatTerm::SetVectorElements(real_t *rhs) {}




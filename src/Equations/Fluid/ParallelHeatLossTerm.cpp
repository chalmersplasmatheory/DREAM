/**
 * Heat flux parallel loss term in the plasma HALO region.
 */

#include "DREAM/DREAMException.hpp"
#include "DREAM/Equations/Fluid/ParallelHeatLossTerm.hpp"
#include "DREAM/Equations/Scalar/WallCurrentTerms.hpp"

using namespace DREAM;

/**
 * Constructor.
 */
HaloRegionHeatLossTerm::HaloRegionHeatLossTerm(
	FVM::Grid *grid, FVM::UnknownQuantityHandler *unknowns, IonHandler *ions,
    real_t sf, bool userGivenPsiEdge_t0, real_t PsiEdge_t0
) : FVM::DiagonalComplexTerm(grid, unknowns, grid),
	unknowns(unknowns), ions(ions), scaleFactor(sf), 
    userGivenPsiEdge_t0(userGivenPsiEdge_t0), psi_edge_t0(PsiEdge_t0), 
    id_psi(unknowns->GetUnknownID(OptionConstants::UQTY_POL_FLUX)){

    SetName("HaloRegionHeatLossTerm");
    this->rGrid = grid->GetRadialGrid();
    this->GridRebuilt();

    this->id_T_cold  = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->id_W_cold = unknowns->GetUnknownID(OptionConstants::UQTY_W_COLD);
    this->id_N_i = unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS);
    this->id_W_i = unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER);
    this->id_jtot  = unknowns->GetUnknownID(OptionConstants::UQTY_J_TOT);

    // Loop through all ions to get the one with the minimum Z
    for (len_t i = 0; i < ions->GetNZ(); i++) {
        int currentZ = ions->GetZ(i); 
        if (currentZ < minZ) {        
            minZ = currentZ;          
            minIndex = i;             
        }
    }

    // Check if a valid minimum Z was found and set the ion quantities
    if (minIndex != -1) {
        this-> m_i = ions->GetIonSpeciesMass(minIndex); 
        this-> Z = minZ;  
                                     
    }
}

/**
 * Destructor
 */
HaloRegionHeatLossTerm::~HaloRegionHeatLossTerm(){
    Deallocate();

}           


/**
 * Called after the grid is rebuilt; (re)allocates memory
 * for all quantities
 */
bool HaloRegionHeatLossTerm::GridRebuilt(){

    this->DiagonalComplexTerm::GridRebuilt();
    Deallocate();
    this->HeatLoss = new real_t[nr];
    
    return true;
}


/**
* Find current radial point of last closed flux surface
* and find psi_edge_t0 if not provided in input
*/
void HaloRegionHeatLossTerm::FindRadiusOfLCFS(){
    
    if(this->userGivenPsiEdge_t0 == 0) {
        
        // Extrapolate psi from last two ir:s to find the approximate edge value at t=0
        real_t Rlow = this->rGrid->GetR(nr-2);
        real_t Rhigh = this->rGrid->GetR(nr-1);
        real_t Psilow = unknowns->GetUnknownInitialData(id_psi)[nr-2];
        real_t Psihigh = unknowns->GetUnknownInitialData(id_psi)[nr-1];
        real_t slope = (Psihigh - Psilow) / (Rhigh - Rlow);
        real_t Redge = this->rGrid->GetR_f(nr);
        
        this->psi_edge_t0 = Psihigh + slope * (Redge - Rhigh); 
    }
    
    // Find if the sign of PsiDiff() should be +1 or -1
    if(this->signFixed == false) { 
        real_t psi_mid_init = unknowns->GetUnknownInitialData(id_psi)[0];
        real_t psi_edge_init = unknowns->GetUnknownInitialData(id_psi)[nr-1];
        if(psi_mid_init-psi_edge_init > 0)  this->sign = 1; 
        this->signFixed = true;
    }
    
    bool exists = false;
    for (len_t ir = 0; ir < nr; ir++){
        if(this->sign * PsiDiff(ir) >= 0){
            this->ir_LCFS = (int_t)ir; // Last ir to pass critera is ir_LCFS
            exists = true;
        }
    }
    if ( exists == false )  this->ir_LCFS = -1;  // No flux surface is closed
}


/**
* Difference between psi_p at ir+1/2 and the edge psi_p
*/
real_t HaloRegionHeatLossTerm::PsiDiff(len_t ir){
    // Interpolate (extrapolate for ir=0) psi for estimate of psi at inner 
    // radial grid cell wall of current grid cell
    real_t psi_f = InterpolatePsi(ir); 
    return psi_f - this->psi_edge_t0;
}


/**
* Interpolate flux between ir:s to get psi at radial grid cell walls 
* (Doesn't work for nr=1, add warning / if...else... with old condition in PsiDiff?)
*/
real_t HaloRegionHeatLossTerm::InterpolatePsi(len_t jr){

    // jr index for grid cell walls
    len_t ir_j; // Helper index, lowest ir for which to find psi 
    if (jr == 0) {
        ir_j = 0; // Extrapolate for jr=0
    } else {
        ir_j = jr-1;
    }
    if (nr==1)
		throw DREAMException("ParallelHeatLossTerm: the halo region loss term requires nr > 1");
    // Need special solution if nr=1(0)?
    real_t Rlow = this->rGrid->GetR(ir_j);
    real_t Rhigh = this->rGrid->GetR(ir_j+1);
    real_t Psilow = unknowns->GetUnknownDataPrevious(id_psi)[ir_j]; // Previous data to skip dependence of psi in Jacobian
    real_t Psihigh = unknowns->GetUnknownDataPrevious(id_psi)[ir_j+1];
    real_t slope = (Psihigh - Psilow) / (Rhigh - Rlow);
    real_t Redge = this->rGrid->GetR_f(jr);
    real_t interpPsi = Psilow + slope * (Redge - Rlow); // Extrapolates when jr=0 (then Redge < Rlow)
    
    return interpPsi;
}


/**
* Currently a step function at ir_LCFS, could be changed to continuous
*/
real_t HaloRegionHeatLossTerm::StepFunction(len_t ir){
    if ((int_t)ir > this->ir_LCFS) { 
        return 1.; 
    } else { 
        return 0.; 
    }
}


/** 
* SETWEIGHTS: Needed for DiagonalTerm, defined similar to AvalancheGrowthTerm 
*/
void HaloRegionHeatLossTerm::SetWeights() {
    FindRadiusOfLCFS(); // Find ir_LCFS
    const real_t *jtot = this->unknowns->GetUnknownData(id_jtot);
	
    real_t *T_cold = unknowns->GetUnknownData(id_T_cold); 
    real_t *N_i = unknowns->GetUnknownData(id_N_i); 
    real_t *W_i = unknowns->GetUnknownData(id_W_i); 

    for (len_t ir = 0; ir < nr; ir++) {
	   	real_t T_i = 2. / 3. * W_i[ir] / N_i[ir]; 
        real_t T_e = T_cold[ir] * Constants::ec; 
        real_t mu0Ip = Constants::mu0 * TotalPlasmaCurrentFromJTot::EvaluateIpInsideR(ir,this->grid->GetRadialGrid(),jtot);
    	real_t qR0 = this->grid->GetRadialGrid()->SafetyFactorNormalized(ir,mu0Ip);

        this->weights[ir] =  StepFunction(ir) * 2. / 3. * kappa  * sqrt((T_e + gamma * T_i) / m_i) / (M_PI * qR0); 
    }
}
/**
 * SETDIFFWEIGHTS: Needed for DiagonalComplexTerm
 */
void HaloRegionHeatLossTerm::SetDiffWeights(len_t derivId, len_t nMultiples) {
    // Retrieve necessary arrays from the unknowns handler
    real_t *T_cold = unknowns->GetUnknownData(id_T_cold);
    real_t *W_i = unknowns->GetUnknownData(id_W_i);
    real_t *N_i = unknowns->GetUnknownData(id_N_i);
    const real_t *jtot = unknowns->GetUnknownData(id_jtot);

    // Calculate the electron temperature at the radial index 'nr'
    real_t T_e = T_cold[nr] * Constants::ec;

    auto calculateMu0Ip = [&](len_t ir) {
        return Constants::mu0 * TotalPlasmaCurrentFromJTot::EvaluateIpInsideR(ir, grid->GetRadialGrid(), jtot);
    };

    auto calculateQ = [&](len_t ir) {
        return grid->GetRadialGrid()->SafetyFactorNormalized(ir, calculateMu0Ip(ir));
    };

    auto calculateWeight = [&](len_t ir, real_t qR0, real_t coeff) {


        return - StepFunction(ir) * coeff / (m_i * N_i[ir] * N_i[ir] * sqrt((T_e + gamma * 2. / 3. * W_i[ir] / N_i[ir]) / m_i)) / (2 * M_PI * qR0);
    };

    real_t coeff = (derivId == id_T_cold) ? (2. / 3.) : (4. / 9. * kappa * gamma);
    
    for (len_t i_ion = 0; i_ion <= nMultiples; i_ion++) {
        for (len_t ir = 0; ir < nr; ir++) {
            real_t qR0 = calculateQ(ir);
            if (derivId == id_N_i) {
                this->diffWeights[i_ion * nr + ir] = calculateWeight(ir, qR0, coeff * W_i[ir]);
            } else if (derivId == id_T_cold) {
                this->diffWeights[ir] = - calculateWeight(ir, qR0, coeff * W_i[ir]);
            } else if (derivId == id_W_i){
                this->diffWeights[i_ion * nr + ir] = - calculateWeight(ir, qR0, coeff);
            }
        }
    }
}



/**
* Return the weights on the radial grid
*/
const real_t* HaloRegionHeatLossTerm::GetHaloRegionHeatLossWeights(){
    return this->HeatLoss;
}


/**
* Deallocate HeatLoss
*/
void HaloRegionHeatLossTerm::Deallocate(){
    if(this->HeatLoss != nullptr) delete [] this->HeatLoss; 
}




/**
 * Heat flux parallel loss term in the plasma HALO region.
 */

#include "DREAM/DREAMException.hpp"
#include "DREAM/Equations/Fluid/ParallelHeatLossTerm.hpp"

using namespace DREAM;

/**
 * Constructor.
 */
ParallelHeatLossTerm::ParallelHeatLossTerm(
	FVM::Grid *grid, FVM::UnknownQuantityHandler *unknowns,
	FVM::Grid *operandGrid, real_t sf,
	bool userGivenPsiEdge_t0, real_t PsiEdge_t0
) : RunawaySourceTerm(grid, unknowns),
	FVM::DiagonalComplexTerm(grid, unknowns, operandGrid),
	unknowns(unknowns), scaleFactor(sf), 
    userGivenPsiEdge_t0(userGivenPsiEdge_t0), psi_edge_t0(PsiEdge_t0), 
    id_psi(unknowns->GetUnknownID(OptionConstants::UQTY_POL_FLUX)) {

    SetName("ParallelHeatLossTerm");
    this->rGrid = grid->GetRadialGrid();
    this->GridRebuilt();

    this->id_T_cold  = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->id_W_cold = unknowns->GetUnknownID(OptionConstants::UQTY_W_COLD);
    this->id_N_i = unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS);
    this->id_W_i = unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER);

    // Assuming the highest possible value for Z at the beginning (or use a specific large constant)
    len_t minIndex = -1;
    int minZ = std::numeric_limits<int>::max(); // Initialize minZ to the maximum possible integer

    real_t *n_i = x->GetUnknownData(id_N_i); // Ions density
    real_t *W_i = x->GetUnknownData(id_W_i); // Ions thermal energy
    real_t *T_i = 2. / 3. * W_i / n_i; // Ions temperature

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
ParallelHeatLossTerm::~ParallelHeatLossTerm(){
    Deallocate();

}           


/**
 * Called after the grid is rebuilt; (re)allocates memory
 * for all quantities
 */
bool ParallelHeatLossTerm::GridRebuilt(){

    this->DiagonalComplexTerm::GridRebuilt();
    Deallocate();
    this->HeatLoss = new real_t[nr];
    
    return true;
}


/**
* Set weights for term to be multiplied with unknown (n_RE). 
*/
void ParallelHeatLossTerm::SetHeatLoss(){
    
    FindRadiusOfLCFS(); // Find ir_LCFS

    real_t T_e = unknowns->GetUnknownData(id_T_cold)[nr] * Constants::ec;
    real_t W = unknowns->GetUnknownData(id_W_cold)[nr];  

    real_t *n_i = x->GetUnknownData(id_N_i)[nr]; // Ions density
    real_t *W_i = x->GetUnknownData(id_W_i)[nr]; // Ions thermal energy
    real_t *T_i = 2. / 3. * W_i / n_i; // Ions temperature  

    for (len_t ir = 0; ir < nr; ir++){
        this-> HeatLoss[ir] = StepFunction(ir) * kappa * 2. / 3. * W(ir) * sqrt((T_e(ir) + gamma * T_i(ir)) / m_i); // Change HeatLoss from 0 if ir>ir_LCFS
    }
}


/**
* Find current radial point of last closed flux surface
* and find psi_edge_t0 if not provided in input
*/
void ParallelHeatLossTerm::FindRadiusOfLCFS(){
    
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
        if(psi_mid_init-psi_edge_init > 0) { this->sign = 1; }
        this->signFixed = true;
    }
    
    bool exists = false;
    for (len_t ir = 0; ir < nr; ir++){
        if(this->sign * PsiDiff(ir) >= 0){
            this->ir_LCFS = (int_t)ir; // Last ir to pass critera is ir_LCFS
            exists = true;
        }
    }
    if ( exists == false ) { this->ir_LCFS = -1; } // No flux surface is closed
}


/**
* Difference between psi_p at ir+1/2 and the edge psi_p
*/
real_t ParallelHeatLossTerm::PsiDiff(len_t ir){
    // Interpolate (extrapolate for ir=0) psi for estimate of psi at inner 
    // radial grid cell wall of current grid cell
    real_t psi_f = InterpolatePsi(ir); 
    return psi_f - this->psi_edge_t0;
}


/**
* Interpolate flux between ir:s to get psi at radial grid cell walls 
* (Doesn't work for nr=1, add warning / if...else... with old condition in PsiDiff?)
*/
real_t ParallelHeatLossTerm::InterpolatePsi(len_t jr){

    // jr index for grid cell walls
    len_t ir_j; // Helper index, lowest ir for which to find psi 
    if (jr == 0) {
        ir_j = 0; // Extrapolate for jr=0
    } else {
        ir_j = jr-1;
    }
    
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
real_t ParallelHeatLossTerm::StepFunction(len_t ir){
    if ((int_t)ir > this->ir_LCFS) { 
        return 1.; 
    } else { 
        return 0.; 
    }
}


/** 
* SETWEIGHTS: Needed for DiagonalTerm, defined similar to AvalancheGrowthTerm 
*/
void ParallelHeatLossTerm::SetWeights() {

    SetHeatLoss(); // Set HeatLoss to be the weights on the radial grid

    Fr(nr,0,0) = 2. / 3. * W * kappa * sqrt((T_e + gamma * 2. / 3 * W_i / n_i)  m_i);
    
    len_t offset = 0;
    for (len_t ir = 0; ir < nr; ir++){
        const real_t V = this->GetVolumeScaleFactor(ir);
        
        for (len_t j = 0; j < n2[ir]; j++){
	        for (len_t i = 0; i < n1[ir]; i++){
                weights[offset + n1[ir]*j+i] = scaleFactor * this->HeatLoss[ir] * V;
                offset += n1[ir]*n2[ir];
            }
        }   
    }
}


/** 
* SETDIFFWEIGHTS: Needed for DiagonalComplexTerm: derivative of cs with respect to Te and Ti
*/
void ParalleHeatLossTerm::SetDiffWeights(len_t derivId, len_t /*nMultiples*/){
    real_t T_e = T_cold[nr] * Constants::ec;
    real_t W = unknowns->GetUnknownData(id_W_cold)[nr];  
    real_t W_i = unknowns->GetUnknownData(id_W_i)[nr];
    

    if (derivId == id_n_i) {
        dFr(nr,0,0,0) += (kappa * W) / (3 . * m_i * sqrt( (T_e + gamma * 2. / 3. * W_i / n_i) / m_i));
    } else if (derivId == id_T_cold) {
        dFr(nr,0,0,0) += (kappa * gamma * W * W_i) / (m_i * n_i^2 * sqrt( (T_e + gamma * 2. / 3. * W_i / n_i) / m_i));
    } else if (derivId == id_W_i) {
        dFr(nr,0,0,0) += (kappa * gamma * W) / (m_i * n_i * sqrt( (T_e + gamma * 2. / 3. * W_i / n_i) / m_i));
    }

}


/**
* Return the weights on the radial grid
*/
const real_t* ParallelHeatLossTerm::GetParallelHeatLossWeights(){
    return this->HeatLoss;
}


/**
* Deallocate HeatLoss
*/
void ParallelHeatLossTerm::Deallocate(){
    if(this->HeatLoss != nullptr){ delete [] this->HeatLoss; }
}




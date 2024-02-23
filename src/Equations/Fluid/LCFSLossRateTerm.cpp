/**
 * An implementation of a simple loss term modeling the loss of runaways when the flux surfaces are broken due to plasma drift.
 * A negative and linear in n_RE runaway generation rate term using a, from the input script, provided timescale t_loss. The term
 * is applied to all radial grid points outside of the last closed flux surface (LCFS). 
 */

#include "DREAM/DREAMException.hpp"
#include "DREAM/Equations/Fluid/LCFSLossRateTerm.hpp"

using namespace DREAM;

/**
 * Constructor.
 */
LCFSLossRateTerm::LCFSLossRateTerm(FVM::Grid *grid, FVM::UnknownQuantityHandler *unknowns, FVM::Grid *operandGrid, real_t sf, real_t* t_loss, len_t userGivenPsiEdge_t0, real_t PsiEdge_t0) : 
    RunawaySourceTerm(grid, unknowns), FVM::DiagonalComplexTerm(grid, unknowns, operandGrid), unknowns(unknowns), scaleFactor(sf), t_loss(t_loss), 
    userGivenPsiEdge_t0(userGivenPsiEdge_t0), psi_edge_t0(PsiEdge_t0), 
    id_psi(unknowns->GetUnknownID(OptionConstants::UQTY_POL_FLUX)) {

    SetName("LCFSLossRateTerm");
    this->rGrid = grid->GetRadialGrid();
    this->GridRebuilt();
}

/**
 * Destructor
 */
LCFSLossRateTerm::~LCFSLossRateTerm(){
    Deallocate();
}           


/**
 * Called after the grid is rebuilt; (re)allocates memory
 * for all quantities
 */
bool LCFSLossRateTerm::GridRebuilt(){

    this->DiagonalComplexTerm::GridRebuilt();
    Deallocate();
    this->GammaLoss = new real_t[nr];
    
    return true;
}


/**
* Set weights for term to be multiplied with unknown (n_RE). 
*/
void LCFSLossRateTerm::SetGammaLoss(){

    for (len_t ir = 0; ir < nr; ir++){ // Empty GammaLoss
        this->GammaLoss[ir] = 0;
    }
    
    FindRadiusOfLCFS(); // Find ir_LCFS
    
    for (len_t ir = 0; ir < nr; ir++){
        this->GammaLoss[ir] = -StepFunction(ir) / this->t_loss[ir]; // Change GammaLoss from 0 if ir>ir_LCFS
    }
}



/**
* Find current radial point of last closed flux surface
* and find psi_edge_t0 if not provided in input
*/
void LCFSLossRateTerm::FindRadiusOfLCFS(){
    
    if(this->userGivenPsiEdge_t0 == 0) {
        //this->psi_edge_t0 = unknowns->GetUnknownInitialData(id_psi)[nr-1];
        
        // Need special solution if nr=1(0)?
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
        real_t psi_mid_init = unknowns->GetUnknownInitialData(id_psi)[0]; // Use interpolate instead?
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
* Difference between the poloidal flux at radial point ir and at the edge at t=0
*/
real_t LCFSLossRateTerm::PsiDiff(len_t ir){
    //real_t psi = unknowns->GetUnknownDataPrevious(id_psi)[ir]; 
    // Interpolate (extrapolate for ir=0) psi for estimate of psi at inner 
    // radial grid cell wall of current grid cell
    real_t psi = InterpolatePsi(ir); 
    return psi - this->psi_edge_t0;
}




/**
* Interpolate flux between ir:s to get psi at radial grid cell walls 
* (Doesn't work for nr=1, add warning / if...else... with old condition in PsiDiff?)
*/
real_t LCFSLossRateTerm::InterpolatePsi(len_t jr){

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
real_t LCFSLossRateTerm::StepFunction(len_t ir){
    if ((int_t)ir > this->ir_LCFS) { 
        return 1.; 
    } else { 
        return 0.; 
    }
}



/** 
* SETWEIGHTS: Needed for DiagonalTerm, defined similar to AvalancheGrowthTerm 
*/
void LCFSLossRateTerm::SetWeights() {

    SetGammaLoss(); // Set GammaLoss to be the weights on the radial grid
    
    len_t offset = 0;
    for (len_t ir = 0; ir < nr; ir++){
        const len_t xiIndex = this->GetXiIndexForEDirection(ir);
        const real_t V = this->GetVolumeScaleFactor(ir);

        weights[offset + n1[ir]*xiIndex] = scaleFactor * this->GammaLoss[ir] * V;
        offset += n1[ir]*n2[ir];
    }
}



/** 
* SETDIFFWEIGHTS: Needed for DiagonalComplexTerm, now no derivatives for the weights but could include psi dependence here 
*/
void LCFSLossRateTerm::SetDiffWeights(len_t /*derivID*/, len_t /*nMultiples*/){}



/**
* Return the weights on the radial grid
*/
const real_t* LCFSLossRateTerm::GetLCFSLossWeights(){
    return this->GammaLoss;
}



/**
* Deallocate GammaLoss
*/
void LCFSLossRateTerm::Deallocate(){
    if(this->GammaLoss != nullptr){ delete [] this->GammaLoss; }
}







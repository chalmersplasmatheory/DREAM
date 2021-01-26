/**
 * Implementation of equation term representing the runaway generation rate
 * due to hottail when using an analytic distribution function
 */

#include "DREAM/Equations/Fluid/HottailRateTerm.hpp"

using namespace DREAM;

/**
 * Constructor.
 * 
 * gsl_altPc* contains parameters and functions needed
 * to evaluate the critical runaway momentum in the hottail
 * calculation using a gsl root-finding algorithm
 * (using the 'alternative' model for pc in Ida's MSc thesis)
 */
HottailRateTerm::HottailRateTerm(           
    FVM::Grid *grid, AnalyticDistributionHottail *dist, FVM::UnknownQuantityHandler *unknowns,
    real_t sf
) : FVM::EquationTerm(grid), RunawaySourceTerm(grid, unknowns), distHT(dist), unknowns(unknowns), scaleFactor(sf){}

/**
 * Destructor
 */
HottailRateTerm::~HottailRateTerm(){
    DeallocateAll();
}           

/**
 * Called after the grid is rebuilt; (re)allocates memory
 * for all quantities
 */
bool HottailRateTerm::GridRebuilt(){
    this->EquationTerm::GridRebuilt();

    DeallocateAll();
    pCrit      = new real_t[nr];
    gamma      = new real_t[nr];
    for(len_t ir=0; ir<nr; ir++){
        pCrit[ir] = 0;
        gamma[ir] = 0;
    }
    len_t Nmax = EquationTerm::GetMaxNumberOfMultiplesJacobian();
    dGammaDPc = new real_t[nr*Nmax];
    for(len_t i=0; i<nr*Nmax; i++)
        dGammaDPc[i] = 0;

    return true;
}

/**
 * Sets the vector elements of this equation term
 */
void HottailRateTerm::SetVectorElements(real_t *vec, const real_t *){
    len_t offset = 0;
    for(len_t ir=0; ir<nr; ir++){
        const len_t xiIndex = this->GetXiIndexForEDirection(ir);
        const len_t np1 = this->grid->GetMomentumGrid(ir)->GetNp1();
        const len_t np2 = this->grid->GetMomentumGrid(ir)->GetNp2();
        real_t V = GetVolumeScaleFactor(ir);

        // Insert at p=0, xi=+1 (or -1 if E is negative)
        vec[offset + np1*xiIndex + 0] += scaleFactor*gamma[ir]*V;

        offset += np1*np2;
    }
}

/**
 * Deallocator
 */
void HottailRateTerm::DeallocateAll(){
    if(pCrit != nullptr){
        delete [] pCrit;
        delete [] gamma;
        delete [] dGammaDPc;
    }
}

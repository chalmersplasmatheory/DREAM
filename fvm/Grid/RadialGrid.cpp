/**
 * Implementation of the radial grid.
 */

#include <algorithm>
#include <vector>
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/RadialGridGenerator.hpp"


using namespace std;
using namespace DREAM::FVM;

/***********************
 * Constructors        *
 ***********************/
/**
 * Initialize an empty grid by only specifying the
 * grid size.
 *
 * rg: Object to use for (re-)generating the radial grid.
 * t0: Time to initialize grid at.
 */
RadialGrid::RadialGrid(RadialGridGenerator *rg, const real_t t0)
    : nr(rg->GetNr()), generator(rg) {

    // Build radial grid for the first time using
    // the given RadialGridGenerator
    rg->Rebuild(t0, this);
}

/**
 * Destructor.
 */
RadialGrid::~RadialGrid() {
    // Delete radial grid quantities as usual
    DeallocateMagneticField();
    DeallocateVprime();
    DeallocateGrid();
    DeallocateFSAvg();

    delete this->generator;
}

/**
 * Deallocators.
 */
void RadialGrid::DeallocateGrid() {
    if (this->r == nullptr)
        return;

    delete [] this->dr_f;
    delete [] this->dr;
    delete [] this->r_f;
    delete [] this->r;
}

void RadialGrid::DeallocateMagneticField() {
    if (this->theta == nullptr)
        return;

    delete [] this->B_f;
    delete [] this->B;
    delete [] this->theta;
    delete [] this->Bmin;
    delete [] this->Bmin_f;
    delete [] this->Jacobian;
    delete [] this->Jacobian_f;
    
}
void RadialGrid::DeallocateVprime() {
    if (this->Vp == nullptr)
        return;

    for (len_t i = 0; i < this->nr; i++) {
        delete [] this->Vp_f2[i];
        delete [] this->Vp_f1[i];
        delete [] this->Vp[i];
    }
    for (len_t i = 0; i < this->nr+1; i++)
        delete [] this->Vp_fr[i];

    delete [] this->Vp_f2;
    delete [] this->Vp_f1;
    delete [] this->Vp_fr;
    delete [] this->Vp;
}

void RadialGrid::DeallocateFSAvg(){
    if (this->effectivePassingFraction == nullptr)
        return;
    
    for (len_t i = 0; i < this->nr; i++) {
        delete [] this->xiBounceAverage_f1[i];
        delete [] this->xiBounceAverage_f2[i];
        delete [] this->xi21MinusXi2OverB2_f1[i];
        delete [] this->xi21MinusXi2OverB2_f2[i];
    }

    delete [] this->magneticFieldMRS;
    delete [] this->effectivePassingFraction;
    delete [] this->nabla_rSq_avg;
    delete [] this->xiBounceAverage_f1;
    delete [] this->xiBounceAverage_f2;
    delete [] this->xi21MinusXi2OverB2_f1;
    delete [] this->xi21MinusXi2OverB2_f2;
}

/***************************
 * PUBLIC METHODS          *
 ***************************/
/**
 * Rebuilds any non-static (i.e. time dependent) grids
 * used. This can be used if, for example, a dynamically
 * evolving magnetic equilibrium is used, or if some
 * grids are adaptive.
 *
 * t: Time to which re-build the grids for.
 */
bool RadialGrid::Rebuild(const real_t t) {
    // Re-build radial grid
    if (this->generator->NeedsRebuild(t))
        return this->generator->Rebuild(t, this);
    else return false;
}


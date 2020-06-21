
#include <algorithm>
#include <vector>
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/RadialGridGenerator.hpp"

using namespace std;
using namespace DREAM::FVM;

RadialGridGenerator::RadialGridGenerator(const len_t nr) : nr(nr) {}


RadialGridGenerator::~RadialGridGenerator(){}


/**
 * Rebuilds magnetic field data and stores all quantities needed for flux surface and bounce averages.
 */
void RadialGridGenerator::RebuildJacobians(RadialGrid *rGrid) {
    nr = rGrid->GetNr();

    CreateMagneticFieldData(rGrid->GetR(),rGrid->GetR_f());

    rGrid->SetReferenceMagneticFieldData(
        ntheta_ref, theta_ref, B_ref, B_ref_f,
        Jacobian_ref, Jacobian_ref_f, ROverR0_ref, ROverR0_ref_f,
        NablaR2_ref, NablaR2_ref_f, 
        Bmin, Bmin_f, Bmax, Bmax_f, BtorGOverR0, BtorGOverR0_f, R0
    );
}



/**
 * Implementation of BounceSurfaceQuantity class which contains data and calculations
 * of poloidal angle-dependent quantities used in bounce averages. Contains pointer
 * to corresponding FluxSurfaceQuantity.
 */

#include "FVM/Grid/BounceSurfaceQuantity.hpp"
#include "FVM/Grid/fluxGridType.enum.hpp"


using namespace DREAM::FVM;

/**
 * Constructor
 */
BounceSurfaceQuantity::BounceSurfaceQuantity(Grid *g, FluxSurfaceQuantity *fluxSurfaceQuantity)
    : fluxSurfaceQuantity(fluxSurfaceQuantity), grid(g)
{
    gsl_acc = gsl_interp_accel_alloc();
}

/**
 * Destructor.
 */
BounceSurfaceQuantity::~BounceSurfaceQuantity(){
    DeallocateData();
    gsl_interp_accel_free(gsl_acc);
}

/**
 * Interpolates data to the trapped poloidal theta grid.
 * XXX: Assumes same momentum grid at all radii
 */
void BounceSurfaceQuantity::InterpolateToBounceGrid(
    real_t ***&bounceData, fluxGridType fluxGridType
){
    len_t nr = this->nr + (fluxGridType==FLUXGRIDTYPE_RADIAL);
    len_t n1 = np1[0] + (fluxGridType==FLUXGRIDTYPE_P1);
    len_t n2 = np2[0] + (fluxGridType==FLUXGRIDTYPE_P2);

    // XXX optimization: assume p-xi grid
    bool isPXiGrid = true;
    if(isPXiGrid){
        real_t *tmp = new real_t[ntheta_interp_trapped];
        for(len_t ir = 0; ir<nr; ir++)
            for(len_t j=0; j<n2; j++)
                if(IsTrapped(ir,0,j,fluxGridType, grid)){
                    for(len_t it=0; it<ntheta_interp_trapped; it++)
                        tmp[it] = fluxSurfaceQuantity->evaluateAtTheta(ir, ThetaBounceAtIt(ir,0,j,it,fluxGridType), fluxGridType);
                    for(len_t i=0; i<n1; i++){
                        bounceData[ir][n1*j+i] = new real_t[ntheta_interp_trapped];
                        for(len_t it=0; it<ntheta_interp_trapped; it++)
                            bounceData[ir][n1*j+i][it] = tmp[it];
                    }
                }
        delete [] tmp;
    } else {
        for(len_t ir = 0; ir<nr; ir++)
            for(len_t i=0; i<n1; i++)
                for(len_t j=0; j<n2; j++)
                    if(IsTrapped(ir,i,j,fluxGridType, grid)){
                        bounceData[ir][n1*j+i] = new real_t[ntheta_interp_trapped];
                        for(len_t it=0; it<ntheta_interp_trapped; it++)
                            bounceData[ir][n1*j+i][it] = fluxSurfaceQuantity->evaluateAtTheta(ir, ThetaBounceAtIt(ir,i,j,it,fluxGridType), fluxGridType);
                    }
    }
}

/**
 * Calculate and store data on the trapped theta grid.
 */
void BounceSurfaceQuantity::SetDataForTrapped(
    len_t ntheta_interp_trapped, real_t *quad_x_ref
){
    this->ntheta_interp_trapped = ntheta_interp_trapped;
    this->quad_x_ref = quad_x_ref;

    InterpolateToBounceGrid(bounceData, FLUXGRIDTYPE_DISTRIBUTION);
    InterpolateToBounceGrid(bounceData_fr, FLUXGRIDTYPE_RADIAL);
    InterpolateToBounceGrid(bounceData_f1, FLUXGRIDTYPE_P1);
    InterpolateToBounceGrid(bounceData_f2, FLUXGRIDTYPE_P2);
}

/**
 * Helper function to get isTrapped from Grid.
 */
bool BounceSurfaceQuantity::IsTrapped(len_t ir, len_t i, len_t j, fluxGridType fluxGridType, Grid *grid){
    switch(fluxGridType){
        case FLUXGRIDTYPE_DISTRIBUTION:
            return grid->IsTrapped(ir,i,j); // isTrapped[ir][np1[ir]*j+i];
        case FLUXGRIDTYPE_RADIAL:
            return grid->IsTrapped_fr(ir,i,j); // isTrapped_fr[ir][np1[0]*j+i];
        case FLUXGRIDTYPE_P1:
            return grid->IsTrapped_f1(ir,i,j); //isTrapped_f1[ir][(np1[ir]+1)*j+i];
        case FLUXGRIDTYPE_P2:
            return grid->IsTrapped_f2(ir,i,j); //isTrapped_f2[ir][np1[ir]*j+i];
        default:
            throw FVMException("Invalid fluxGridType: '%d' called in BounceSurfaceQuantity.", fluxGridType);
            return false;
    }
}

/**
 * Helper function to get thetaBounce1 from Grid.
 */
real_t BounceSurfaceQuantity::Theta_B1(len_t ir, len_t i, len_t j, fluxGridType fluxGridType, Grid *grid){
    if(!IsTrapped(ir,i,j,fluxGridType,grid))
        return 0;
    switch(fluxGridType){
        case FLUXGRIDTYPE_DISTRIBUTION:
            return grid->GetThetaBounce1(ir,i,j); // isTrapped[ir][np1[ir]*j+i];
        case FLUXGRIDTYPE_RADIAL:
            return grid->GetThetaBounce1_fr(ir,i,j); // isTrapped_fr[ir][np1[0]*j+i];
        case FLUXGRIDTYPE_P1:
            return grid->GetThetaBounce1_f1(ir,i,j); //isTrapped_f1[ir][(np1[ir]+1)*j+i];
        case FLUXGRIDTYPE_P2:
            return grid->GetThetaBounce1_f2(ir,i,j); //isTrapped_f2[ir][np1[ir]*j+i];
        default:
            throw FVMException("Invalid fluxGridType: '%d' called in BounceSurfaceQuantity.", fluxGridType);
            return false;
    }
}
/**
 * Helper function to get thetaBounce2 from Grid.
 */
real_t BounceSurfaceQuantity::Theta_B2(len_t ir, len_t i, len_t j, fluxGridType fluxGridType, Grid *grid){
    if(!IsTrapped(ir,i,j,fluxGridType,grid))
        return 2*M_PI;
    switch(fluxGridType){
        case FLUXGRIDTYPE_DISTRIBUTION:
            return grid->GetThetaBounce2(ir,i,j); // isTrapped[ir][np1[ir]*j+i];
        case FLUXGRIDTYPE_RADIAL:
            return grid->GetThetaBounce2_fr(ir,i,j); // isTrapped_fr[ir][np1[0]*j+i];
        case FLUXGRIDTYPE_P1:
            return grid->GetThetaBounce2_f1(ir,i,j); //isTrapped_f1[ir][(np1[ir]+1)*j+i];
        case FLUXGRIDTYPE_P2:
            return grid->GetThetaBounce2_f2(ir,i,j); //isTrapped_f2[ir][np1[ir]*j+i];
        default:
            throw FVMException("Invalid fluxGridType: '%d' called in BounceSurfaceQuantity.", fluxGridType);
            return false;
    }
}

/**
 * Helper function to get data on the trapped grid.
 */
const real_t *BounceSurfaceQuantity::GetBounceData(len_t ir, len_t i, len_t j, fluxGridType fluxGridType) const {
    switch(fluxGridType){
        case FLUXGRIDTYPE_DISTRIBUTION:
            return bounceData[ir][np1[ir]*j+i];
        case FLUXGRIDTYPE_RADIAL:
            return bounceData_fr[ir][np1[0]*j+i];
        case FLUXGRIDTYPE_P1:
            return bounceData_f1[ir][(np1[ir]+1)*j+i];
        case FLUXGRIDTYPE_P2:
            return bounceData_f2[ir][np1[ir]*j+i];
        default:
            throw FVMException("Invalid fluxGridType: '%d' called in BounceSurfaceQuantity.", fluxGridType);
            return nullptr;
    }

}

/**
 * Helper function to get stored quantity data.
 */
const real_t *BounceSurfaceQuantity::GetData(len_t ir, len_t i, len_t j, fluxGridType fluxGridType) const {
    if(IsTrapped(ir,i,j,fluxGridType,grid)){
        return GetBounceData(ir,i,j,fluxGridType);
    } else 
        return fluxSurfaceQuantity->GetData(ir, fluxGridType);
}

/**
 * Evaluates the quantity at any poloidal angle theta.
 */
const real_t BounceSurfaceQuantity::evaluateAtTheta(len_t ir, real_t theta, fluxGridType fluxGridType) const {
    return fluxSurfaceQuantity->evaluateAtTheta(ir,theta,fluxGridType);
}

/**
 * Maps the reference quadrature theta_trapped grid (defined on [0,1]) to poloidal angles.
 */
real_t BounceSurfaceQuantity::ThetaBounceAtIt(len_t ir, len_t i, len_t j, len_t it, fluxGridType fluxGridType){ 
    real_t t1 = Theta_B1(ir,i,j,fluxGridType,grid);
    real_t t2 = Theta_B2(ir,i,j,fluxGridType,grid);

	real_t tval;
	if (t1 > t2) {
		tval = (t1-2*M_PI) + (t2-t1+2*M_PI) * quad_x_ref[it];
		if (tval < 0)
			tval += 2*M_PI;
	} else
		tval = t1 + (t2-t1) * quad_x_ref[it]; 
	
	return tval;
}

/**
 * Deallocate one bounceData.
 * XXX: assumes same momentum grid at all radii
 */
void BounceSurfaceQuantity::DeleteData(real_t ***&data, len_t nr, len_t np1, len_t np2, fluxGridType fluxgrid){
    for(len_t ir=0; ir < nr; ir++){
        for (len_t j = 0; j < np2; j++)
            for(len_t i = 0; i < np1; i++) {
                bool del = false;
                switch (fluxgrid) {
                    case FLUXGRIDTYPE_DISTRIBUTION: del = grid->IsTrapped(ir, i, j); break;
                    case FLUXGRIDTYPE_RADIAL: del = grid->IsTrapped_fr(ir, i, j); break;
                    case FLUXGRIDTYPE_P1: del = grid->IsTrapped_f1(ir, i, j); break;
                    case FLUXGRIDTYPE_P2: del = grid->IsTrapped_f2(ir, i, j); break;
                    default: del = false; break;
                }

                if (del)
                    delete [] data[ir][j*np2+i];
            }
        delete [] data[ir];
    }
    delete [] data;
}

/**
 * Deallocate all trapped data
 */
void BounceSurfaceQuantity::DeallocateData(){
    if(bounceData == nullptr)
        return;
    DeleteData(bounceData,    nr,   np1[0],   np2[0], FLUXGRIDTYPE_DISTRIBUTION);
    DeleteData(bounceData_fr, nr+1, np1[0],   np2[0], FLUXGRIDTYPE_RADIAL);
    DeleteData(bounceData_f1, nr,   np1[0]+1, np2[0], FLUXGRIDTYPE_P1);
    DeleteData(bounceData_f2, nr,   np1[0],   np2[0]+1, FLUXGRIDTYPE_P2);
    
    trappedAllocated = false;
    passingAllocated = false;

}

/**
 * Allocate bounceData.
 */
void BounceSurfaceQuantity::AllocateSingle(real_t ***&bounceData, len_t nr, len_t n1, len_t n2){
    bounceData = new real_t**[nr];    
    for(len_t ir = 0; ir<nr; ir++)
        bounceData[ir] = new real_t*[n1*n2];
}

/**
 * Allocate all trapped data.
 */
void BounceSurfaceQuantity::AllocateData(){
    AllocateSingle(bounceData, nr, np1[0], np2[0]);
    AllocateSingle(bounceData_fr, nr+1, np1[0], np2[0]);
    AllocateSingle(bounceData_f1, nr, np1[0]+1, np2[0]);
    AllocateSingle(bounceData_f2, nr, np1[0], np2[0]+1);
}

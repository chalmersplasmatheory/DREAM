#include "FVM/Grid/BounceSurfaceQuantity.hpp"


using namespace DREAM::FVM;


BounceSurfaceQuantity::BounceSurfaceQuantity(Grid *g, FluxSurfaceQuantity *fluxSurfaceQuantity)
    : grid(g)  
{
    quantityData    = fluxSurfaceQuantity->GetData();
    quantityData_fr = fluxSurfaceQuantity->GetData_fr();
    quantitySpline  = fluxSurfaceQuantity->GetInterpolator();
    quantitySpline_fr = fluxSurfaceQuantity->GetInterpolator_fr();

    nr = grid->GetNr();
    np1 = new len_t[nr];
    np2 = new len_t[nr];
    for(len_t ir=0; ir<nr; ir++){
        np1[ir] = grid->GetMomentumGrid(ir)->GetNp1();
        np2[ir] = grid->GetMomentumGrid(ir)->GetNp2();
    }
    gsl_acc = gsl_interp_accel_alloc();
}

/**
 * Destructor.
 */
BounceSurfaceQuantity::~BounceSurfaceQuantity(){
    DeallocateBounceData();
    gsl_interp_accel_free(gsl_acc);

    delete [] np1;
    delete [] np2;
}


void BounceSurfaceQuantity::Initialize(
    bool **isTrapped, bool **isTrapped_fr, 
    bool **isTrapped_f1, bool **isTrapped_f2 )
{
    this->isTrapped = isTrapped;
    this->isTrapped_fr = isTrapped_fr;
    this->isTrapped_f1 = isTrapped_f1;
    this->isTrapped_f2 = isTrapped_f2;
}


/**
 * Deallocate bounceData.
 * XXX: assumes same momentum grid at all radii
 */
void DeleteData(real_t ***&data, len_t nr, len_t np1, len_t np2){
    for(len_t ir=0; ir<nr; ir++){
        for(len_t i = 0; i<np1*np2; i++)
            delete [] data[ir][i];
        delete [] data[ir];
    }
    delete [] data;
}
void BounceSurfaceQuantity::DeallocateBounceData(){
    if(bounceData==nullptr)
        return;
    DeleteData(bounceData,    nr,   np1[0],   np2[0]);
    DeleteData(bounceData_fr, nr+1, np1[0],   np2[0]);
    DeleteData(bounceData_f1, nr,   np1[0]+1, np2[0]);
    DeleteData(bounceData_f2, nr,   np1[0],   np2[0]+1);

}


/**
 * Interpolates data to the poloidal bouncegrid bounceTheta.
 */
void InterpolateToBounceGrid(
    real_t ***&bounceData, real_t ***bounceTheta, FluxSurfaceQuantity *fluxSurfaceQuantity, bool **isTrapped,
    len_t nr, len_t np1, len_t np2, len_t ntheta, gsl_interp_accel *gsl_acc, fluxGridType fluxGridType){
    bounceData = new real_t**[nr];
    for(len_t ir = 0; ir<nr; ir++){
        bounceData[ir] = new real_t*[np1*np2];
        for(len_t i=0; i<np1*np2; i++){
            bounceData[ir][i] = new real_t[ntheta];
            if(isTrapped[ir][i])
                for(len_t it=0; it<ntheta; it++)
                    bounceData[ir][i][it] = fluxSurfaceQuantity->evaluateAtTheta(ir, bounceTheta[ir][i][it], fluxGridType);
        }
    }
}
void BounceSurfaceQuantity::InterpolateMagneticDataToBounceGrids(
    len_t ntheta_interp_trapped, real_t ***bounceTheta, real_t ***bounceTheta_fr,
     real_t ***bounceTheta_f1, real_t ***bounceTheta_f2
){
    DeallocateBounceData();

    InterpolateToBounceGrid(
        bounceData, bounceTheta, fluxSurfaceQuantity, isTrapped,
        nr, np1[0], np2[0], ntheta_interp_trapped, gsl_acc, FLUXGRIDTYPE_DISTRIBUTION);
    InterpolateToBounceGrid(
        bounceData_fr, bounceTheta_fr, fluxSurfaceQuantity, isTrapped_fr,
        nr+1, np1[0], np2[0], ntheta_interp_trapped, gsl_acc, FLUXGRIDTYPE_RADIAL);
    InterpolateToBounceGrid(
        bounceData_f1, bounceTheta_f1, fluxSurfaceQuantity, isTrapped_f1,
        nr, np1[0]+1, np2[0], ntheta_interp_trapped, gsl_acc, FLUXGRIDTYPE_P1);
    InterpolateToBounceGrid(
        bounceData_f2, bounceTheta_f2, fluxSurfaceQuantity, isTrapped_f2,
        nr, np1[0], np2[0]+1, ntheta_interp_trapped, gsl_acc, FLUXGRIDTYPE_P2);
}


const bool BounceSurfaceQuantity::IsTrapped(len_t ir, len_t i, len_t j, fluxGridType fluxGridType) const {
    switch(fluxGridType){
        case FLUXGRIDTYPE_DISTRIBUTION:
            return isTrapped[ir][np1[ir]*j+i];
        case FLUXGRIDTYPE_RADIAL:
            return isTrapped_fr[ir][np1[0]*j+i];
        case FLUXGRIDTYPE_P1:
            return isTrapped_f1[ir][(np1[ir]+1)*j+i];
        case FLUXGRIDTYPE_P2:
            return isTrapped_f2[ir][np1[ir]*j+i];
        default:
            throw FVMException("Invalid fluxGridType: '%d' called in BounceSurfaceQuantity.", fluxGridType);
            return false;
    }
}
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
const real_t BounceSurfaceQuantity::evaluateAtTheta(len_t ir, real_t theta, fluxGridType fluxGridType) const {
    return fluxSurfaceQuantity->evaluateAtTheta(ir,theta,fluxGridType);
}

const real_t *BounceSurfaceQuantity::GetData(len_t ir, len_t i, len_t j, fluxGridType fluxGridType) const {
    if(IsTrapped(ir,i,j,fluxGridType)){
        return GetBounceData(ir,i,j,fluxGridType);
    } else 
        return fluxSurfaceQuantity->GetData(ir, fluxGridType);
}

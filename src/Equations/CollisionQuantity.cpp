/**
 * Implementation of a class from which the Coulomb logarithms and collision frequencies are derived.
 */

#include "DREAM/Equations/CollisionQuantity.hpp"
//#include "DREAM/Constants.hpp"
//#include "DREAM/Settings/OptionConstants.hpp"
//#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/NotImplementedException.hpp"

using namespace DREAM;

/**
 * Constructor.
 */
CollisionQuantity::CollisionQuantity(FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ih,  
                enum OptionConstants::momentumgrid_type mgtype,  struct collqty_settings *cqset){
    mg = g->GetMomentumGrid(0);
    rGrid = g->GetRadialGrid();

    collQtySettings = cqset;
    ionHandler = ih;
    unknowns = u;

    // Various settings that appear in many places in the calculations
    isPXiGrid = (mgtype==OptionConstants::MOMENTUMGRID_TYPE_PXI);
    isPartiallyScreened = (collQtySettings->collfreq_type==OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_PARTIALLY_SCREENED);
    isNonScreened = (collQtySettings->collfreq_type==OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_NON_SCREENED);
    isNonlinear = (collQtySettings->nonlinear_mode == OptionConstants::EQTERM_NONLINEAR_MODE_NON_REL_ISOTROPIC);
    isBrems = (collQtySettings->bremsstrahlung_mode != OptionConstants::EQTERM_BREMSSTRAHLUNG_MODE_NEGLECT);
    // ID of quantities that contribute to collision frequencies
    id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    id_ni    = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);

    /**
     * Set buildOnlyF1F2=false if quantities need to be evaluated on the distribution 
     * and radial flux grids. For now hardcoded to true because it isn't expected to 
     * be needed. In fact the calculation on the radial flux grid is not supported,
     * since we do not yet interpolate in the unknown quantities.
     */
    buildOnlyF1F2 = true;

    /**
     * This is the ad-hoc k-parameter appearing in Linneas paper that sets the transition
     * from superthermal limits to p->0. For now hardcoded to the recommended value k=5.
     */
    kInterpolate = 5;
}

/**
 * Destructor.
 */
CollisionQuantity::~CollisionQuantity(){
    DeallocateCollisionQuantities();
}

/**
 * Rebuilds collision quantities; the quantities are split into 
 * multiple partial contributions, and only those that change are
 * rebuilt. 
 */
void CollisionQuantity::Rebuild(){
    if (gridRebuilt){ // Reallocate and rebuild everything
        nr  = rGrid->GetNr();
        nZ  = ionHandler->GetNZ();
        nzs = ionHandler->GetNzs();
        np1 = mg->GetNp1();
        np2 = mg->GetNp2();
        if(isPXiGrid)
            np2_store = 1;
        else
            np2_store = mg->GetNp2();
        // We should Deallocate before we update nr, nZ etc. so that deletions occurs with parameters of previous iteration
        AllocateCollisionQuantities(); 
        
        RebuildConstantTerms();
        RebuildPlasmaDependentTerms();
        AssembleQuantity();
        gridRebuilt = false;
    } else if (parametersHaveChanged()){
        RebuildPlasmaDependentTerms();
        AssembleQuantity();
    }    
}


/** 
 * Assembles the collision frequency from the various partial contributions, on all grids.
 */
void CollisionQuantity::AssembleQuantity(){
    if(!buildOnlyF1F2){
        AssembleQuantity(collisionQuantity,nr,np1,np2,FVM::FLUXGRIDTYPE_DISTRIBUTION);
        AssembleQuantity(collisionQuantity_fr,nr /* +1 this case should be treated properly somehow*/,np1,np2,FVM::FLUXGRIDTYPE_RADIAL);
    }
    AssembleQuantity(collisionQuantity_f1,nr,np1+1,np2,FVM::FLUXGRIDTYPE_P1);
    AssembleQuantity(collisionQuantity_f2,nr,np1,np2+1,FVM::FLUXGRIDTYPE_P2);

}

/** 
 * Returns true if any unknown quantities that affect collision quantities have changed. 
 */
bool CollisionQuantity::parametersHaveChanged(){
    return unknowns->HasChanged(id_ncold) || unknowns->HasChanged(id_Tcold) || unknowns->HasChanged(id_ni);
}


void CollisionQuantity::AllocateCollisionQuantities(){
    DeallocateCollisionQuantities();

    if(!buildOnlyF1F2){
        AllocateCollisionQuantity(collisionQuantity,nr,np1,np2);
        AllocateCollisionQuantity(collisionQuantity_fr,nr+1,np1,np2);
    }
    AllocateCollisionQuantity(collisionQuantity_f1,nr,np1+1,np2);
    AllocateCollisionQuantity(collisionQuantity_f2,nr,np1,np2+1);


    AllocatePartialQuantities();
}
void CollisionQuantity::AllocateCollisionQuantity(real_t **&collisionQuantity, len_t nr, len_t np1, len_t np2){
    collisionQuantity = new real_t*[nr];
    for(len_t ir = 0; ir<nr; ir++){
        collisionQuantity[ir] = new real_t[np1*np2]; 
    }
}

void CollisionQuantity::DeallocateCollisionQuantity(real_t **&collisionQuantity, len_t nr){
    if(collisionQuantity != nullptr){
        for(len_t ir = 0; ir<nr; ir++)
            delete [] collisionQuantity[ir];
        delete [] collisionQuantity;
    }
}

void CollisionQuantity::DeallocateCollisionQuantities(){
    if(!buildOnlyF1F2){
        DeallocateCollisionQuantity(collisionQuantity,nr);
        DeallocateCollisionQuantity(collisionQuantity_fr,nr+1);
    }
    DeallocateCollisionQuantity(collisionQuantity_f1,nr);
    DeallocateCollisionQuantity(collisionQuantity_f2,nr);

    
}







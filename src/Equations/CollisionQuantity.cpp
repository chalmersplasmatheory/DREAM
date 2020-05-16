
#include "DREAM/Equations/CollisionQuantity.hpp"
#include "DREAM/Constants.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/NotImplementedException.hpp"

using namespace DREAM;

CollisionQuantity::CollisionQuantity(FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ih,  
                enum OptionConstants::momentumgrid_type mgtype,  struct CollisionQuantityHandler::collqtyhand_settings *cqset){
    mg = g->GetMomentumGrid(0);
    rGrid = g->GetRadialGrid();

    collQtySettings = cqset;
    ionHandler = ih;
    unknowns = u;

    isPXiGrid = (mgtype==OptionConstants::MOMENTUMGRID_TYPE_PXI);
    isPartiallyScreened = (collQtySettings->collfreq_type==OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_PARTIALLY_SCREENED);
    isNonScreened = (collQtySettings->collfreq_type==OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_NON_SCREENED);
    isNonlinear = (collQtySettings->nonlinear_mode == OptionConstants::EQTERM_NONLINEAR_MODE_NON_REL_ISOTROPIC);

    // ID of quantities that contribute to collision frequencies
    id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    id_ni    = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    id_fhot = unknowns->GetUnknownID(OptionConstants::UQTY_F_HOT);
}


CollisionQuantity::~CollisionQuantity(){
    DeallocateCollisionQuantities();
    
}

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
        AllocateCollisionQuantities();
        
        RebuildConstantTerms();
        RebuildPlasmaDependentTerms();
    } else if (parametersHaveChanged()){
        RebuildPlasmaDependentTerms();
    }    
    AssembleQuantity();

    gridRebuilt = false;
}


/** 
 * Assembles the collision frequency from the various partial contributions, on all grids.
 */
void CollisionQuantity::AssembleQuantity(){
    if(!buildOnlyF1F2){
        AssembleQuantity(collisionQuantity,nr,np1,np2,0);
        AssembleQuantity(collisionQuantity_fr,nr+1,np1,np2,1);
    }
    AssembleQuantity(collisionQuantity_f1,nr,np1+1,np2,2);
    AssembleQuantity(collisionQuantity_f2,nr,np1,np2+1,3);

}


bool CollisionQuantity::parametersHaveChanged(){
    if(unknowns->HasChanged(id_ncold) || unknowns->HasChanged(id_Tcold) || unknowns->HasChanged(id_ni))
        return true;
    else
        return false;    
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







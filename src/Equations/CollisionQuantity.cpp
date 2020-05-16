
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

    isPXiGrid = (mgtype==OptionConstants::MOMENTUMGRID_TYPE_PXI);
    collQtySettings = cqset;
    ionHandler = ih;
    unknowns = u;

    isPartiallyScreened = (collQtySettings->collfreq_type==OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_PARTIALLY_SCREENED);
    isNonScreened = (collQtySettings->collfreq_type==OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_NON_SCREENED);
    isNonlinear = (collQtySettings->nonlinear_mode == OptionConstants::EQTERM_NONLINEAR_MODE_NON_REL_ISOTROPIC);

    id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    id_ni    = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    id_fhot = unknowns->GetUnknownID(OptionConstants::UQTY_F_HOT);
//    AllocateCollisionQuantities();    
}


CollisionQuantity::~CollisionQuantity(){
    DeallocateCollisionQuantities();
    DeallocateGSL();
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
        if(collQtySettings->collfreq_mode==OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL)
            InitializeGSLWorkspace();
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

/**
 * Assembles collision frequency on one of the grids.
 */
void CollisionQuantity::AssembleQuantity(real_t **&collisionQuantity, len_t nr, len_t np1, len_t np2, len_t fluxGridType){
    real_t *nColdContribution = new real_t[1];
    real_t *niContribution = new real_t[nzs];
    real_t collQty;
    real_t *ncold = unknowns->GetUnknownData(id_ncold);
    const len_t *Zs = ionHandler->GetZs();

    len_t indZ;
    for(len_t ir=0; ir<nr; ir++){
        for(len_t i=0; i<np1; i++){
            for(len_t j=0; j<np2; j++){
                nColdContribution = GetUnknownPartialContribution(id_ncold,ir,i,j,fluxGridType,nColdContribution);
                niContribution = GetUnknownPartialContribution(id_ni,ir,i,j,fluxGridType,niContribution);
                collQty = ncold[ir]*nColdContribution[0];
                for(len_t iz = 0; iz<nZ; iz++){
                    for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                        indZ = ionHandler->GetIndex(iz,Z0);
                        collQty += ionHandler->GetIonDensity(ir,iz,Z0)*niContribution[indZ];
                    }
                }
                collisionQuantity[ir][j*np1+i] = collQty; 
            }
        }
    }
    delete [] nColdContribution;
    delete [] niContribution;
}


void CollisionQuantity::AddNonlinearContribution(const real_t *lnLc){
    real_t *fHot = unknowns->GetUnknownData(id_fhot);
    for (len_t ir=0;ir<nr;ir++)
        for(len_t i=0; i<np1; i++)
            for(len_t j=0; j<np2; j++)
                for(len_t ip=0; ip<np1; ip++)
                    collisionQuantity[ir][np1*j+i] += lnLc[ir]*nonlinearMat[i][ip]*fHot[np1*ir + ip];
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


    if (isNonlinear){
        nonlinearMat = new real_t*[np1+1]; // multiply matrix by f lnLc to get p*nu_s on p flux grid
        for (len_t i = 0; i<np1+1; i++){
            nonlinearMat[i] = new real_t[np1];
        }
        const real_t *p = mg->GetP1();
        trapzWeights = new real_t[np1];
        for (len_t i = 1; i<np1-1; i++){
            trapzWeights[i] = (p[i+1]-p[i-1])/2;
        }
    }
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

    if(nonlinearMat != nullptr){
        for(len_t i = 0; i<np1+1;i++){
            delete [] nonlinearMat[i];
        }
        delete [] nonlinearMat;
        delete [] trapzWeights;
    }
}









real_t CollisionQuantity::psi0Integrand(real_t x, void *params){
    real_t gamma = *(real_t *) params;
    return 1/sqrt( (x+gamma)*(x+gamma)-1 );
} 
real_t CollisionQuantity::psi1Integrand(real_t x, void *params){
    real_t gamma = *(real_t *) params;
    return (x+gamma)/sqrt((x+gamma)*(x+gamma)-1); // integrated with weight w(x) = exp(-(x-gamma)/Theta) 
} 
/** 
 * Evaluates integral appearing in relativistic test-particle operator
 * Psi0 = int_0^p exp( -(sqrt(1+s^2)-1)/Theta) / sqrt(1+s^2) ds;
 */
real_t CollisionQuantity::evaluatePsi0(len_t ir, real_t p) {
    real_t gamma = sqrt(1+p*p);
    real_t *T_cold = unknowns->GetUnknownData(id_Tcold);
    
    gsl_function F;
    F.function = &(CollisionQuantity::psi0Integrand); 
    F.params = &gamma;
    real_t psi0int; 
    gsl_integration_fixed(&F, &psi0int, gsl_w[ir]);

    
    real_t Theta = T_cold[ir] / Constants::mc2inEV;
    return evaluateExp1OverThetaK(Theta,0) - exp( -(gamma-1)/Theta ) * psi0int;

}
real_t CollisionQuantity::evaluatePsi1(len_t ir, real_t p) {
    
    real_t *T_cold = unknowns->GetUnknownData(id_Tcold);
    real_t gamma = sqrt(1+p*p);
    gsl_function F;
    F.function = &(CollisionQuantity::psi1Integrand); 
    F.params = &gamma;
    real_t psi1int; 
    gsl_integration_fixed(&F, &psi1int, gsl_w[ir]);

    real_t Theta = T_cold[ir] / Constants::mc2inEV;
    return evaluateExp1OverThetaK(Theta,1) - exp( -(gamma-1)/Theta ) * psi1int;
    

}


real_t CollisionQuantity::evaluateExp1OverThetaK(real_t Theta, real_t n) {
    real_t ThetaThreshold = 0.002;
    /**
     * Since cyl_bessel_k ~ exp(-1/Theta), for small Theta you get precision issues.
     * Instead using asymptotic expansion for bessel_k for small Theta.
     */
    if (Theta > ThetaThreshold)
        return exp(1/Theta)*std::cyl_bessel_k(n,1/Theta);
    else {
//        return sqrt(M_PI*Theta/2)*(1 + 15*Theta/8 + 105*Theta*Theta/128 - 945*Theta*Theta*Theta/3072);
        real_t n2 = n*n;
        return sqrt(M_PI*Theta/2)*(1 + (4*n2-1)/8 * Theta + (4*n2-1)*(4*n2-9)*Theta*Theta/128 + (4*n2-1)*(4*n2-9)*(4*n2-25)*Theta*Theta*Theta/3072);
    }
}


/**
 * Initializes a GSL workspace for each radius (used for relativistic test particle operator evaluation),
 * using a T_cold-dependent fixed quadrature. 
 */
void CollisionQuantity::InitializeGSLWorkspace(){
 /** 
  * (consider using a single regular dynamic quadrature instead as the integral is somewhat tricky, 
  * since in the limit p/mc -> 0 the integral is sharply peaked at p_min -- goes as int 1/sqrt(x) dx,0,inf --
  * and may be challenging to resolve using a fixed point quadrature)
  */

    real_t *T_cold = unknowns->GetUnknownData(id_Tcold);
    DeallocateGSL();
    gsl_w = new gsl_integration_fixed_workspace*[nr];
    const real_t lowerLim = 0; // integrate from 0 to inf
    const gsl_integration_fixed_type *T = gsl_integration_fixed_laguerre;
    const len_t Npoints = 20; // play around with this number -- may require larger, or even sufficient with lower
    const real_t alpha = 0.0;
    real_t b;
    real_t Theta;
    for (len_t ir = 0; ir<nr; ir++){
        Theta = T_cold[ir]/Constants::mc2inEV;
        b = 1/Theta;
        gsl_w[ir] = gsl_integration_fixed_alloc(T, Npoints, lowerLim, b, alpha, 0.0);
    }
}


void CollisionQuantity::DeallocateGSL(){
    if (this->gsl_w == nullptr)
        return;

    for (len_t ir=0; ir<this->nr; ir++)
        gsl_integration_fixed_free(gsl_w[ir]);
}


real_t *CollisionQuantity::GetUnknownPartialContribution(len_t id_unknown, len_t ir, len_t i, len_t j, len_t fluxGridMode, real_t *&partQty){
    if(fluxGridMode == 0){
        GetPartialContribution(id_unknown, ir, i, j, partQty);
    } else if(fluxGridMode == 1){
        GetPartialContribution_fr( id_unknown,  ir,  i, j,partQty);
    } else if(fluxGridMode == 2){
        GetPartialContribution_f1(id_unknown,  ir, i, j, partQty);
    } else if(fluxGridMode == 3){
        GetPartialContribution_f2(id_unknown, ir, i, j, partQty);
    } else 
        throw FVM::FVMException("Invalid fluxGridMode: Set 0 (distribution grid), 1 (radial flux grid), 2 (p1 flux grid) or 3 (p2 flux grid)");
    
    return partQty;
}


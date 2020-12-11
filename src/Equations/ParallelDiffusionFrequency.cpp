/**
 * Implementation of a class which handles the calculation of the 
 * parallel diffusion coefficient nu_||. The D^pp component of the 
 * collision operator is defined as (m_e c)^2 nu_|| (and not the 
 * definition that one often encounters in the literature, where it 
 * has been multiplied by p^2/2).
 * It should be rebuilt after the SlowingDownFrequency has been rebuilt.
 */

/**
 * For the linearised case with a known temperature, it is uniquely prescribed by nu_s,
 * by requiring the preservation of the Maxwell-Jüttner steady-state distribution. 
 * The method is described in Appendix A.1 of doc/notes/theory.pdf.
 * The non-linear contribution corresponds to the isotropic component of the
 * non-relativistic operator following Rosenbluth, Macdonald & Judd, Phys Rev (1957),
 * and is described in doc/notes/theory.pdf Appendix B.
 */

// TODO: Implement real_t *GetUnknownPartialContribution(len_t id_unknown,len_t fluxGridMode, real_t *&partQty);

#include "DREAM/Equations/ParallelDiffusionFrequency.hpp"
#include "DREAM/NotImplementedException.hpp"

using namespace DREAM;


/**
 * Constructor.
 */
ParallelDiffusionFrequency::ParallelDiffusionFrequency(FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ih,  
                SlowingDownFrequency *nuS, CoulombLogarithm *lnLee,
                enum OptionConstants::momentumgrid_type mgtype,  struct collqty_settings *cqset)
                : CollisionQuantity(g,u,ih,mgtype,cqset){
    this->lnLambdaEE = lnLee;
    this->nuS = nuS;
    includeDiffusion = (cqset->collfreq_mode==OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL);
}

/**
 * Destructor.
 */
ParallelDiffusionFrequency::~ParallelDiffusionFrequency(){
    DeallocateCollisionQuantities();    
    DeallocatePartialQuantities();      
}

/**
 * Calculates the parallel diffusion frequency from the values of the slowing-down frequency.
 */
void ParallelDiffusionFrequency::AssembleQuantity(real_t **&collisionQuantity, len_t nr, len_t np1, len_t np2, FVM::fluxGridType fluxGridType){
    if(!includeDiffusion){
        for(len_t ir=0;ir<nr;ir++)
            for(len_t it=0;it<np1*np2;it++)
                collisionQuantity[ir][it] = 0;
        return;
    }

    const real_t *gammaVec = mg->GetGamma(fluxGridType);
    if(collQtySettings->screened_diffusion==OptionConstants::COLLQTY_SCREENED_DIFFUSION_MODE_MAXWELLIAN){
        real_t *const* nuSQty = nuS->GetValue(fluxGridType);
        for(len_t ir=0; ir<nr; ir++)
            for(len_t i=0; i<np1; i++)
                for(len_t j=0; j<np2; j++){
                    len_t pind = np1*j+i;
                    collisionQuantity[ir][pind] = nuSQty[ir][pind] * rescaleFactor(ir,gammaVec[pind]);
                }
    } else if(collQtySettings->screened_diffusion==OptionConstants::COLLQTY_SCREENED_DIFFUSION_MODE_ZERO){
        const real_t *partialNuS = nuS->GetUnknownPartialContribution(id_ncold,fluxGridType);
        const real_t *ncold = unknowns->GetUnknownData(id_ncold);
        len_t offset = 0;
        for(len_t ir=0; ir<nr; ir++){
            for(len_t i=0; i<np1; i++)
                for(len_t j=0; j<np2; j++){
                    len_t pind = np1*j+i;
                    collisionQuantity[ir][pind] = ncold[ir]*partialNuS[offset + pind] * rescaleFactor(ir,gammaVec[pind]);
                }
            offset += np1*np2;
        }
    }
        
    if(isNonlinear && (fluxGridType == FVM::FLUXGRIDTYPE_P1))
        SetNonlinearPartialContribution(lnLambdaEE, fHotPartialContribution_f1);
}


/**
 * Allocator
 */
void ParallelDiffusionFrequency::AllocatePartialQuantities(){
    DeallocatePartialQuantities();
    Theta = new real_t[nr];
    partContrib = new real_t[nzs*(nr+1)*(np1+1)*(np2+1)];
    if (isNonlinear){
        nonlinearMat = new real_t*[np1+1]; 
        for (len_t i = 0; i<np1+1; i++)
            nonlinearMat[i] = new real_t[np1];
        const real_t *p = mg->GetP1();
        trapzWeights = new real_t[np1];
        for (len_t i = 1; i<np1-1; i++)
            trapzWeights[i] = (p[i+1]-p[i-1])/2;
        fHotPartialContribution_f1 = new real_t[nr*np1*(np1+1)];
    }    
}


/**
 * Deallocator
 */
void ParallelDiffusionFrequency::DeallocatePartialQuantities(){
    if(Theta !=nullptr)
        delete [] Theta;

    if(partContrib != nullptr)
        delete [] partContrib;

    if(nonlinearMat != nullptr){
        for(len_t i = 0; i<np1+1;i++)
            delete [] nonlinearMat[i];
        delete [] nonlinearMat;
        delete [] trapzWeights;
        delete [] fHotPartialContribution_f1;
    }

}

/**
 * Rebuilds the part of the calculation that depends on plasma parameters.
 */
void ParallelDiffusionFrequency::RebuildPlasmaDependentTerms(){
    if(!includeDiffusion)
        return;
    real_t *Tcold = unknowns->GetUnknownData(id_Tcold);
    for(len_t ir=0; ir<nr;ir++)
        Theta[ir] = Tcold[ir]/Constants::mc2inEV;
}

/**
 * The factor by which the slowing-down frequency should be multiplied in order
 * to yield the parallel diffusion frequency.
 */
real_t ParallelDiffusionFrequency::rescaleFactor(len_t ir, real_t gamma){
    return Theta[ir]*gamma;
}

/**
 * Evaluates the frequency at radial grid point ir and momentum p.
 */
real_t ParallelDiffusionFrequency::evaluateAtP(len_t ir, real_t p, struct collqty_settings *inSettings){
    if(inSettings->collfreq_mode != OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL)
        return 0;
    return rescaleFactor(ir,sqrt(1+p*p))*nuS->evaluateAtP(ir,p,inSettings);
}


/**
 * Evaluates the jacobian of the frequency with respect to 
 * the unknown derivId at radial grid point ir and momentum p
 */
real_t ParallelDiffusionFrequency::evaluatePartialAtP(len_t ir, real_t p, len_t derivId, len_t n,struct collqty_settings *inSettings){
    if(inSettings->collfreq_mode != OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL)
        return 0;
    real_t dnuPar0 = rescaleFactor(ir,sqrt(1+p*p))*nuS->evaluatePartialAtP(ir,p,derivId,n,inSettings);
    if(derivId==id_Tcold){ // and contribution from derivative of rescaleFactor \propto Tcold
        real_t *Tcold = unknowns->GetUnknownData(id_Tcold);
        dnuPar0 += evaluateAtP(ir,p,inSettings) / Tcold[ir];
    }
    return dnuPar0;
}



/** Adds the non-linear contribution to the collision frequency. For now, only supports 
 * hot-tail grids where np2=1 and using a pxi-grid, and only updates the p flux grid 
 * component.
 */
void ParallelDiffusionFrequency::AddNonlinearContribution(){
    real_t *fHot = unknowns->GetUnknownData(OptionConstants::UQTY_F_HOT);
    const real_t* const fHotPartialContribution_f1 = GetNonlinearPartialContribution(FVM::FLUXGRIDTYPE_P1);

    for (len_t ir=0;ir<nr;ir++)
        for(len_t i=0; i<np1+1; i++)
            for(len_t ip=0; ip<np1; ip++)
                collisionQuantity_f1[ir][i] += fHotPartialContribution_f1[ip*(np1+1)*nr + ir*(np1+1) + i] * fHot[np1*ir+ip];
}

void ParallelDiffusionFrequency::SetNonlinearPartialContribution(CoulombLogarithm *lnLambda, real_t *&partQty){
    if(partQty==nullptr)
        partQty = new real_t[np1*(np1+1)*nr];
    for(len_t it=0; it < np1*(np1+1)*nr; it++)
        partQty[it] = 0;

    for(len_t i=0; i<np1+1; i++)
        for(len_t ir=0;ir<nr;ir++)
            for(len_t ip=0; ip<np1; ip++)
                partQty[ip*(np1+1)*nr + (np1+1)*ir + i] = lnLambda->GetLnLambdaT(ir)*nonlinearMat[i][ip];
}

const real_t* ParallelDiffusionFrequency::GetNonlinearPartialContribution(FVM::fluxGridType fluxGridType) const{
    if(fluxGridType==FVM::FLUXGRIDTYPE_P1)
        return fHotPartialContribution_f1;
    else {
//        throw FVM::FVMException("Invalid fluxGridType. Nonlinear contribution only supported for p1 flux grid.");
        return nullptr;
    }
}


/**
 * Rebuilds partial contributions that only depend on the grid.
 */
void ParallelDiffusionFrequency::RebuildConstantTerms(){
    if(isNonlinear)
        calculateIsotropicNonlinearOperatorMatrix();
}

/**
 * Calculates a Rosenbluth potential matrix defined such that when it is muliplied
 * by the f_hot distribution vector, yields the parallel diffusion frequency.
 */
void ParallelDiffusionFrequency::calculateIsotropicNonlinearOperatorMatrix(){

    if( !(isPXiGrid && (mg->GetNp2() == 1)) )
        throw NotImplementedException("Nonlinear collisions only implemented for hot tails (np2=1) and p-xi grid");

    
    const real_t *p_f = mg->GetP1_f();
    const real_t *p = mg->GetP1();

    // See doc/notes/theory.pdf appendix B for details on discretization of integrals;
    // uses a trapezoidal rule
    real_t p2, p2f;
    real_t weightsIm1, weightsI;
    for (len_t i = 1; i<np1+1; i++){
        p2f = p_f[i]*p_f[i];
        p2 = p[0]*p[0];
        nonlinearMat[i][0]   = (4*M_PI/3) * constPreFactor*( (p[1]-p[0])/2  + p[0]/5 )* p2*p2/(p_f[i]*p2f);
        for (len_t ip = 1; ip < i-1; ip++){
            p2 = p[ip]*p[ip];
            nonlinearMat[i][ip]   = (4*M_PI/3) * constPreFactor* trapzWeights[ip]*p2*p2 / (p_f[i]*p2f);
        } 
        p2 = p[i-1]*p[i-1];
        weightsIm1 = (p[i-1]-p[i-2])/2 + (p_f[i]-p[i-1])/(p[i]-p[i-1])*( (2*p[i]-p_f[i]-p[i-1])/2 );
        nonlinearMat[i][i-1]   = (4*M_PI/3) * constPreFactor * weightsIm1*p2*p2 / (p_f[i]*p2f);
        p2 = p[i]*p[i];
        weightsI = (p_f[i]-p[i-1])*(p_f[i]-p[i-1])/(p[i]-p[i-1]);
        nonlinearMat[i][i]   = (4*M_PI/3) * constPreFactor * weightsI*p2*p2 / (p_f[i]*p2f);

        // add contribution from p'>p terms near p'=p
        p2 = p[i-1]*p[i-1];
        weightsIm1 = (1.0/2)*(p[i]-p_f[i])*(p[i]-p_f[i])/(p[i]-p[i-1]);
        nonlinearMat[i][i-1]   += (4*M_PI/3) * constPreFactor * weightsIm1*p[i-1];
        p2 = p[i]*p[i];
        weightsI = (p[i+1]-p[i])/2 + (1.0/2)*(p[i]-p_f[i])*(p_f[i]+p[i]-2*p[i-1])/(p[i]-p[i-1]);
        nonlinearMat[i][i]   += (4*M_PI/3) * constPreFactor * weightsI * p[i];


        for (len_t ip = i+1; ip < np1-1; ip++){
            nonlinearMat[i][ip] = (4*M_PI/3) * constPreFactor* trapzWeights[ip]*p[ip];
        } 
        real_t weightsEnd = (p[np1-1]-p[np1-2])/2;
        nonlinearMat[i][np1-1] = (4*M_PI/3) * constPreFactor* weightsEnd*p[np1-1];
        
    }


}




/**
 * Calculation of the partial contribution to the collision frequency from the unknown quantity
 * with ID id_unknown. Returns the partial derivative of the term with respect to that quantity 
 * (ignoring variations with lnLambda). 
 */
const real_t* ParallelDiffusionFrequency::GetUnknownPartialContribution(len_t id_unknown, FVM::fluxGridType fluxGridType){
    if( (id_unknown == id_ncold) || (id_unknown == id_ni) || (id_unknown == id_Tcold)){ // for ncold and ni, simply rescale the values from nuS
        const real_t *gammaVec = mg->GetGamma(fluxGridType);
        
        len_t nr  = this->nr  + (fluxGridType == FVM::FLUXGRIDTYPE_RADIAL);
        len_t np1 = this->np1 + (fluxGridType == FVM::FLUXGRIDTYPE_P1);
        len_t np2 = this->np2 + (fluxGridType == FVM::FLUXGRIDTYPE_P2);
    
        len_t numZs = 1; // number of multiples
        if(id_unknown == id_ni)
            numZs = nzs;

        const real_t *partContribNuS = nuS->GetUnknownPartialContribution(id_unknown,fluxGridType);
        //real_t *partContrib = new real_t[numZs*nr*np1*np2];
        for(len_t i=0; i<nzs*(this->nr+1)*(this->np1+1)*(this->np2+1); i++) // reset entire partContrib array
            partContrib[i] = 0;
        if(!includeDiffusion)
            return partContrib;

        for(len_t pind=0; pind<np1*np2; pind++)
            for(len_t ir=0; ir<nr; ir++){
                real_t rescaleFact =  rescaleFactor(ir,gammaVec[pind]);
                for(len_t indZ = 0; indZ<numZs; indZ++){
                    len_t i = (indZ*nr + ir)*np1*np2 + pind;
                    partContrib[i] += partContribNuS[i]  * rescaleFact ;
                }
            }
        if(id_unknown == id_Tcold){ // add jacobian of the rescale factor
            const real_t *const *nuS0 = nuS->GetValue(fluxGridType);
            for(len_t pind=0; pind<np1*np2; pind++){
                real_t partRescaleFact = gammaVec[pind] / Constants::mc2inEV;
                for(len_t ir=0; ir<nr; ir++)
                    partContrib[ir*np1*np2 + pind] += nuS0[ir][pind] * partRescaleFact;
            }
        }
        return partContrib;
    } 
    else if(id_unknown == unknowns->GetUnknownID(OptionConstants::UQTY_F_HOT)){
        if(!( (fluxGridType==FVM::FLUXGRIDTYPE_P1)&&(np2==1)&&(isPXiGrid) ) ){
            throw FVM::FVMException("Nonlinear contribution to collision frequencies is only implemented for hot-tails, with p-xi grid and np2=1 and evaluated on the p flux grid.");
            return nullptr;
        }
        return GetNonlinearPartialContribution(fluxGridType);
    } else {
        return nullptr;
//        throw FVM::FVMException("Invalid id_unknown: %s does not contribute to the collision frequencies",unknowns->GetUnknown(id_unknown)->GetName());
    }
}
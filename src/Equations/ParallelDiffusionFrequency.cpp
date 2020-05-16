
#include "DREAM/Equations/ParallelDiffusionFrequency.hpp"
#include "DREAM/NotImplementedException.hpp"

using namespace DREAM;



ParallelDiffusionFrequency::ParallelDiffusionFrequency(FVM::Grid *g, FVM::UnknownQuantityHandler *u, IonHandler *ih,  
                SlowingDownFrequency *nuS,
                enum OptionConstants::momentumgrid_type mgtype,  struct CollisionQuantityHandler::collqtyhand_settings *cqset)
                : CollisionQuantity(g,u,ih,mgtype,cqset){
    this->nuS = nuS;
}


void ParallelDiffusionFrequency::GetPartialContribution(len_t id_unknown, len_t ir, len_t i, len_t j, real_t *&partQty){
    nuS->GetPartialContribution(id_unknown,ir,i,j,partQty);
    rescaleFrequency(id_unknown,ir,mg->GetP(i,j),partQty);
}
void ParallelDiffusionFrequency::GetPartialContribution_fr(len_t id_unknown, len_t ir, len_t i, len_t j, real_t *&partQty){
    nuS->GetPartialContribution_fr(id_unknown,ir,i,j,partQty);
    rescaleFrequency(id_unknown,ir,mg->GetP(i,j),partQty);
}
void ParallelDiffusionFrequency::GetPartialContribution_f1(len_t id_unknown, len_t ir, len_t i, len_t j, real_t *&partQty){
    nuS->GetPartialContribution_f1(id_unknown,ir,i,j,partQty);
    rescaleFrequency(id_unknown,ir,mg->GetP_f1(i,j),partQty);
}
void ParallelDiffusionFrequency::GetPartialContribution_f2(len_t id_unknown, len_t ir, len_t i, len_t j, real_t *&partQty){
    nuS->GetPartialContribution_f2(id_unknown,ir,i,j,partQty);
    rescaleFrequency(id_unknown,ir,mg->GetP_f2(i,j),partQty);
}



void ParallelDiffusionFrequency::rescaleFrequency(len_t id_unknown, len_t ir, real_t p,real_t *&partQty){
    len_t n;
    if(id_unknown == id_ncold)
        n = 1;
    else if (id_unknown == id_ni)
        n = nzs;
    else if (id_unknown == id_fhot)
        n = nr*np1;

    for(len_t it = 0; it<n; it++){
        partQty[it] *= rescaleFactor(ir,p);
    }

}

real_t ParallelDiffusionFrequency::rescaleFactor(len_t ir, real_t p){
    real_t *Tcold = unknowns->GetUnknownData(id_Tcold);
    return (Tcold[ir]/Constants::mc2inEV)*sqrt(1+p*p);
}

real_t ParallelDiffusionFrequency::evaluateAtP(len_t ir, real_t p){
    return rescaleFactor(ir,p)*nuS->evaluateAtP(ir,p);
}



/**
 * Rebuilds partial contributions that only depend on the grid.
 */
void ParallelDiffusionFrequency::RebuildConstantTerms(){
    if(isNonlinear)
        calculateIsotropicNonlinearOperatorMatrix();
}



// Calculates Rosenbluth potential matrices defined such that when they are muliplied
// by the f_hot distribution vector, yields the three collision frequencies.
// XXX assuming for now same grid at all radii, and hot tail grid (PXi, nxi=1). 
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
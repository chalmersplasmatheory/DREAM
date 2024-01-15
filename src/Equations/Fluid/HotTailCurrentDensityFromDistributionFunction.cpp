/**
 * Implementation of an operator which evaluates the parallel current density
 * moment j_||/(B/Bmin) normalized to the local magnetic field strength
 * carried by the distribution function when using the hot-tail approximation
 * (nxi=1). The method is documented in doc/notes/theory; for low speeds it
 * uses the Lorentz approximation where j_|| ~ -E df/dp, and for high speeds
 * takes the limit j_|| ~ sign(E)*v*f.
 */

#include "DREAM/Equations/Fluid/HotTailCurrentDensityFromDistributionFunction.hpp"
#include <limits>

using namespace DREAM;


/**
 * Constructor.
 */
HotTailCurrentDensityFromDistributionFunction::HotTailCurrentDensityFromDistributionFunction(
    FVM::Grid *fluidGrid, FVM::Grid *hottailGrid, FVM::UnknownQuantityHandler *u, PitchScatterFrequency *nuD,
    enum OptionConstants::collqty_collfreq_mode collfreq_mode, bool withFullJacobian
) : EquationTerm(fluidGrid), fluidGrid(fluidGrid), hottailGrid(hottailGrid), unknowns(u), nuD(nuD) {

    SetName("HotTailCurrentDensityFromDistributionFunction");

    id_fhot  = unknowns->GetUnknownID(OptionConstants::UQTY_F_HOT);
    id_Eterm = unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD);
    id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    id_ni    = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);

    AddUnknownForJacobian(u,id_fhot);
    AddUnknownForJacobian(u,id_Eterm);
    AddUnknownForJacobian(u,id_ncold);
    AddUnknownForJacobian(u,id_Tcold);
    if(withFullJacobian)
        AddUnknownForJacobian(u,id_ni);

    
    
    isCollFreqModeFULL = (collfreq_mode == OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL);
}


/**
 * Destructor.
 */
HotTailCurrentDensityFromDistributionFunction::~HotTailCurrentDensityFromDistributionFunction() { 
    Deallocate();
}



/**
 * Rebuilds the quadrature for this integral term.
 */
void HotTailCurrentDensityFromDistributionFunction::Rebuild(const real_t,const real_t, FVM::UnknownQuantityHandler*) {
    if(!hasBeenInitialised)
        hasBeenInitialised = GridRebuilt();

    const real_t *Eterm =  unknowns->GetUnknownData(id_Eterm);
    
    for(len_t ir=0; ir<nr; ir++)
        for(len_t i=0; i<np[ir]; i++)
            nuD_vec[ir][i] = nuD->evaluateAtP(ir,hottailGrid->GetMomentumGrid(ir)->GetP1(i));

    SetJ1Weights(Eterm, nuD_vec, this->J1Weights);

    const real_t *f = unknowns->GetUnknownData(id_fhot);
    real_t fDist;
    len_t offset = 0;
    for(len_t ir=0; ir<nr; ir++){
        j1Vec[ir] = 0;
        j2Vec[ir] = 0;
        real_t sgnE = (Eterm[ir]>0) - (Eterm[ir]<0);
        for(len_t i=0; i<np[ir]; i++){
            fDist = f[offset+i];
            // If using collfreq_mode FULL, subtract an analytic Maxwellian of 
            // density ncold and temperature Tcold from the distribution
            if(isCollFreqModeFULL){
                const real_t p = hottailGrid->GetMomentumGrid(ir)->GetP1(i);
                const real_t ncold = unknowns->GetUnknownData(id_ncold)[ir];
                const real_t Tcold = unknowns->GetUnknownData(id_Tcold)[ir];
                fDist -= Constants::RelativisticMaxwellian(p, ncold, Tcold);
            }
            
            j1Vec[ir] += J1Weights[ir][i]*fDist;
            j2Vec[ir] += sgnE * J2Weights[ir][i]*fDist;
        }
        offset += np[ir];
    }
}

/**
 * Sets the matrix elements in the matrix. The quadrature is described in
 * doc/notes/theory in the section Discretization of hot-tail current.
 */
void HotTailCurrentDensityFromDistributionFunction::SetJ1Weights(const real_t *Eterm, const real_t *const*nu_D, real_t **J1Weights){
    const real_t *EffPass = fluidGrid->GetRadialGrid()->GetEffPassFrac();
    const real_t *Bavg = fluidGrid->GetRadialGrid()->GetFSA_B2();

    real_t g0 = 2*M_PI/3.0 * Constants::ec * Constants::c; // prefactor
    for(len_t ir = 0; ir<nr; ir++){
        FVM::MomentumGrid *mg = hottailGrid->GetMomentumGrid(ir);
        // Reset weights
        for(len_t i=0; i<np[ir]; i++)
            J1Weights[ir][i] = 0;

        real_t E = Constants::ec * Eterm[ir] / (Constants::me * Constants::c * sqrt(Bavg[ir]));
        real_t gConst = g0 * E * EffPass[ir];
        // set weights
        for(len_t i=0; i<np[ir]-1; i++){ // assume df/dp=0 at p=0 and at p=pMax
            real_t p = mg->GetP1(i);
            real_t p2 = p*p;
            real_t vp2 = p*p2/sqrt(1+p2);
            real_t Delta_p = mg->GetDp1(i);
            real_t g_i = Delta_p * gConst * vp2 / nu_D[ir][i];
            real_t delta_p;
            // sets j1 = int -g*df/dp*Delta_p ~ sum_i J1Weights_i * f_i.
            if(i==0){ // assume df/dp=0 at p=0
                delta_p = mg->GetP1(0) + mg->GetP1(1); // symmetry about origin: p(-1) = -p(0)
                J1Weights[ir][1] += -g_i / delta_p;
                J1Weights[ir][0] += g_i / delta_p; // f[-1] = f[0] by assumed df/dp=0 
            } else {
                delta_p = mg->GetP1(i+1) - mg->GetP1(i-1);
                J1Weights[ir][i+1] += -g_i / delta_p;
                J1Weights[ir][i-1] += g_i / delta_p;
            }
        }
        
    }  
}



/**
 * Method that is called whenever the grid is rebuilt. 
 * Allocates memory and stores plasma-independent quantities. 
 */
bool HotTailCurrentDensityFromDistributionFunction::GridRebuilt() {
    Deallocate();

    J1Weights = new real_t*[nr];
    J2Weights = new real_t*[nr];
    diffWeights = new real_t*[nr];
    dNuDmat = new real_t*[nr];
    nuD_vec = new real_t*[nr];
    dEterm  = new real_t[nr];
    j1Vec   = new real_t[nr];
    j2Vec   = new real_t[nr];
    np      = new len_t[nr];

    for(len_t ir=0; ir<nr; ir++){
        FVM::MomentumGrid *mg = hottailGrid->GetMomentumGrid(ir);
        // storing hottail grid np since EquationTerm::n1 corresponds 
        // to the fluid grid and = 1
        np[ir] = mg->GetNp1(); 
    
        J1Weights[ir] = new real_t[np[ir]];
        J2Weights[ir] = new real_t[np[ir]];
        diffWeights[ir] = new real_t[np[ir]];
        dNuDmat[ir]     = new real_t[np[ir]];
        nuD_vec[ir]     = new real_t[np[ir]];

        dEterm[ir] = 1;
        
        // initialise J2Weights: current density quadrature in theta<<1 limit
        real_t h0 = 4*M_PI * Constants::ec * Constants::c;
        // set weights
        for(len_t i=0; i<np[ir]; i++){
            real_t p = mg->GetP1(i);
            real_t p2 = p*p;
            real_t vp2 = p*p2/sqrt(1+p2);
            J2Weights[ir][i] = h0 * vp2 * mg->GetDp1(i);
        }
    }

    return true;
}


/**
 * Deallocator
 */
void HotTailCurrentDensityFromDistributionFunction::Deallocate() {
    if(J1Weights == nullptr)
        return;

    for(len_t ir=0; ir<nr; ir++){
        delete [] J1Weights[ir];
        delete [] J2Weights[ir];
        delete [] diffWeights[ir];
        delete [] dNuDmat[ir];
    }
    delete [] J1Weights;
    delete [] J2Weights;
    delete [] diffWeights;
    delete [] dNuDmat;
    delete [] dEterm;
    delete [] j1Vec;
    delete [] j2Vec; 
    delete [] np;
}


/**
 * Set the jacobian elements for this term. 
 *
 * derivId: Unknown ID of derivative with respect to which differentiation
 *          should be done.
 * unknId:  ID of the unknown to differentiate.
 * jac:     Jacobian matrix to set elements of.
 * x:       Value of the unknown quantity.
 */
bool HotTailCurrentDensityFromDistributionFunction::SetJacobianBlock(
    const len_t /*unknId*/, const len_t derivId, FVM::Matrix *jac, const real_t* f
) {
    // return unless derivId corresponds to a quantity that J1Weights depends on
    len_t nMultiples;
    if( !HasJacobianContribution(derivId, &nMultiples) )
        return false;

    // treat f_hot block separately 
    const real_t *E = unknowns->GetUnknownData(id_Eterm); 
    if(derivId == id_fhot){
        len_t offset = 0;
        for(len_t ir=0; ir<nr; ir++){
            real_t j1 = j1Vec[ir];
            real_t j2 = j2Vec[ir];
            real_t r = hypot(j1,j2);
            real_t sgn_E = (E[ir]>0) - (E[ir]<0);

            // Keep r^3 from underflowing and causing a
            // floating-point error
            if (r > cbrt(std::numeric_limits<real_t>::min())) {
                for(len_t i=0; i<np[ir]; i++)
                    jac->SetElement(ir, offset + i, 
                        (j1*J2Weights[ir][i]*sgn_E + j2*J1Weights[ir][i] ) / r
                        - (j1*J1Weights[ir][i] + j2*J2Weights[ir][i]*sgn_E)*j1*j2 / (r*r*r)
                    );
            }
            
            offset += np[ir];
        }
    }

    // sum over multiples (e.g. ion species)
    for(len_t n=0; n<nMultiples; n++){    
        // set Jacobian of the quadrature weights
        if (derivId == id_Eterm)
            SetJ1Weights(dEterm, nuD_vec, diffWeights);
        
        else {
            // set (inverse) partial deflection frequency for this nMultiple
            for(len_t ir=0; ir<nr; ir++) 
                for(len_t i=0; i<np[ir]; i++){
                    const real_t p = hottailGrid->GetMomentumGrid(ir)->GetP1(i);
                    const real_t dNuD = nuD->evaluatePartialAtP(ir,p,derivId,n);
                    if(dNuD == 0)
                        dNuDmat[ir][i] = std::numeric_limits<real_t>::infinity();
                    else
                        dNuDmat[ir][i] = - nuD_vec[ir][i]*nuD_vec[ir][i]
                                        / dNuD; 
                }
            
            SetJ1Weights(E, dNuDmat, diffWeights);
        }

        // set elements of the Jacobian matrix
        len_t offset_r = 0;
        for(len_t ir=0; ir<nr; ir++){
            real_t j1 = j1Vec[ir];
            real_t j2 = j2Vec[ir];
            real_t r = hypot(j1,j2);

            if (std::min(j1, j2) > cbrt(std::numeric_limits<real_t>::min())) {
                real_t dj1 = 0;
                for(len_t i=0; i<np[ir];i++)
                    dj1 += diffWeights[ir][i]*f[ offset_r + i ];

                jac->SetElement(ir, nr*n + ir, j2*j2*j2*dj1/(r*r*r));
            }
            offset_r += np[ir];
        }
    }    

    if(isCollFreqModeFULL)
        AddJacobianBlockMaxwellian(derivId, jac);

    return true;
}

/**
 * Adds the contribution to JacobianBlock from the Maxwellian distribution that we
 * subtract from f_hot when collfreq_mode is FULL.
 */
void HotTailCurrentDensityFromDistributionFunction::AddJacobianBlockMaxwellian(const len_t derivId, FVM::Matrix *jac){
    if((derivId == id_ncold)||(derivId == id_Tcold)){
        const real_t *ncold = unknowns->GetUnknownData(id_ncold);
        const real_t *Tcold = unknowns->GetUnknownData(id_Tcold);
        const real_t *Eterm = unknowns->GetUnknownData(id_Eterm);
        for(len_t ir=0; ir<nr; ir++){
            real_t j1 = j1Vec[ir];
            real_t j2 = j2Vec[ir];
            real_t r = hypot(j1,j2);
            real_t sgn_E = (Eterm[ir]>0) - (Eterm[ir]<0);

            for(len_t i=0; i<np[ir]; i++){
                const real_t p = hottailGrid->GetMomentumGrid(ir)->GetP1(i);
                real_t dF=0;
                if(derivId==id_ncold)
                    Constants::RelativisticMaxwellian(p,ncold[ir], Tcold[ir],&dF, nullptr);
                else
                    Constants::RelativisticMaxwellian(p,ncold[ir], Tcold[ir],nullptr, &dF);

                jac->SetElement(ir, ir, 
                    -dF * ((j1*J2Weights[ir][i]*sgn_E + j2*J1Weights[ir][i] ) / r
                    - (j1*J1Weights[ir][i] + j2*J2Weights[ir][i]*sgn_E)*j1*j2 / (r*r*r) )
                );
            }
        }
    }
}


/**
 * Set the elements of the function vector 'F' in the non-linear
 * solver.
 *
 * vec: Vector to set elements of.
 * f:   Current value of the unknown quantity for which to evaluate
 *      this operator.
 */
void HotTailCurrentDensityFromDistributionFunction::SetVectorElements(real_t *vec, const real_t */*f*/) {
    for(len_t ir=0; ir<nr; ir++){
        real_t j1 = j1Vec[ir];
        real_t j2 = j2Vec[ir];

        if (std::min(j1, j2) > sqrt(std::numeric_limits<real_t>::min()))
            vec[ir] += j1*j2 / sqrt(j1*j1 + j2*j2);
    }
}

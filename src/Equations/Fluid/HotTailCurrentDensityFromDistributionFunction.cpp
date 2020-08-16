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
    FVM::Grid *fluidGrid, FVM::Grid *hottailGrid, FVM::UnknownQuantityHandler *u, PitchScatterFrequency *nuD
) : EquationTerm(fluidGrid), fluidGrid(fluidGrid), hottailGrid(hottailGrid), unknowns(u), nuD(nuD) {

    id_fhot  = unknowns->GetUnknownID(OptionConstants::UQTY_F_HOT);
    id_Eterm = unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD);
    id_ni    = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);

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
    const real_t *const*nu_D = nuD->GetValue();
    SetGWeights(Eterm, nu_D, this->gWeights);

    len_t offset = 0;
    for(len_t ir=0; ir<nr; ir++){
        for(len_t i=0; i<np[ir]; i++){
            // useLorentzLimit is true if the Lorentz limit j ~ df/dp
            // predicts a smaller current than the high-energy limit j ~ vf
            if(abs(gWeights[ir][i]) < abs(hWeights[ir][i]) )
                useLorentzLimit[ir][i] = true;
            else
                useLorentzLimit[ir][i] = false;
        }
        offset += np[ir];
    }
}

/**
 * Sets the matrix elements in the matrix. The quadrature is described in
 * doc/notes/theory in the section Discretization of hot-tail current.
 */
void HotTailCurrentDensityFromDistributionFunction::SetGWeights(const real_t *Eterm, const real_t *const*nu_D, real_t **gWeights){
    const real_t *EffPass = fluidGrid->GetRadialGrid()->GetEffPassFrac();
    const real_t *Bavg = fluidGrid->GetRadialGrid()->GetFSA_B2();
    for(len_t ir = 0; ir<nr; ir++){
        FVM::MomentumGrid *mg = hottailGrid->GetMomentumGrid(ir);
        // Reset weights
        for(len_t i=0; i<np[ir]; i++)
            gWeights[ir][i] = 0;

        real_t E = Constants::ec * Eterm[ir] / (Constants::me * Constants::c * sqrt(Bavg[ir]));
        real_t gConst = 4*M_PI/3.0 * Constants::ec * Constants::c * E * EffPass[ir];
        // set weights
        for(len_t i=0; i<np[ir]; i++){
            real_t p = mg->GetP1(i);
            real_t vp2 = p*p*p/sqrt(1+p*p);
            real_t g_i = gConst * vp2 / (nu_D[ir][i] * delta_p[ir][i]);
            
            // df/dp region p<pcut: sets -g*df/dp*Delta_p.
            if(i==0){ // assume df/dp=0 at p=0
                gWeights[ir][i+1] += -g_i * Delta_p[ir][i];
                gWeights[ir][i]   += g_i * Delta_p[ir][i];
            } else if(i==(np[ir]-1)) { // assume f=0 at pmax
                gWeights[ir][i] += g_i * Delta_p[ir][i];
            } else {
                gWeights[ir][i+1] += -g_i * Delta_p[ir][i];
                gWeights[ir][i-1] += g_i * Delta_p[ir][i];
            }
        }
        
    }  
}



/**
 * Method that is called whenever the grid is rebuilt. 
 * Allocates memory. 
 */
bool HotTailCurrentDensityFromDistributionFunction::GridRebuilt() {
    Deallocate();

    nr = fluidGrid->GetNr();
    np = new len_t[nr];
    Delta_p  = new real_t*[nr];
    delta_p  = new real_t*[nr];
    gWeights = new real_t*[nr];
    hWeights = new real_t*[nr];
    diffWeights = new real_t*[nr];
    dEterm  = new real_t[nr];
    dNuDmat = new real_t*[nr];
    useLorentzLimit = new bool*[nr];
    
    for(len_t ir=0; ir<nr; ir++){
        FVM::MomentumGrid *mg = hottailGrid->GetMomentumGrid(ir);
        np[ir] = mg->GetNp1();
        const real_t *p   = mg->GetP1();

        delta_p[ir] = new real_t[np[ir]];
        Delta_p[ir] = new real_t[np[ir]];
        gWeights[ir] = new real_t[np[ir]];
        hWeights[ir] = new real_t[np[ir]];
        diffWeights[ir] = new real_t[np[ir]];
        dNuDmat[ir] = new real_t[np[ir]];
        useLorentzLimit[ir] = new bool[np[ir]];

        dEterm[ir] = 1;
        delta_p[ir][0] = 2*p[0];
        Delta_p[ir][0] = mg->GetDp1(0);
        delta_p[ir][np[ir]-1] = mg->GetP1_f(np[ir]) - p[np[ir]-1];
        Delta_p[ir][np[ir]-1] = mg->GetDp1(np[ir]-1);
        for(len_t i=1; i<np[ir]-1; i++){
            delta_p[ir][i] = p[i+1] - p[i-1];
            Delta_p[ir][i] = mg->GetDp1(i);
        }

        // initialise hWeights: current density quadrature in theta<<1 limit
        real_t hConst = 4*M_PI * Constants::ec * Constants::c;
        // set weights
        for(len_t i=0; i<np[ir]; i++){
            real_t p = mg->GetP1(i);
            real_t vp2 = p*p*p/sqrt(1+p*p);
            real_t h_i = hConst * vp2;
            hWeights[ir][i] = h_i*Delta_p[ir][i];
        }
    }

    return true;
}


/**
 * Deallocator
 */
void HotTailCurrentDensityFromDistributionFunction::Deallocate() {
    if(np==nullptr)
        return;

    delete [] np;
    for(len_t ir=0; ir<nr; ir++){
        delete [] delta_p[ir];
        delete [] Delta_p[ir];
        delete [] gWeights[ir];
        delete [] hWeights[ir];
        delete [] diffWeights[ir];
        delete [] dNuDmat[ir];
    }
    delete [] delta_p;
    delete [] Delta_p;
    delete [] gWeights;
    delete [] hWeights;
    delete [] diffWeights;
    delete [] dEterm;
    delete [] dNuDmat;

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
void HotTailCurrentDensityFromDistributionFunction::SetJacobianBlock(
    const len_t unknId, const len_t derivId, FVM::Matrix *jac, const real_t* f
) {
    if (derivId == unknId)
        this->SetMatrixElements(jac,nullptr);

    // For now, ignore Jacobian with respect to p_cut.
    if( !((derivId==id_Eterm) || (derivId==id_ni) || (derivId==id_ncold) || (derivId==id_Tcold)) )
        return;
    
    len_t offset_n = 0;
    // sum over multiples (e.g. ion species)
    for(len_t n=0; n<unknowns->GetUnknown(derivId)->NumberOfMultiples(); n++){    

        // set Jacobian of the quadrature weights
        if (derivId == id_Eterm) {
            SetGWeights(dEterm, nuD->GetValue(), diffWeights);
        } else if ((derivId == id_ni) || (derivId == id_ncold) || (derivId == id_Tcold)){
            FVM::fluxGridType fgType = FVM::FLUXGRIDTYPE_DISTRIBUTION;
            const real_t *dNuD = nuD->GetUnknownPartialContribution(derivId,fgType);
            // set (inverse) partial deflection frequency for this nMultiple
            for(len_t ir=0; ir<nr; ir++) 
                for(len_t i=0; i<np[ir]; i++){
                    if(dNuD[offset_n+i] == 0)
                        dNuDmat[ir][i] = std::numeric_limits<real_t>::infinity();
                    else
                        dNuDmat[ir][i] = - nuD->GetValue(ir,i,0)*nuD->GetValue(ir,i,0)
                                        / dNuD[offset_n+i]; 
                }
            
            SetGWeights(unknowns->GetUnknownData(id_Eterm), dNuDmat, diffWeights);

            offset_n += hottailGrid->GetNCells();            
        }

        // set elements of the Jacobian matrix
        len_t offset_r = 0;
        for(len_t ir=0; ir<nr; ir++){
            real_t F = 0;
            for(len_t i=0; i<np[ir];i++)
                if(useLorentzLimit[ir][i]) // only Lorentz limit part that depends on plasma parameters
                    F += diffWeights[ir][i]*f[ offset_r + i ];

            jac->SetElement(ir, nr*n + ir, F);
            offset_r += np[ir];
        }
    }    
}



/**
 * Set the elements of the linear operator matrix corresponding to
 * this operator.
 *
 * mat: Linear operator matrix to set elements of.
 * rhs: Equation right-hand-side.
 */
void HotTailCurrentDensityFromDistributionFunction::SetMatrixElements(FVM::Matrix *mat, real_t*) {
    len_t offset = 0;
    real_t weight;
    const real_t *E = unknowns->GetUnknownData(id_Eterm);
    for(len_t ir=0; ir<nr; ir++){
        real_t sgnE = (E[ir]>0) - (E[ir]<0);
        for(len_t i=0; i<np[ir]; i++){
            if(useLorentzLimit[ir][i])
                weight = gWeights[ir][i];
            else 
                weight = hWeights[ir][i] * sgnE;

            mat->SetElement(ir, offset + i, weight );
        }
        offset += np[ir];
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
void HotTailCurrentDensityFromDistributionFunction::SetVectorElements(real_t *vec, const real_t *f) {
    len_t offset = 0;
    real_t weight;
    const real_t *E = unknowns->GetUnknownData(id_Eterm);
    for(len_t ir=0; ir<nr; ir++){
        real_t sgnE = (E[ir]>0) - (E[ir]<0);
        for(len_t i=0; i<np[ir]; i++){
            if(useLorentzLimit[ir][i])
                weight = gWeights[ir][i];
            else 
                weight = hWeights[ir][i] * sgnE;

            vec[ir] += weight * f[offset+i];
        }
        offset += np[ir];
    }
}


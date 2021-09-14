/**
 * Implementation of the kinetic ionization term
 *
 *   d n_i / dt = int( v*sigma_ion_i*f dp )
 *
 * where
 *
 *   sigma_ion_i = total ionization cross-section for ion i 
 *                 and incident electron of momentum p
 *             f = electron distribution function 
 *
 * Note that this equation is applied to a single _ion species_,
 * (and to all its charge states).
 */

#include "DREAM/ADAS.hpp"
#include "DREAM/Constants.hpp"
#include "DREAM/Equations/Fluid/IonKineticIonizationTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/NotImplementedException.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM;

// import kinetic ionization data
#include "../../Atomics/kineticionizationdata.cpp"


/**
 * Constructor.
 */
IonKineticIonizationTerm::IonKineticIonizationTerm(
    FVM::Grid *momentGrid, FVM::Grid *fGrid, len_t momentId, len_t fId, 
    FVM::UnknownQuantityHandler *u, IonHandler *ihdl, const len_t iIon, 
    OptionConstants::eqterm_ionization_mode im, bool isPXiGrid, const len_t id_nf,
    real_t pThreshold, FVM::MomentQuantity::pThresholdMode pMode
) : IonEquationTerm<FVM::MomentQuantity>(momentGrid, fGrid, momentId, fId, u, pThreshold, pMode, ihdl, iIon), 
    ionization_mode(im), isPXiGrid(isPXiGrid), id_nfast(id_nf) 
{
    SetName("IonKineticIonizationTerm");

    this->id_ions = u->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    this->FVM::MomentQuantity::AddUnknownForJacobian(u, id_ions);

    if(im==OptionConstants::EQTERM_IONIZATION_MODE_KINETIC_APPROX_JAC)
        // if approximate jacobian, here sets only a correction using the fast density (rather than full distribution)
        this->FVM::MomentQuantity::AddUnknownForJacobian(u, id_nfast);
    

    this->tableIndexIon = GetTableIndex(Zion);
    
    this->GridRebuilt();
    this->FVM::MomentQuantity::AllocateDiffIntegrand();
}

/**
 * Destructor.
 */
IonKineticIonizationTerm::~IonKineticIonizationTerm() {
    Deallocate();
}


/**
 * Allocate memory for storing the ionization rate coefficients.
 */
void IonKineticIonizationTerm::Allocate() {
    Deallocate();
    // XXX: assumes same momentum grid at all radii
    const len_t n1n2 = this->fGrid->GetNp1(0)*this->fGrid->GetNp2(0);

    this->IntegrandAllCS = new real_t*[Zion+1];
    for(len_t Z0=0; Z0<=Zion; Z0++)
        this->IntegrandAllCS[Z0] = new real_t[n1n2];

    tmpVec = new real_t[nr];
}


/**
 * Deallocate memory for the ionization rate coefficients.
 */
void IonKineticIonizationTerm::Deallocate() {
    if(IntegrandAllCS == nullptr)
        return;
    for(len_t Z0=0; Z0<=Zion; Z0++)
        delete [] this->IntegrandAllCS[Z0];
    delete [] this->IntegrandAllCS;

    delete [] tmpVec;
}


/**
 * Method called whenever the grid is rebuilt.
 */
bool IonKineticIonizationTerm::GridRebuilt() {
    this->IonEquationTerm<FVM::MomentQuantity>::GridRebuilt();
    Allocate();

    RebuildIntegrand();

    return true;
}


/**
 * Sets integrand and diffIntegrand (wrt n_i) to  
 * the appropriate values for charge number Z0
 */
void IonKineticIonizationTerm::SetIntegrand(const len_t Z0){
    ResetIntegrand();
    len_t offset = 0;

    for(len_t ir=0; ir<nr; ir++){
        len_t N = fGrid->GetMomentumGrid(ir)->GetNCells();
        real_t ni = ions->GetIonDensity(ir,iIon,Z0);
        real_t ni1 = (Z0>0) ? ions->GetIonDensity(ir,iIon,Z0-1) : 0;
        for(len_t pind=0; pind<N; pind++){
            integrand[offset+pind] = -ni * IntegrandAllCS[Z0][pind];
            if(Z0>0)
                integrand[offset+pind] += ni1 * IntegrandAllCS[Z0-1][pind];
        }
        offset += N;
    }
}

void IonKineticIonizationTerm::SetDiffIntegrand(len_t){
    len_t Z0 = this->Z0ForDiffIntegrand;
    len_t rOffset = this->rOffsetForDiffIntegrand;
    len_t offset = 0;
    for(len_t ir=0; ir<nr; ir++){
        len_t N = fGrid->GetMomentumGrid(ir)->GetNCells();
        for(len_t pind=0; pind<N; pind++){
            len_t diffOffset = rOffset * N;
            diffIntegrand[diffOffset+offset+pind] = -IntegrandAllCS[Z0][pind];
            if(Z0>0){
                diffOffset -= this->nIntegrand;
                diffIntegrand[diffOffset+offset+pind] = IntegrandAllCS[Z0-1][pind];
            }
        }
        offset += N;
    }
}


/**
 * Evaluates and stores the cross section for all charge states and momenta.
 */
void IonKineticIonizationTerm::RebuildIntegrand(){
    // XXX: assumes same momentum grid at all radii
    FVM::MomentumGrid *mg = this->fGrid->GetMomentumGrid(0);
    len_t np1 = mg->GetNp1();
    len_t np2 = mg->GetNp2();
    if(isPXiGrid)
        for(len_t Z0=0; Z0<Zion; Z0++){
            const real_t *params = kinetic_rate_table[tableIndexIon].params + Z0*this->nParamsForFit;
            for(len_t i=0; i<np1; i++){
                real_t in=0;
                const real_t p = mg->GetP1(i);
                const real_t v = Constants::c * p/sqrt(1+p*p);
                in = v * EvaluateIonizationCrossSection(p, params);
                for(len_t j=0; j<np2; j++)
                    IntegrandAllCS[Z0][np1*j+i] = in;
            }
        }
    else
        for(len_t Z0=0; Z0<Zion; Z0++){
            const real_t *params = kinetic_rate_table[tableIndexIon].params + Z0*this->nParamsForFit;
            for(len_t i=0; i<np1; i++)
                for(len_t j=0; j<np2; j++){
                    const real_t p = mg->GetP(i,j);
                    const real_t v = Constants::c * p/sqrt(1+p*p);
                    IntegrandAllCS[Z0][np1*j+i] = v * EvaluateIonizationCrossSection(p, params);
                }
        }
    // zero ionization rate for fully ionized ion
    len_t N = np1*np2;
    for(len_t i=0; i<N; i++)
        IntegrandAllCS[Zion][i]=0;
}

/**
 * Evaluate electron impact ionization cross section [m^2]
 * using fitted parameters.
 * */
real_t IonKineticIonizationTerm::EvaluateIonizationCrossSection(real_t p, const real_t *params){
    return EvaluateBCGSingleSubshell(p,params[0],params[1],params[2]);
}

/**
 * Evaluates the Burgess-Chidichimo-Garland cross section for a single subshell.
 * Described in 
 *      NA Garland et al., Impact of a minority relativistic electron tail 
 *      interacting with a thermal plasma containing high-atomic-number impurities,
 *      Phys. Plasmas 27, 040702 (2020)
 *  
 *  p:        Incident electron momentum normalized to mc
 *  C:        Prefactor of cross section formula
 *  I_pot_eV: Effective ionization potential of the subshell in eV
 *  betaStar: Near-threshold modification factor
 */
real_t IonKineticIonizationTerm::EvaluateBCGSingleSubshell(real_t p, real_t C, real_t I_pot_eV, real_t betaStar){
    real_t I_pot = I_pot_eV / DREAM::Constants::mc2inEV; // convert to units of mc2
    real_t p2 = p*p;
    real_t gamma = sqrt(1+p2);
    real_t E = p2/(gamma+1); // electron kinetic energy in units of mc2

    real_t U = E/I_pot; // electron kinetic energy per ionization threshold energy
    if(U<1)
        return 0;

    real_t preFactor = C*M_PI*Constants::a0*Constants::a0;

    // non-relativistic Burgess-Chidichimo formula
    real_t I_nonRel = preFactor * Constants::Ry*Constants::Ry/(I_pot_eV*I_pot_eV)
                    * pow( log(U), 1 + betaStar/U ) / U;

    real_t beta = p/gamma;
    // relativistic Bethe-type formula from Garland
    real_t I_rel = preFactor * Constants::alpha * Constants::alpha
                * Constants::Ry / I_pot_eV * ( log(p2/(2*I_pot)) - beta*beta );
    real_t S = 1.0 / (1 + exp( 1 - E*Constants::mc2inEV*1e-5 ) );

    // return matched formula approximately valid for all energies
    return (1-S)*I_nonRel + S*I_rel;
}


/**
 * Returns the index in the kinetic_ionization_rate struct 
 * corresponding to the ion with index iZ
 */
int_t IonKineticIonizationTerm::GetTableIndex(len_t Z){
    int_t ind  = -1;
    for(len_t i=0; i<IonKineticIonizationTerm::kinetic_rate_n; i++)
        if(IonKineticIonizationTerm::kinetic_rate_table[i].Z == Z)
            ind = i;
    if(ind==-1)
        throw ADASException(
            "IonKineticIonizationTerm: Element with charge '" LEN_T_PRINTF_FMT "' not in database (kineticionizationdata.cpp).",
            Z
        );
    
    return ind;
}


/**
 * Sets the jacobian block for this equation term 
 * utilizing the functionality in MomentQuantity.
 */
bool IonKineticIonizationTerm::SetCSJacobianBlock(
    const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t *f,
    const len_t iIon, const len_t Z0, const len_t rOffset
) {
    bool contributes = false;
    if(uqtyId==derivId && ionization_mode != OptionConstants::EQTERM_IONIZATION_MODE_KINETIC_APPROX_JAC) {
        // set distribution jacobian (uqtyId corresponds to f_hot or f_re)
        this->SetCSMatrixElements(jac,nullptr,iIon,Z0,rOffset);
        contributes = true;
    }

    if(!HasJacobianContribution(derivId))
        return contributes;
    
    len_t rowOffset0 = jac->GetRowOffset();
    len_t colOffset0 = jac->GetColOffset();
    jac->SetOffset(rowOffset0+rOffset,colOffset0);
    // set n_i jacobian
    if (derivId==id_ions){
        contributes = true;

        this->Z0ForDiffIntegrand = Z0;
        this->rOffsetForDiffIntegrand = rOffset;
        this->FVM::MomentQuantity::SetJacobianBlock(uqtyId, derivId, jac, f);
    }

    if(derivId == id_nfast){
        contributes = true;
        // Set approximate fast electron jacobian under the assumption that the kinetic
        // ionization equation term is directly proportional to the fast density:
        // if hot, integrate over hot region and divide by fast density (n_hot)
        // if re, integrate over entire distribution and divide by fast density (n_re)
        SetIntegrand(Z0);
        const real_t *n = unknowns->GetUnknownData(id_nfast);

        // Reset column vector
        for (len_t ir=0; ir < nr; ir++)
            tmpVec[ir] = 0;

        this->MomentQuantity::SetVectorElements(tmpVec, f);
        for(len_t ir=0; ir<nr; ir++)
            if (n[ir] != 0)
                jac->SetElement(ir, ir, tmpVec[ir] / n[ir]);
    } 
    
    jac->SetOffset(rowOffset0,colOffset0);

    return contributes;
}


/**
 * Sets the matrix elements of this equation term
 */
void IonKineticIonizationTerm::SetCSMatrixElements(
    FVM::Matrix *mat, real_t *rhs, const len_t /*iIon*/, const len_t Z0, const len_t rOffset
) {
    SetIntegrand(Z0); 
    len_t rowOffset0 = mat->GetRowOffset();
    len_t colOffset0 = mat->GetColOffset();
    mat->SetOffset(rowOffset0+rOffset,colOffset0);
    this->FVM::MomentQuantity::SetMatrixElements(mat,rhs);
    mat->SetOffset(rowOffset0,colOffset0);
}


/**
 * Sets vector elements for this ion and charge state 
 */
void IonKineticIonizationTerm::SetCSVectorElements(
    real_t *vec, const real_t *f, const len_t /*iIon*/, const len_t Z0, const len_t rOffset
) {
    SetIntegrand(Z0); 
    this->FVM::MomentQuantity::SetVectorElements(vec+rOffset, f);
}

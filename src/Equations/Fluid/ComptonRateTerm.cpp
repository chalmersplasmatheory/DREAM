#include "DREAM/Equations/Fluid/ComptonRateTerm.hpp"
#include "DREAM/NotImplementedException.hpp"
#include "DREAM/DREAMException.hpp"

/**
 * Implementation of a class which represents the Gamma_compton contribution to the n_re equation.
 * Employs the analytical growth rate calculated by RunawayFluid.
 */
using namespace DREAM;
/**
 * Constructor.
 */
ComptonRateTerm::ComptonRateTerm(
    FVM::Grid *g, FVM::UnknownQuantityHandler *uqn,
    RunawayFluid *rf, real_t scaleFactor
) : EquationTerm(g), REFluid(rf), scaleFactor(scaleFactor) {

    this->AllocateGamma();

    this->id_n_tot   = uqn->GetUnknownID(OptionConstants::UQTY_N_TOT);
}

/**
 * Destructor.
 */
ComptonRateTerm::~ComptonRateTerm() {
    DeallocateGamma();
}

/**
 * Allocate memory for the runaway rate.
 */
void ComptonRateTerm::AllocateGamma() {
    this->gamma = new real_t[this->grid->GetNr()];
}

/**
 * Free memory for the runaway rate.
 */
void ComptonRateTerm::DeallocateGamma() {
    delete [] this->gamma;
}

/**
 * Method called when the grid has been rebuilt.
 */
bool ComptonRateTerm::GridRebuilt() {
    DeallocateGamma();
    AllocateGamma();

    return true;
}

void ComptonRateTerm::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler *uqn) {
    const len_t nr = this->grid->GetNr();

    for (len_t ir = 0; ir < nr; ir++)
        this->gamma[ir] = REFluid->GetComptonRunawayRate(ir);

    this->data_n_tot   = uqn->GetUnknownData(id_n_tot);
}
 
void ComptonRateTerm::SetJacobianBlock(const len_t, const len_t derivId, FVM::Matrix *jac, const real_t*){
    const len_t nr  = this->grid->GetNr();
    len_t offset = 0;

    for(len_t ir=0;ir<nr;ir++){
        const len_t np1 = this->grid->GetMomentumGrid(ir)->GetNp1();
        const len_t np2 = this->grid->GetMomentumGrid(ir)->GetNp2();

        if(derivId==id_n_tot){
            jac->SetElement(offset+np1*(np2-1)+0,offset+np1*(np2-1)+0,this->gamma[ir]/this->data_n_tot[ir]);
        }
        offset += np1*np2;
    }
}

/**
 * Set the linear operator matrix elements corresponding
 * to this term.
 */
void ComptonRateTerm::SetMatrixElements(FVM::Matrix*, real_t *rhs) {
    this->SetVectorElements(rhs, nullptr);
}

/**
 * Set the non-linear function vector for this term.
 */
void ComptonRateTerm::SetVectorElements(real_t *vec, const real_t*) {
    len_t offset = 0;
    
    for (len_t ir = 0; ir < nr; ir++) {
        const len_t np1 = this->grid->GetMomentumGrid(ir)->GetNp1();
        const len_t np2 = this->grid->GetMomentumGrid(ir)->GetNp2();

        // Insert at p=0, xi=1
        vec[offset + np1*(np2-1) + 0] += scaleFactor*gamma[ir];

        offset += np1*np2;
    }
}


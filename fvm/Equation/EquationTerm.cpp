/**
 * Implementation of the 'EquationTerm' base class.
 */

#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"

using namespace DREAM::FVM;


/**
 * Constructor.
 */
EquationTerm::EquationTerm(Grid *rg)
    : grid(rg) {
    
    this->AllocateMemory();
}

/**
 * Destructor.
 */
EquationTerm::~EquationTerm() {
    this->DeallocateMemory();
}

/**
 * Allocate memory for cached grid sizes.
 */
void EquationTerm::AllocateMemory() {
    DeallocateMemory();
    this->nr = this->grid->GetNr();

    this->n1 = new len_t[nr];
    this->n2 = new len_t[nr];

    for (len_t i = 0; i < nr; i++) {
        this->n1[i] = this->grid->GetMomentumGrid(i)->GetNp1();
        this->n2[i] = this->grid->GetMomentumGrid(i)->GetNp2();
    }
}

/**
 * De-allocate memory for the cached grid sizes.
 */
void EquationTerm::DeallocateMemory() {
    if (n2 != nullptr)
        delete [] this->n2;

    if (n1 != nullptr)
        delete [] this->n1;
}

/**
 * Function called when any of the grids have been re-built.
 */
bool EquationTerm::GridRebuilt() {
    this->AllocateMemory();

    return true;
}

/**
 * Returns true if derivId is in the derivIdsJacobian vector.
 *
 * derivId:    ID of unknown quantity to differentiate w.r.t.
 * nMultiples: OUTPUT. If not 'nullptr', will contain the number
 *             of multiples of the quantity with respect to which
 *             the derivative was taken.
 */
bool EquationTerm::HasJacobianContribution(len_t derivId, len_t *nMultiples) {
    bool hasDerivIdContribution = false;
    for(len_t i_deriv = 0; i_deriv < derivIdsJacobian.size(); i_deriv++){
        if (derivId == derivIdsJacobian[i_deriv]){
            if(nMultiples != nullptr)
                *nMultiples = derivNMultiplesJacobian[i_deriv];
            hasDerivIdContribution = true;
        }
    }
    return hasDerivIdContribution;
}


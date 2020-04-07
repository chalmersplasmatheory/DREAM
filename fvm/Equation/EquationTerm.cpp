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
 * Allocate memory for cahced grid sizes.
 */
void EquationTerm::AllocateMemory() {
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


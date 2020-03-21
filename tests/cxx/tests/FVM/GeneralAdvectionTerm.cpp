/**
 * Implementation of a general advection term.
 * The aim of this module is to help verify that the 'AdvectionTerm'
 * in the TQS FVM library is implemented correctly.
 */

#include "FVM/Equation/AdvectionTerm.hpp"
#include "GeneralAdvectionTerm.hpp"


using namespace TQSTESTS::FVM;

/**
 * Constructor.
 */
GeneralAdvectionTerm::GeneralAdvectionTerm(TQS::FVM::RadialGrid *rg)
    : TQS::FVM::AdvectionTerm(rg) {
    
}

/**
 * Build the coefficients of this advection term.
 */
void GeneralAdvectionTerm::Rebuild(const real_t) {
    const len_t nr = this->grid->GetNr();
    len_t offset = 0;

    for (len_t ir = 0; ir < nr; ir++) {
        auto *mg = this->grid->GetMomentumGrid(ir);
        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();

        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1; i++) {
                Fr(ir, i, j) = F1(ir, i, j) = F2(ir, i, j) = offset + j*np1 + i;
            }
        }

        offset += np1*np2;
    }
}


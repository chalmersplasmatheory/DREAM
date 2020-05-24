/**
 * Implementation of a general diffusion term.
 * The aim of this module is to help verify that the 'DiffusionTerm'
 * in the DREAM FVM library is implemented correctly.
 */

#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "GeneralDiffusionTerm.hpp"


using namespace DREAMTESTS::FVM;

/**
 * Constructor.
 */
GeneralDiffusionTerm::GeneralDiffusionTerm(DREAM::FVM::Grid *g, const real_t v)
    : DREAM::FVM::DiffusionTerm(g, true), value(v) {
    
}

/**
 * Build the coefficients of this diffusion term.
 */
void GeneralDiffusionTerm::Rebuild(
    const real_t t, const real_t, DREAM::FVM::UnknownQuantityHandler*
) {
    const len_t nr = this->grid->GetNr();
    len_t offset = 0;

    for (len_t ir = 0; ir < nr; ir++) {
        auto *mg = this->grid->GetMomentumGrid(ir);
        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();

        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1; i++) {
                real_t v;
                if (this->value == 0)
                    v = offset + j*np1 + i + 1;
                else
                    v = this->value;

                if (t == 0) {
                    Drr(ir, i, j) = v;
                    D11(ir, i, j) = 0;
                    D22(ir, i, j) = 0;
                    D12(ir, i, j) = 0;
                    D21(ir, i, j) = 0;
                } else if (t == 1) {
                    Drr(ir, i, j) = 0;
                    D11(ir, i, j) = v;
                    D22(ir, i, j) = 0;
                    D12(ir, i, j) = 0;
                    D21(ir, i, j) = 0;
                } else if (t == 2) {
                    Drr(ir, i, j) = 0;
                    D11(ir, i, j) = 0;
                    D22(ir, i, j) = v;
                    D12(ir, i, j) = 0;
                    D21(ir, i, j) = 0;
                } else if (t == 3) {
                    Drr(ir, i, j) = 0;
                    D11(ir, i, j) = 0;
                    D22(ir, i, j) = 0;
                    D12(ir, i, j) = v;
                    D21(ir, i, j) = 0;
                } else if (t == 4) {
                    Drr(ir, i, j) = 0;
                    D11(ir, i, j) = 0;
                    D22(ir, i, j) = 0;
                    D12(ir, i, j) = 0;
                    D21(ir, i, j) = v;
                } else {
                    Drr(ir, i, j) = v;
                    D11(ir, i, j) = v + 0.5;
                    D22(ir, i, j) = v + 1.0;
                    D12(ir, i, j) = v + 1.5;
                    D21(ir, i, j) = v + 2.0;
                }
            }
        }

        offset += np1*np2;
    }

    // Evaluate on flux grids
    // ir = nr
    auto *mg = this->grid->GetMomentumGrid(nr-1);
    for (len_t j = 0; j < mg->GetNp2(); j++) {
        for (len_t i = 0; i < mg->GetNp1(); i++) {
            real_t v;
            if (this->value == 0)
                v = offset + j*mg->GetNp1() + i + 1;
            else
                v = this->value;

            if (t == 0 || t > 4) Drr(nr, i, j) = v;
            else Drr(nr, i, j) = 0;
        }
    }

    // j = np2
    for (len_t ir = 0; ir < nr; ir++) {
        mg = this->grid->GetMomentumGrid(ir);
        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();

        for (len_t i = 0; i < np1; i++) {
            real_t v;
            if (this->value == 0)
                v = offset + np2*np1 + i + 1;
            else
                v = this->value;
            
            if (t == 2) {
                D22(ir, i, np2) = v;
                D21(ir, i, np2) = 0;
            } else if (t == 4) {
                D22(ir, i, np2) = 0;
                D21(ir, i, np2) = v;
            } else if (t > 4) {
                D22(ir, i, np2) = v;
                D21(ir, i, np2) = v;
            } else {
                D22(ir, i, np2) = 0;
                D21(ir, i, np2) = 0;
            }
        }
    }

    // i = np1
    for (len_t ir = 0; ir < nr; ir++) {
        mg = this->grid->GetMomentumGrid(ir);
        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();

        for (len_t j = 0; j < np2; j++) {
            real_t v;
            if (this->value == 0)
                v = offset + j*np1 + np1 + 1;
            else
                v = this->value;
            
            if (t == 1) {
                D11(ir, np1, j) = v;
                D12(ir, np1, j) = 0;
            } else if (t == 3) {
                D11(ir, np1, j) = 0;
                D12(ir, np1, j) = v;
            } else if (t > 4) {
                D11(ir, np1, j) = v;
                D12(ir, np1, j) = v;
            } else {
                D11(ir, np1, j) = 0;
                D12(ir, np1, j) = 0;
            }
        }
    }
}


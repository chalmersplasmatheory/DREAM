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
    const real_t t, const real_t dt, DREAM::FVM::UnknownQuantityHandler*
) {
    const len_t nr = this->grid->GetNr();
    len_t offset = 0;

    len_t buildIndex = static_cast<len_t>(t);
    len_t gridIndex  = static_cast<len_t>(dt);

    #define SETDRR(V) if (i < np1 && j < np2) Drr(ir,i,j) = (V)
    #define SETD11(V) if (j < np2) D11(ir,i,j) = (V)
    #define SETD12(V) if (j < np2) D12(ir,i,j) = (V)
    #define SETD21(V) if (i < np1) D21(ir,i,j) = (V)
    #define SETD22(V) if (i < np1) D22(ir,i,j) = (V)

    for (len_t ir = 0; ir < nr; ir++) {
        auto *mg = this->grid->GetMomentumGrid(ir);
        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();

        for (len_t j = 0; j < np2+1; j++) {
            for (len_t i = 0; i < np1+1; i++) {
                real_t v;
                if (this->value == 0)
                    v = offset + j*np1 + i + 1;
                else
                    v = this->value;

                if (buildIndex == 0) {
                    SETDRR(v);
                    SETD11(0);
                    SETD22(0);
                    SETD12(0);
                    SETD21(0);
                } else if (buildIndex == 1) {
                    SETDRR(0);
                    SETD11(v);
                    SETD22(0);
                    SETD12(0);
                    SETD21(0);
                } else if (buildIndex == 2) {
                    SETDRR(0);
                    SETD11(0);
                    SETD22(v);
                    SETD12(0);
                    SETD21(0);
                } else if (buildIndex == 3) {
                    SETDRR(0);
                    SETD11(0);
                    SETD22(0);
                    SETD12(v);
                    SETD21(0);
                } else if (buildIndex == 4) {
                    SETDRR(0);
                    SETD11(0);
                    SETD22(0);
                    SETD12(0);
                    SETD21(v);
                } else if (buildIndex == 101) { // D11 is non-zero in specified grid point only
                    SETDRR(0);
                    SETD22(0);
                    SETD12(0);
                    SETD21(0);

                    if (i == gridIndex) { SETD11(v); }
                    else { SETD11(0); }
                } else if (buildIndex == 102) { // D22 is non-zero in specified grid point only
                    SETDRR(0);
                    SETD11(0);
                    SETD12(0);
                    SETD21(0);

                    if (i == gridIndex) { SETD22(v); }
                    else { SETD22(0); }
                } else if (buildIndex == 103) { // D12 is non-zero in specified grid point only
                    SETDRR(0);
                    SETD11(0);
                    SETD22(0);
                    SETD21(0);

                    if (i == gridIndex) { SETD12(v); }
                    else { SETD12(0); }
                } else if (buildIndex == 104) { // D21 is non-zero in specified grid point only
                    SETDRR(0);
                    SETD11(0);
                    SETD22(0);
                    SETD12(0);

                    if (i == gridIndex) { SETD21(v); }
                    else { SETD21(0); }
                } else {
                    SETDRR(v);
                    SETD11(v + 0.5);
                    SETD22(v + 1.0);
                    SETD12(v + 1.5);
                    SETD21(v + 2.0);
                }
            }
        }

        offset += np1*np2;
    }

    // Evaluate on flux grids
    // ir = nr
    /*auto *mg = this->grid->GetMomentumGrid(nr-1);
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
    }*/
}


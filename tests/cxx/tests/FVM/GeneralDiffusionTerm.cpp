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
    this->DiffusionTerm::ResetCoefficients();
    const len_t nr = this->grid->GetNr();
    len_t offset = 0;

    len_t buildIndex = static_cast<len_t>(t);
    len_t gridIndex  = static_cast<len_t>(dt);

    #define SETDRR(V) if (ir < nr && i < np1 && j < np2 && ir != 0) Drr(ir,i,j) = (V)
    #define SETD11(V) if (ir < nr && i < np1 && j < np2 && i  != 0) D11(ir,i,j) = (V)
    #define SETD12(V) if (ir < nr && i < np1 && j < np2 && i  != 0) D12(ir,i,j) = (V)
    #define SETD21(V) if (ir < nr && i < np1 && j < np2 && j  != 0) D21(ir,i,j) = (V)
    #define SETD22(V) if (ir < nr && i < np1 && j < np2 && j  != 0) D22(ir,i,j) = (V)

    for (len_t ir = 0; ir < nr+1; ir++) {
        DREAM::FVM::MomentumGrid *mg;
        if (ir < nr)
            mg = this->grid->GetMomentumGrid(ir);
        else
            mg = this->grid->GetMomentumGrid(nr-1);
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
                } else if (buildIndex == 1) {
                    SETD11(v);
                } else if (buildIndex == 2) {
                    SETD22(v);
                } else if (buildIndex == 3) {
                    SETD12(v);
                } else if (buildIndex == 4) {
                    SETD21(v);
                } else if (buildIndex == 100 && i<np1 && j<np2) { // DRR is non-zero in specified grid point only
                    if (ir == gridIndex)
                        Drr(ir,i,j) = v;
                } else if (buildIndex == 101 && ir<nr && j<np2) { // D11 is non-zero in specified grid point only
                    if (i == gridIndex)
                        D11(ir,i,j) = v; 
                } else if (buildIndex == 102 && ir<nr && i<np1) { // D22 is non-zero in specified grid point only
                    if (i == gridIndex)
                        D22(ir,i,j) = v;
                } else if (buildIndex == 103 && ir<nr && j<np2) { // D12 is non-zero in specified grid point only
                    if (i == gridIndex) 
                        D12(ir,i,j) = v;
                } else if (buildIndex == 104 && ir<nr && i<np1) { // D21 is non-zero in specified grid point only
                    if (i == gridIndex) 
                        D21(ir,i,j) = v;
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
}


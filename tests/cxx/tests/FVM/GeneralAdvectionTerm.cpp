/**
 * Implementation of a general advection term.
 * The aim of this module is to help verify that the 'AdvectionTerm'
 * in the DREAM FVM library is implemented correctly.
 */

#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "GeneralAdvectionTerm.hpp"

using namespace DREAMTESTS::FVM;

/**
 * Constructor.
 */
GeneralAdvectionTerm::GeneralAdvectionTerm(DREAM::FVM::Grid *g, const real_t v)
    : DREAM::FVM::AdvectionTerm(g, true), value(v) {}

/**
 * Build the coefficients of this advection term.
 */
void GeneralAdvectionTerm::Rebuild(
    const real_t t, const real_t dt, DREAM::FVM::UnknownQuantityHandler*
) {
    this->AdvectionTerm::ResetCoefficients();
    const len_t nr = this->grid->GetNr();
    len_t offset = 0;

    len_t buildIndex = static_cast<len_t>(t);
    len_t gridIndex  = static_cast<len_t>(dt);

    #define SETFR(V) if (i < np1 && j < np2 && ir != 0 && ir<nr) Fr(ir,i,j) = (V)
    #define SETF1(V) if (ir < nr && j < np2 && i  != 0 && i<np1) F1(ir,i,j) = (V)
    #define SETF2(V) if (ir < nr && i < np1 && j  != 0 && j<np2) F2(ir,i,j) = (V)

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
                if (this->value == 0) // set different v in all elements
                    v = sqrt(offset + j*(np1+1) + i + 1);
                else
                    v = this->value;

                if (buildIndex == 0) {
                    SETFR(v);
                } else if (buildIndex == 1) {
                    SETF1(v);
                } else if (buildIndex == 2) {
                    SETF2(v);
                } else if (buildIndex == 100 && i<np1 && j<np2){  // Fr on index 'gridIndex' only
                    if (ir == gridIndex)
                        Fr(ir,i,j) = v; 
                } else if (buildIndex == 101 && ir<nr && j<np2){  // F1 on index 'gridIndex' only
                    if (i == gridIndex)
                        F1(ir,i,j) = v; 
                } else if (buildIndex == 102 && ir<nr && i<np1){  // F2 on index 'gridIndex' only
                    if (j == gridIndex) 
                        F2(ir,i,j) = v; 
                } else {
                    SETFR(v);
                    SETF1(v + 0.5);
                    SETF2(v + 1.0);
                }
            }
        }
        offset += np1*np2;
    }

    this->deltar->SetCoefficient(fr);
    this->delta1->SetCoefficient(f1);
    this->delta2->SetCoefficient(f2);
}

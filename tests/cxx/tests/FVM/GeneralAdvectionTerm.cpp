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
    : DREAM::FVM::AdvectionTerm(g, true), value(v) {
    
}

/**
 * Build the coefficients of this advection term.
 */
void GeneralAdvectionTerm::Rebuild(
    const real_t t, const real_t dt, DREAM::FVM::UnknownQuantityHandler*
) {
    
    const len_t nr = this->grid->GetNr();
    len_t offset = 0;

    len_t buildIndex = static_cast<len_t>(t);
    len_t gridIndex  = static_cast<len_t>(dt);

    #define SETFR(V) if (i < np1 && j < np2) Fr(ir,i,j) = (V)
    #define SETF1(V) if (ir < nr && j < np2) F1(ir,i,j) = (V)
    #define SETF2(V) if (ir < nr && i < np1) F2(ir,i,j) = (V)

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
                    SETFR(v);
                    SETF1(0);
                    SETF2(0);
                } else if (buildIndex == 1) {
                    SETFR(0);
                    SETF1(v);
                    SETF2(0);
                } else if (buildIndex == 2) {
                    SETFR(0);
                    SETF1(0);
                    SETF2(v);
                } else if (buildIndex == 101) {   // F1 on grid boundary only
                    SETFR(0);
                    SETF2(0);

                    if (i == gridIndex) { SETF1(v); }
                    else { SETF1(0); }
                } else if (buildIndex == 102) {  // F2 on grid boundary only
                    SETFR(0);
                    SETF1(0);

                    if (j == gridIndex) { SETF2(v); }
                    else { SETF2(0); }
                } else {
                    SETFR(v);
                    SETF1(v + 0.5);
                    SETF2(v + 1.0);
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
            Fr(0,i,j) = Fr(nr,i,j) = 0;
//
            real_t v;
            if (this->value == 0)
                v = offset + j*mg->GetNp1() + i + 1;
            else
                v = this->value;
            if (t == 0 || t > 2) Fr(nr, i, j) = v;
            else Fr(nr, i, j) = 0;
//
        }
    }

    // j = np2
    for (len_t ir = 0; ir < nr; ir++) {
        mg = this->grid->GetMomentumGrid(ir);
        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();

        for (len_t i = 0; i < np1; i++) {
            F2(ir,i,0) = F2(ir,i,np2) = 0;
            //
            real_t v;
            if (this->value == 0)
                v = offset + np2*np1 + i + 1;
            else
                v = this->value;
            
            if (t >= 2) F2(ir, i, np2) = v;
            else F2(ir, i, np2) = 0;
            //
        }
    }

    // i = np1
    for (len_t ir = 0; ir < nr; ir++) {
        mg = this->grid->GetMomentumGrid(ir);
        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();

        for (len_t j = 0; j < np2; j++) {
            F1(ir,0,j) = F1(ir,np1,j) = 0;
            //
            real_t v;
            if (this->value == 0)
                v = offset + j*np1 + np1 + 1;
            else
                v = this->value;
            
            if (t == 1 || t > 2) F1(ir, np1, j) = v;
            else F1(ir, np1, j) = 0;
            //
        }
    }*/
    this->deltar->SetCoefficient(fr);
    this->delta1->SetCoefficient(f1);
    this->delta2->SetCoefficient(f2);
}


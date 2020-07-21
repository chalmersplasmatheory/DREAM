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
void GeneralAdvectionTerm::Rebuild(const real_t t, const real_t, DREAM::FVM::UnknownQuantityHandler*) {
    
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
                    Fr(ir, i, j) = v;
                    F1(ir, i, j) = 0;
                    F2(ir, i, j) = 0;
                } else if (t == 1) {
                    Fr(ir, i, j) = 0;
                    F1(ir, i, j) = v;
                    F2(ir, i, j) = 0;
                } else if (t == 2) {
                    Fr(ir, i, j) = 0;
                    F1(ir, i, j) = 0;
                    F2(ir, i, j) = v;
                } else {
                    Fr(ir, i, j) = v;
                    F1(ir, i, j) = v + 0.5;
                    F2(ir, i, j) = v + 1.0;
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
            Fr(0,i,j) = Fr(nr,i,j) = 0;
/*
            real_t v;
            if (this->value == 0)
                v = offset + j*mg->GetNp1() + i + 1;
            else
                v = this->value;
            if (t == 0 || t > 2) Fr(nr, i, j) = v;
            else Fr(nr, i, j) = 0;
*/
        }
    }

    // j = np2
    for (len_t ir = 0; ir < nr; ir++) {
        mg = this->grid->GetMomentumGrid(ir);
        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();

        for (len_t i = 0; i < np1; i++) {
            F2(ir,i,0) = F2(ir,i,np2) = 0;
            /*
            real_t v;
            if (this->value == 0)
                v = offset + np2*np1 + i + 1;
            else
                v = this->value;
            
            if (t >= 2) F2(ir, i, np2) = v;
            else F2(ir, i, np2) = 0;
            */
        }
    }

    // i = np1
    for (len_t ir = 0; ir < nr; ir++) {
        mg = this->grid->GetMomentumGrid(ir);
        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();

        for (len_t j = 0; j < np2; j++) {
            F1(ir,0,j) = F1(ir,np1,j) = 0;
            /*
            real_t v;
            if (this->value == 0)
                v = offset + j*np1 + np1 + 1;
            else
                v = this->value;
            
            if (t == 1 || t > 2) F1(ir, np1, j) = v;
            else F1(ir, np1, j) = 0;
            */
        }
    }
    this->deltar->SetCoefficient(fr);
    this->delta1->SetCoefficient(f1);
    this->delta2->SetCoefficient(f2);
    
}


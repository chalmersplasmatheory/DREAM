/**
 * This file implements member methods of 'EquationSystem' that
 * provide information about the 'EquationSystem' object.
 */

#include <cstdio>
#include "DREAM/EquationSystem.hpp"


using namespace DREAM;
using namespace std;

/**
 * Prints a list of the non-trivial unknowns appearing in
 * the equation system.
 */
void EquationSystem::PrintNonTrivialUnknowns() {
    printf("LIST OF NON-TRIVIAL UNKNOWNS\n");
    printf("ID   NAME          # ELEMENTS\n");

    for (auto it = nontrivial_unknowns.begin(); it != nontrivial_unknowns.end(); it++) {
        FVM::UnknownQuantity *uq = unknowns[*it];

        printf("%3zu  %-12s  " LEN_T_PRINTF_FMT "\n", *it, uq->GetName().c_str(), uq->NumberOfElements());
    }

    printf("\n");
}

/**
 * Prints a list of the trivial unknowns (i.e. those not
 * appearing in the equation system)
 */
void EquationSystem::PrintTrivialUnknowns() {
    printf("LIST OF TRIVIAL UNKNOWNS\n");
    printf("ID   NAME\n");

    for (len_t i = 0; i < unknowns.Size(); i++) {
        if (std::find(nontrivial_unknowns.begin(), nontrivial_unknowns.end(), i) != nontrivial_unknowns.end())
            continue;

        FVM::UnknownQuantity *uq = unknowns[i];
        printf(
            "%3" LEN_T_PRINTF_FMT_PART "  %-12s\n", 
            i, uq->GetName().c_str()
        );
    }

    printf("\n");
}


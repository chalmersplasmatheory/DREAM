/**
 * This file implements member methods of 'EquationSystem' that
 * provide information about the 'EquationSystem' object.
 */

#include <cstdio>
#include "DREAM/EquationSystem.hpp"
#include "DREAM/IO.hpp"
#include "DREAM/UnknownQuantityEquation.hpp"


using namespace DREAM;
using namespace std;

/**
 * Prints a list of the externally iterated unknowns appearing in
 * the equation system.
 */
void EquationSystem::PrintExternallyIteratedUnknowns() {
	if (this->external_unknowns.size() == 0)
		return;
	
	IO::PrintInfo("LIST OF EXTERNALLY ITERATED UNKNOWNS");
    IO::PrintInfo("ID   NAME          # ELEMENTS   DESCRIPTION");

	for (len_t uqnId : external_unknowns) {
		FVM::UnknownQuantity *uq = unknowns[uqnId];
		UnknownQuantityEquation *eq = unknown_equations[uqnId];

		IO::PrintInfo(
			"%3" LEN_T_PRINTF_FMT_PART "  %-12s  %10" LEN_T_PRINTF_FMT_PART "   (%s)",
			uqnId, uq->GetName().c_str(), uq->NumberOfElements(),
			eq->GetDescription().c_str()
		);
	}

	IO::PrintInfo();
}

/**
 * Prints a list of the non-trivial unknowns appearing in
 * the equation system.
 */
void EquationSystem::PrintNonTrivialUnknowns() {
    IO::PrintInfo("LIST OF NON-TRIVIAL UNKNOWNS");
    IO::PrintInfo("ID   NAME          # ELEMENTS   DESCRIPTION");

    for (auto it = nontrivial_unknowns.begin(); it != nontrivial_unknowns.end(); it++) {
        FVM::UnknownQuantity *uq = unknowns[*it];
        UnknownQuantityEquation *eq = unknown_equations[*it];

        IO::PrintInfo(
            "%3" LEN_T_PRINTF_FMT_PART "  %-12s  %10" LEN_T_PRINTF_FMT_PART "   (%s)",
            *it, uq->GetName().c_str(), uq->NumberOfElements(),
            eq->GetDescription().c_str()
        );
    }

    IO::PrintInfo();
}

/**
 * Prints a list of the trivial unknowns (i.e. those not
 * appearing in the equation system)
 */
void EquationSystem::PrintTrivialUnknowns() {
    IO::PrintInfo("LIST OF TRIVIAL UNKNOWNS");
    IO::PrintInfo("ID   NAME           DESCRIPTION");

    for (len_t i = 0; i < unknowns.Size(); i++) {
        if (std::find(nontrivial_unknowns.begin(), nontrivial_unknowns.end(), i) != nontrivial_unknowns.end())
            continue;

        FVM::UnknownQuantity *uq = unknowns[i];
        UnknownQuantityEquation *eq = unknown_equations[i];

        IO::PrintInfo(
            "%3" LEN_T_PRINTF_FMT_PART "  %-12s   (%s)", 
            i, uq->GetName().c_str(),
            eq->GetDescription().c_str()
        );
    }

    IO::PrintInfo();
}


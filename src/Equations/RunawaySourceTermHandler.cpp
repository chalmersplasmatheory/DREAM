/**
 * A collection of fluid runaway source terms.
 */

#include <algorithm>
#include "DREAM/DREAMException.hpp"
#include "DREAM/Equations/RunawaySourceTermHandler.hpp"
#include "FVM/Equation/Operator.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"


using namespace DREAM;
using namespace std;


/**
 * Constructor.
 */
RunawaySourceTermHandler::RunawaySourceTermHandler() { }

/**
 * Destructor.
 */
RunawaySourceTermHandler::~RunawaySourceTermHandler() {
    this->applyToAll([this](FVM::EquationTerm *eqterm) {
        delete eqterm;
    });
}

/**
 * Apply the given operator to all source terms that
 * have been initialized.
 */
void RunawaySourceTermHandler::applyToAll(std::function<void(FVM::EquationTerm*)> op) {
    if (this->avalanche != nullptr)
        op(this->avalanche);
    if (this->compton != nullptr)
        op(this->compton);
    if (this->dreicer != nullptr)
        op(this->dreicer);
    if (!this->tritium.empty()) {
        for (auto t : this->tritium)
            op(t);
    }
}
void RunawaySourceTermHandler::applyToAll(const std::function<void(FVM::EquationTerm*)> op) const {
    if (this->avalanche != nullptr)
        op(this->avalanche);
    if (this->compton != nullptr)
        op(this->compton);
    if (this->dreicer != nullptr)
        op(this->dreicer);
    if (!this->tritium.empty()) {
        for (auto t : this->tritium)
            op(t);
    }
}

/**
 * Add the available source terms to the given operators.
 *
 * op_nRE:  Operator operating on n_re (runaway electron density).
 * op_nTot: Operator operating on n_tot (total electron density).
 * op_ni:   Operator operating on n_i (ions).
 */
void RunawaySourceTermHandler::AddToOperators(
    FVM::Operator *op_nRE, FVM::Operator *op_nTot,
    FVM::Operator *op_ni
) {
    // n_re
    if (this->avalanche != nullptr) {
        if (op_nRE == nullptr)
            throw DREAMException(
                "RunawaySourceTermHandler: Avalanche generation enabled, but no operator "
                "for n_re provided."
            );
        else
            op_nRE->AddTerm(this->avalanche);
    }
    if (this->dreicer != nullptr) {
        if (op_nRE == nullptr)
            throw DREAMException(
                "RunawaySourceTermHandler: Fluid Dreicer generation enabled, but no operator "
                "for n_re provided."
            );
        else
            op_nRE->AddTerm(this->dreicer);
    }

    // n_tot
    if (this->compton != nullptr) {
        if (op_nTot == nullptr)
            throw DREAMException(
                "RunawaySourceTermHandler: Compton source term enabled, but no operator "
                "for n_tot provided."
            );
        else
            op_nTot->AddTerm(this->compton);
    }

    // n_i
    if (!this->tritium.empty()) {
        if (op_ni == nullptr)
            throw DREAMException(
                "RunawaySourceTermHandler: Tritium source term enabled, but no operator "
                "for n_i provided."
            );
        else
            for (auto t : this->tritium)
                op_ni->AddTerm(t);
    }
}


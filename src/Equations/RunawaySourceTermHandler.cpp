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
    this->applyToAll([](FVM::EquationTerm *eqterm) {
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
    if (this->hottail != nullptr)
        op(this->hottail);
    if (this->lcfs_loss != nullptr)
        op(this->lcfs_loss);
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
    if (this->hottail != nullptr)
        op(this->hottail);
    if (this->lcfs_loss != nullptr)
        op(this->lcfs_loss);
    if (!this->tritium.empty()) {
        for (auto t : this->tritium)
            op(t);
    }
}

/**
 * Add the available source terms to the given operators.
 *
 * op_nRE:     Operator operating on n_re (runaway electron density).
 * op_nTot:    Operator operating on n_tot (total electron density).
 * op_ni:      Operator operating on n_i (ions).
 * op_nRE_neg: Operator operating on n_re_neg (density of runaways travelling in the xi0<0 direction).
 */
void RunawaySourceTermHandler::AddToOperators(
    FVM::Operator *op_nRE, FVM::Operator *op_nTot,
    FVM::Operator *op_ni, FVM::Operator *op_nRE_neg
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

        if (op_nRE_neg != nullptr &&
            this->avalanche_neg != nullptr &&
            this->avalanche_negpos != nullptr) {
            
            op_nRE_neg->AddTerm(this->avalanche_neg);
            op_nRE_neg->AddTerm(this->avalanche_negpos);
        }
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
    if (this->hottail != nullptr) {
        if (op_nRE == nullptr)
            throw DREAMException(
                "RunawaySourceTermHandler: Fluid hottail generation enabled, but no operator "
                "for n_re provided."
            );
        else
            op_nRE->AddTerm(this->hottail);
    }
    if (this->lcfs_loss != nullptr) {
        if (op_nRE == nullptr)
            throw DREAMException(
                "RunawaySourceTermHandler: LCFS loss term enabled, but no operator "
                "for n_re provided."
            );
        else
            op_nRE->AddTerm(this->lcfs_loss);
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


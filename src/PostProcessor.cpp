/**
 * POST PROCESSING OF UNKNOWNS
 *
 * The purpose of this class is to calculate and store the values of various
 * quantities derived from the unknown quantities solved for by DREAM. This
 * includes, for example, the runaway rate dn_RE/dt
 */

#include "DREAM/PostProcessor.hpp"
#include "DREAM/Settings/OptionConstants.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
PostProcessor::PostProcessor(
    FVM::Grid *fluidGrid, FVM::UnknownQuantityHandler *uqh, real_t p0,
    FVM::MomentQuantity::pThresholdMode pMode
) : fluidGrid(fluidGrid), unknowns(uqh), pThreshold(p0), pThresholdMode(pMode) 
{
    const len_t nr = fluidGrid->GetNr();

    // Allocate arrays
    this->runawayRate = new real_t[nr];

    // Get IDs of unknown quantities
    this->id_n_re = uqh->GetUnknownID(OptionConstants::UQTY_N_RE);
}

/**
 * Destructor.
 */
PostProcessor::~PostProcessor() {
    delete [] this->runawayRate;
}

/**
 * Process the equation system and calculate the
 * various post-processing quantities.
 *
 * t: Time of the most recently completed time step.
 */
void PostProcessor::Process(const real_t t) {
    const len_t nr = this->fluidGrid->GetNr();

    // RUNAWAY RATE
    // The runaway rate is defined as the time derivative of the runaway
    // density. It is _exactly_ \Delta n_RE / \Delta t, where \Delta n_RE is
    // the difference of the unknown quantity n_RE in the current and previous
    // time steps.
    // (Note that n_RE is _NOT_ simply the density moment of f_RE, as n_RE also
    // includes particles which have left the kinetic runaway grid)
    real_t *n_RE1 = unknowns->GetUnknownData(this->id_n_re);
    real_t *n_RE0 = unknowns->GetUnknownDataPrevious(this->id_n_re);
    real_t t0     = unknowns->GetUnknownDataPreviousTime(this->id_n_re);

    if (t > t0)
        for (len_t ir = 0; ir < nr; ir++)
            this->runawayRate[ir] = (n_RE1[ir] - n_RE0[ir]) / (t-t0);
}


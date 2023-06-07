/**
 * Class for comparing solution vector to tolerances (both absolute and
 * relative).
 */

#include <string>
#include <unordered_map>
#include <vector>
#include "DREAM/Constants.hpp"
#include "DREAM/ConvergenceChecker.hpp"
#include "DREAM/DREAMException.hpp"
#include "DREAM/IO.hpp"
#include "DREAM/Constants.hpp"
#include "DREAM/Settings/OptionConstants.hpp"


using namespace DREAM;
using namespace std;


/**
 * Constructor.
 */
ConvergenceChecker::ConvergenceChecker(
    FVM::UnknownQuantityHandler *uqh, const vector<len_t> &nontrivials,
    const real_t reltol
)
    : NormEvaluator(uqh, nontrivials) {

    this->nNontrivials = nontrivials.size();
    this->x_2norm  = new real_t[nNontrivials];
    this->dx_2norm = new real_t[nNontrivials];

    // Set relative tolerance for all
    this->SetRelativeTolerance(reltol);

    // Default absolute tolerances
    this->DefineAbsoluteTolerances();
}

/**
 * Destructor.
 */
ConvergenceChecker::~ConvergenceChecker() {
    DeallocateBuffer();

    delete [] this->x_2norm;
    delete [] this->dx_2norm;
}


/**
 * Allocate memory for the 'dx' buffer.
 */
void ConvergenceChecker::AllocateBuffer(const len_t size) {
    if (this->dx_buffer != nullptr)
        DeallocateBuffer();

    this->dx_buffer = new real_t[size];
}

/**
 * Deallocate memory for the 'dx' buffer.
 */
void ConvergenceChecker::DeallocateBuffer() {
    if (this->dx_buffer != nullptr)
        delete [] this->dx_buffer;
}

/**
 * Define absolute tolerances for unknown quantities.
 */
void ConvergenceChecker::DefineAbsoluteTolerances() {
    for (len_t it : this->nontrivials) {
        const string nname = this->unknowns->GetUnknown(it)->GetName();
        this->absTols[it] = GetDefaultAbsTol(nname);
    }
}

/**
 * Return the default absolute tolerance for the specified
 * unknown quantity.
 *
 * name: Name of unknown quantity to get default absolute tolerance for.
 */
real_t ConvergenceChecker::GetDefaultAbsTol(const string &name) {
    real_t abstol_nre = 1e-10; // roughly the minimum re density that can convert the full current in ITER
    if (name == OptionConstants::UQTY_N_RE)
        return abstol_nre;
    else if ( name == OptionConstants::UQTY_J_RE)
        return Constants::ec*Constants::c*abstol_nre;
    else if ( name == OptionConstants::UQTY_F_RE )
        return abstol_nre; 
    else
        return 0;   // No absolute tolerance check
}

/**
 * Returns the error scale (assumes that 'x_2norm' has been
 * updated prior to calling this method). The error scale is
 *
 *   scale = abstol + reltol*x_2norm
 *
 * uqty: ID of unknown quantity to get error scale for.
 */
const real_t ConvergenceChecker::GetErrorScale(const len_t uqty) {
    // Convert from "unknown quantity ID" to "nontrivial quantity index"...
    auto it = std::find(this->nontrivials.begin(), this->nontrivials.end(), uqty);
    if (it == this->nontrivials.end())
        throw DREAMException(
            "The specified unknown quantity is not a non-trivial quantity: " LEN_T_PRINTF_FMT,
            uqty
        );

    len_t idx = std::distance(this->nontrivials.begin(), it);

    return
        this->absTols[idx] +
        this->relTols[idx]*x_2norm[idx];
}

/**
 * Check if the solution is converged. This method returns true if, for
 * every non-trivial unknown quantity,
 *
 *   |dx| <= abstol + reltol*|x|
 *
 * where |.| denotes the 2-norm.
 *
 * x:  Solution vector.
 * dx: Last step 
 */
bool ConvergenceChecker::IsConverged(const real_t *x, const real_t *x1, const real_t *x2, bool verbose) {
    const len_t SSIZE = this->unknowns->GetLongVectorSize(this->nontrivials);
    if (this->dx_size != SSIZE)
        AllocateBuffer(SSIZE);

    for (len_t i = 0; i < SSIZE; i++)
        this->dx_buffer[i] = x1[i]-x2[i];

    return IsConverged(x, this->dx_buffer, verbose);
}
bool ConvergenceChecker::IsConverged(const real_t *x, const real_t *dx, bool verbose) {
    this->Norm2(x, this->x_2norm);
    this->Norm2(dx, this->dx_2norm);

    // Iterate over norms and ensure that all are small
    const len_t N = this->nontrivials.size();
    bool converged = true;

    for (len_t i = 0; i < N; i++) {
		bool conv = true;

        const real_t epsr = this->relTols[this->nontrivials[i]];
        const real_t epsa = this->absTols[this->nontrivials[i]];

        // Is tolerance checking disabled for this quantity?
        if (epsr == 0 && epsa == 0)
            continue;

		//if(x_2norm[i]>0)
        conv = (dx_2norm[i] <= (epsa + epsr*x_2norm[i])); 

        // Guard against infinity...
        if (std::isinf(dx_2norm[i]) || std::isinf(x_2norm[i]))
            conv = false;

        if (verbose) {
#ifdef COLOR_TERMINAL
            if (conv)
                DREAM::IO::PrintInfo(
                    "   \x1B[32m%10s  |x| = %e, |dx| = %e\x1B[0m",
                    this->GetNonTrivialName(i).c_str(),
                    x_2norm[i], dx_2norm[i]
                );
             else
                DREAM::IO::PrintInfo(
                    "   \x1B[1;31m%10s  |x| = %e, |dx| = %e\x1B[0m",
                    this->GetNonTrivialName(i).c_str(),
                    x_2norm[i], dx_2norm[i]
                );
#else
            DREAM::IO::PrintInfo(
                "   %10s  |x| = %e, |dx| = %e",
                this->GetNonTrivialName(i).c_str(),
                x_2norm[i], dx_2norm[i]
            );
#endif
        }

        converged = converged && conv;
    }

	return converged;
}

/**
 * Set the absolute tolerance for the specified unknown.
 *
 * uqty:   ID of the unknown quantity to set absolute tolerance for.
 * abstol: Absolute tolerance to set.
 */
void ConvergenceChecker::SetAbsoluteTolerance(const len_t uqty, const real_t abstol) {
    auto nt = this->nontrivials;
    if (find(nt.begin(), nt.end(), uqty) != nt.end())
		this->absTols[uqty] = abstol;
}

/**
 * Set the relative tolerance for all unknowns. This method
 * overwrites any previously set relative tolerances.
 *
 * reltol: Relative tolerance to require for each unknown.
 */
void ConvergenceChecker::SetRelativeTolerance(const real_t reltol) {
    for (len_t it : this->nontrivials){
        // S_particle should always be ignored since it is expected to often fluctuate at around +/- epsilon levels 
        if(this->unknowns->GetUnknown(it)->GetName() == OptionConstants::UQTY_S_PARTICLE)
            this->relTols[it] = 0;
        else
            this->relTols[it] = reltol;
    }
}

/**
 * Set the relative tolerance for the specified unknown.
 *
 * uqty:   ID of the unknown quantity to set relative tolerance for.
 * reltol: Relative tolerance.
 */
void ConvergenceChecker::SetRelativeTolerance(const len_t uqty, const real_t reltol) {
    auto nt = this->nontrivials;
    if (find(nt.begin(), nt.end(), uqty) != nt.end())
		this->relTols[uqty] = reltol;
}


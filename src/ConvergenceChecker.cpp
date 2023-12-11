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
#include "FVM/UnknownQuantity.hpp"


using namespace DREAM;
using namespace std;


/**
 * Constructor.
 */
ConvergenceChecker::ConvergenceChecker(
    FVM::UnknownQuantityHandler *uqh, vector<UnknownQuantityEquation*> *eqns,
	const vector<len_t> &nontrivials, DiagonalPreconditioner *precond,
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

	this->unknown_eqns = eqns;

	// Set the preconditioner (to use for appropriately
	// scaling the residual vector)
	if (precond)
		this->precond = precond;
	else
		this->precond = new DiagonalPreconditioner(this->unknowns, nontrivials);
	
	// Initialize residual convergence map
	for (len_t id : nontrivials) {
		std::vector<bool> &c = this->residual_conv[id];
		c.push_back(true);

		std::vector<real_t> &e = this->residual_conv_maxerr[id];
		e.push_back(0);
	}
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
bool ConvergenceChecker::IsConverged(
	const real_t *x, const real_t *dx, bool verbose
) {
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
 * Check if the residual vector is close to zero (as must be the case at a
 * solution of the equation system). This check has to make assumptions about
 * what the typical scale of different residual elements are, and so is not
 * suited as a convergence condition. Nevertheless, it can be used to verify
 * that the system of equations are satisfied in the given solution and that
 * the solver has not diverged.
 *
 * iTimeStep: Time step index of current step.
 * dt:        Length of time step.
 * F:         Residual vector.
 */
bool ConvergenceChecker::IsResidualConverged(
	const len_t iTimeStep, const real_t dt, const real_t *F
) {
	bool converged = true;
    const len_t N = this->nontrivials.size();

    for (len_t i = 0, idx = 0; i < N; i++) {
		bool conv = true;
        const real_t epsr = this->relTols[this->nontrivials[i]];
        const real_t epsa = this->absTols[this->nontrivials[i]];

        // Is tolerance checking disabled for this quantity?
        if (epsr == 0 && epsa == 0) {
			SetResidualConverged(this->nontrivials[i], iTimeStep, true, 0);
            continue;
		}

		FVM::UnknownQuantity *uq = this->unknowns->GetUnknown(this->nontrivials[i]);
		UnknownQuantityEquation *eqn = this->unknown_eqns->at(this->nontrivials[i]);

		// Determine equation scale length
		real_t scale = this->precond->GetEquationScale(this->nontrivials[i]);
		if (eqn->HasTransientTerm())
			scale /= dt;

		// Iterate over elements
		real_t mx = 0;
		for (len_t j = 0; j < uq->NumberOfElements(); j++, idx++) {
			real_t sF = std::abs(F[idx]);
			conv = conv && (sF <= (epsa + epsr*scale));

			real_t dev = sF / scale;
			if (dev > mx)
				mx = dev;
		}

		SetResidualConverged(this->nontrivials[i], iTimeStep, conv, mx);

		converged = converged && conv;
	}

	return converged;
}

/**
 * Set a flag indicating whether or not the residual was converged in the
 * specified time step.
 */
void ConvergenceChecker::SetResidualConverged(
	const len_t uid, const len_t iTimeStep,
	const bool conv, const real_t maxerr
) {
	std::vector<bool> &c = this->residual_conv[uid];
	std::vector<real_t> &e = this->residual_conv_maxerr[uid];

	if (c.size()-1 == iTimeStep) {
		c[iTimeStep] = conv;
		e[iTimeStep] = maxerr;
	} else {
		c.push_back(conv);
		e.push_back(maxerr);
	}
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

/**
 * Save stored results from the convergence checker to an output SFile.
 */
void ConvergenceChecker::SaveData(SFile *sf, const string& path) {
	string name = path + "/convergence";
	sf->CreateStruct(name);

	// Convert "residual_converged" to a 2D array
	len_t nt = 0;
	len_t N = nNontrivials;
	uint32_t *rc = nullptr;
	real_t *re = nullptr;
	for (len_t i = 0; i < nNontrivials; i++) {
		vector<bool> &c = this->residual_conv[this->nontrivials[i]];
		vector<real_t> &e = this->residual_conv_maxerr[this->nontrivials[i]];

		if (rc == nullptr) {
			nt = c.size();
			rc = new uint32_t[N * nt];
			re = new real_t[N * nt];
		}
		
		for (len_t j = 0; j < nt; j++) {
			rc[i*nt + j] = c[j];
			re[i*nt + j] = e[j];
		}
	}

	sfilesize_t dims[2] = {N, nt};
	sf->WriteMultiUInt32Array(name + "/residual", rc, 2, dims);
	sf->WriteMultiArray(name + "/residualmaxerror", re, 2, dims);

	delete [] rc;
	delete [] re;
}


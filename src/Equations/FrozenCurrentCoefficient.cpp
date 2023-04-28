/**
 * Implementation of the equation for the frozen current coefficient.
 */

#include "DREAM/Equations/FrozenCurrentCoefficient.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/Matrix.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
FrozenCurrentCoefficient::FrozenCurrentCoefficient(
	FVM::Grid *grid, FVM::Grid *fluidGrid, FVM::Interpolator1D *I_p_presc,
	FVM::UnknownQuantityHandler *unknowns, const real_t D_I_min,
	const real_t D_I_max
) : EquationTerm(grid), fluidGrid(fluidGrid), I_p_presc(I_p_presc),
	id_D_I(unknowns->GetUnknownID(OptionConstants::UQTY_D_I)),
	id_I_p(unknowns->GetUnknownID(OptionConstants::UQTY_I_P)),
	id_j_tot(unknowns->GetUnknownID(OptionConstants::UQTY_J_TOT)),
	D_I_min(D_I_min), D_I_max(D_I_max) {
	
	this->D_I = new real_t[nMaxPoints];
	this->Ip = new real_t[nMaxPoints];
}


/**
 * Destructor.
 */
FrozenCurrentCoefficient::~FrozenCurrentCoefficient() {
	if (this->I_p_presc != nullptr)
		delete this->I_p_presc;
	
	delete [] this->D_I;
	delete [] this->Ip;
}

/**
 * Clear the memory of previous solutions.
 */
void FrozenCurrentCoefficient::ClearSolutions() {
	this->nPoints = 0;
}

/**
 * Push a value Ip(D_I) onto the stack for future estimation of
 * the diffusion coefficient D_I.
 */
void FrozenCurrentCoefficient::PushSolution(
	const real_t D_I, const real_t Ip
) {
	for (int i = std::min(nPoints-1, nMaxPoints-2); i >= 0; i--) {
		this->D_I[i+1] = this->D_I[i];
		this->Ip[i+1] = this->Ip[i];
	}

	this->D_I[0] = D_I;
	this->Ip[0] = Ip;

	nPoints++;
	if (nPoints > nMaxPoints)
		nPoints = nMaxPoints;
}

/**
 * Estimate solution with no prior information.
 */
real_t FrozenCurrentCoefficient::EstimateSolutionZeroth(
	FVM::UnknownQuantityHandler *unknowns, const real_t dt, const real_t Ip
) {
	// Estimate D_I based on formula
	const len_t nr = this->fluidGrid->GetNr();
	const real_t *dr = this->fluidGrid->GetRadialGrid()->GetDr();

	real_t Ip_old = unknowns->GetUnknownDataPrevious(this->id_I_p)[0];
	real_t dIpdt = (Ip - Ip_old) / dt;
	real_t VpVol = this->fluidGrid->GetVpVol(nr-1);
	real_t ja = unknowns->GetUnknownData(this->id_j_tot)[nr-1];

	// assume j(r>a) = 0
	real_t djdr = ja / dr[nr-1];

	return std::abs(2*M_PI / VpVol * dIpdt / djdr);
}

/**
 * Estimate solution using a linear fit.
 */
real_t FrozenCurrentCoefficient::EstimateSolutionLinear(
	const real_t Ipresc
) {
	return
		this->D_I[0] +
			(this->D_I[0] - this->D_I[1]) *
			(Ipresc - this->Ip[0]) / (this->Ip[0] - this->Ip[1]);
	/*return
		this->D_I[0] -
		(this->D_I[0] - this->D_I[1]) *
		Ipresc * (Irat-1) / (this->Ip[0] - Ip_prev);*/
	//return DI - (DIk - this->D_I_prev) * Ipresc * (Irat - 1) / (Ip - Ip_prev);
}

/**
 * Estimate solution using a quadratic fit.
 */
real_t FrozenCurrentCoefficient::EstimateSolutionQuadratic(
	const real_t Ipresc
) {
	// Since a second derivative evaluated using finite differences
	// is only accurate to ~sqrt(sqrt(eps))~1e-4, and since fast
	// convergence can be achieved close to the true solution, we switch
	// to the linear estimate very close to the solution
	if (std::abs(this->Ip[0]-Ipresc) < 1e-4*std::abs(Ipresc))
		return EstimateSolutionLinear(Ipresc);

	real_t Dk = this->D_I[0], Dk1 = this->D_I[1], Dk2 = this->D_I[2];

	real_t Ip  = (this->Ip[0]-this->Ip[1]) / (Dk-Dk1);
	real_t Ip2 = (this->Ip[1]-this->Ip[2]) / (Dk1-Dk2);
	real_t Ipp = (Ip-Ip2) / (Dk-Dk1);

	real_t sq = Ip*Ip + 2*Ipp*(this->Ip[0]-Ipresc);

	// No real solutions or double root
	if (sq <= 0)
		return EstimateSolutionLinear(Ipresc);
	
	real_t rt = sqrt(sq);
	real_t Dp = Dk + (Ip + rt) / Ipp;
	real_t Dm = Dk + (Ip - rt) / Ipp;

	// No positive solutions
	if (Dp < 0 && Dm < 0)
		return EstimateSolutionLinear(Ipresc);
	else if (Dm <= 0)
		return Dp;
	else if (Dp <= 0)
		return Dm;
	else
		return std::min(Dm, Dp);
}

/**
 * Rebuild this equation term.
 */
void FrozenCurrentCoefficient::Rebuild(
	const real_t t, const real_t dt, FVM::UnknownQuantityHandler *unknowns
) {
	real_t Ipresc = this->I_p_presc->Eval(t)[0];
	real_t Ip = unknowns->GetUnknownData(this->id_I_p)[0];
	real_t DI = unknowns->GetUnknownData(this->id_D_I)[0];

	// New time step?
	if (this->prevTime != t)
		ClearSolutions();
	else
		PushSolution(DI, Ip);

	real_t Irat;
	if (Ipresc != 0)
		Irat = Ip/Ipresc;
	else
		Irat = 2;

	if (DI == 0 && (Irat <= 1)) {
		DI = 0;
	} else if (this->nPoints <= 1) {
		if (DI == 0)
			DI = this->EstimateSolutionZeroth(unknowns, dt, Ip);
		else /*if (this->prevTime != t)*/
			DI *= Irat;
	} else {
		if (this->nPoints == 2)
			DI = this->EstimateSolutionLinear(Ipresc);
		else /*if (this->nPoints >= 3)*/
			DI = this->EstimateSolutionQuadratic(Ipresc);
	}

	// Bound lower value
	if (DI < this->D_I_min)
		DI = this->D_I_min;
	// Bound upper value
	if (DI > this->D_I_max)
		DI = this->D_I_max;
	
	this->DI_sol = DI;
	this->prevTime = t;
}

/**
 * Set the jacobian elements for this equation term.
 */
bool FrozenCurrentCoefficient::SetJacobianBlock(
	const len_t, const len_t, FVM::Matrix*, const real_t*
) {
	throw FVM::EquationTermException(
		"FrozenCurrentCoefficient: This term should be added to the external "
		"iterator."
	);
}

/**
 * Set the matrix elements for this equation term.
 */
void FrozenCurrentCoefficient::SetMatrixElements(
	FVM::Matrix*, real_t*
) {
	throw FVM::EquationTermException(
		"FrozenCurrentCoefficient: This term should be added to the external "
		"iterator."
	);
}

/**
 * Set the non-linear vector elements for this equation term.
 */
void FrozenCurrentCoefficient::SetVectorElements(
	real_t *F, const real_t*
) {
	F[0] = this->DI_sol;
}


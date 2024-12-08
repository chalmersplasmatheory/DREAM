/**
 * Diffusion coefficient for the frozen current mode, determined explicitly
 * from the runaway electron density equation.
 */

#include "DREAM/Equations/FrozenCurrent/FrozenCurrentNreCoefficient.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
FrozenCurrentNreCoefficient::FrozenCurrentNreCoefficient(
	FVM::Grid *scalarGrid, FVM::Grid *fluidGrid, FVM::Interpolator1D *I_p,
	FVM::UnknownQuantityHandler *unknowns, RunawayFluid *REfluid,
	RunawaySourceTermHandler *rsth
) : EquationTerm(scalarGrid), fluidGrid(fluidGrid),
	I_p_presc(I_p), REfluid(REfluid) {
	
	this->id_D_I   = unknowns->GetUnknownID(OptionConstants::UQTY_D_I);
	this->id_I_p   = unknowns->GetUnknownID(OptionConstants::UQTY_I_P);
	this->id_n_re  = unknowns->GetUnknownID(OptionConstants::UQTY_N_RE);
	this->id_n_tot = unknowns->GetUnknownID(OptionConstants::UQTY_N_TOT);
	this->id_n_i   = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);

	this->op_n_re = new FVM::Operator(fluidGrid);
	this->op_n_tot = new FVM::Operator(fluidGrid);
	this->op_n_i = new FVM::Operator(fluidGrid);

	// Collect RE source terms
	rsth->AddToOperators(this->op_n_re, this->op_n_tot, this->op_n_i);

	this->gamma_gen = new real_t[fluidGrid->GetNCells()];
	this->deriv = new real_t[fluidGrid->GetNCells()+1];
	this->d2ndr2 = new real_t[fluidGrid->GetNCells()];

	const len_t nr = fluidGrid->GetNCells();
	this->unit = new real_t*[nr+1];
	this->unit[0] = new real_t[nr+1];
	for (len_t ir = 1; ir < nr+1; ir++)
		this->unit[ir] = this->unit[ir-1]+1;
	for (len_t ir = 0; ir < nr+1; ir++)
		this->unit[0][ir] = 1;
}

/**
 * Destructor.
 */
FrozenCurrentNreCoefficient::~FrozenCurrentNreCoefficient() {
	delete this->I_p_presc;
	delete this->op_n_i;
	delete this->op_n_tot;
	delete this->op_n_re;

	delete [] this->gamma_gen;
	delete [] this->deriv;
	delete [] this->d2ndr2;

	delete [] this->unit[0];
	delete [] this->unit;
}

/**
 * Rebuild the diffusion coefficient.
 */
void FrozenCurrentNreCoefficient::Rebuild(
	const real_t t, const real_t dt, FVM::UnknownQuantityHandler *unknowns
) {
	real_t prevTime = unknowns->GetUnknownDataPreviousTime(this->id_n_re);
	if (this->coeffTime < prevTime) {
		const real_t *n_re = unknowns->GetUnknownDataPrevious(this->id_n_re);
		const real_t *n_tot = unknowns->GetUnknownDataPrevious(this->id_n_tot);
		const real_t *n_i = unknowns->GetUnknownDataPrevious(this->id_n_i);

		// Evaluate runaway source rate
		this->op_n_re->RebuildTerms(t, dt, unknowns);
		this->op_n_tot->RebuildTerms(t, dt, unknowns);
		this->op_n_i->RebuildTerms(t, dt, unknowns);

		// Reset vectors
		for (len_t ir = 0; ir < this->fluidGrid->GetNCells(); ir++) {
			this->gamma_gen[ir] = 0;
			this->d2ndr2[ir] = 0;
		}

		this->op_n_re->SetVectorElements(this->gamma_gen, n_re);
		this->op_n_tot->SetVectorElements(this->gamma_gen, n_tot);
		this->op_n_i->SetVectorElements(this->gamma_gen, n_i);
		
		// Integrate source term
		real_t S_gen = 0;

		FVM::RadialGrid *rgrid = this->fluidGrid->GetRadialGrid();
		const real_t *Vp = this->fluidGrid->GetVpVol();
		//const real_t *Vp_f = this->fluidGrid->GetVpVol_f();
		const real_t *dr = rgrid->GetDr();
		//const real_t *dr_f = rgrid->GetDr_f();
		const real_t *GR0 = rgrid->GetBTorG();
		const real_t *Bmin = rgrid->GetBmin();
		const real_t *FSA_R02OverR2 = rgrid->GetFSA_1OverR2();
		const len_t nr = rgrid->GetNr();

		for (len_t ir = 0; ir < nr; ir++) {
			// Transient term is defined with + sign, so
			// "net positive RE generation" = minus sign
			S_gen -=
				gamma_gen[ir] * Vp[ir] * dr[ir]
				* std::abs(
					GR0[ir] / Bmin[ir] * FSA_R02OverR2[ir]
				);
		}

		/*
		// Evaluate n_re derivative on flux grid
		deriv[0] = 0;
		for (len_t ir = 1; ir < nr; ir++)
			deriv[ir] = (n_re[ir] - n_re[ir-1]) / dr_f[ir-1];
		deriv[nr] = -n_re[nr-1] / dr[nr-1];

		// Evaluate second derivative on distribution grid
		for (len_t ir = 0; ir < nr; ir++)
			d2ndr2[ir] =
				(Vp_f[ir+1]*deriv[ir+1] - Vp_f[ir]*deriv[ir])
				/ (Vp[ir]*dr[ir]);
		*/
		this->op_transport->SetVectorElements(
			d2ndr2, n_re, unit,
			this->op_transport->GetDiffusionCoeff11(),
			this->op_transport->GetDiffusionCoeff12(),
			this->op_transport->GetDiffusionCoeff21(),
			this->op_transport->GetDiffusionCoeff22()
		);
		this->bc_transport->AddToVectorElements_c(d2ndr2, n_re, unit, TransportDiffusiveBC::NO_JACOBIAN);

		// Evaluate current loss rate
		real_t S_loss = 0;
		for (len_t ir = 0; ir < nr; ir++) {
			S_loss +=
				d2ndr2[ir] * Vp[ir] * dr[ir]
				* std::abs(
					GR0[ir] / Bmin[ir] * FSA_R02OverR2[ir]
				);
		}

		if (S_loss > 0) {
			real_t D0 = S_gen / S_loss;

			real_t Ip = unknowns->GetUnknownDataPrevious(this->id_I_p)[0];
			real_t dIp = Ip - this->I_p_presc->Eval(prevTime)[0];

			real_t j_timescale = 1e-2;
			real_t dD = 2*M_PI * dIp / (Constants::ec * Constants::c * j_timescale * S_loss);

			if (dD > D0)
				this->D_I = 2*D0;
			else if (dD < -D0)
				this->D_I = 0;
			else
				this->D_I = D0 + dD;
		} else
			this->D_I = 0;
		
		this->coeffTime = prevTime;
	}
}


/**
 * Set the operators which describe the transport
 * of the runaway electrons.
 */
void FrozenCurrentNreCoefficient::SetTransportOperators(
	FVM::DiffusionTerm *op, TransportDiffusiveBC *bc
) {
	this->op_transport = op;
	this->bc_transport = bc;
}


/**
 * Set the matrix elements for this equation term.
 */
void FrozenCurrentNreCoefficient::SetMatrixElements(
	FVM::Matrix *mat, real_t *rhs
) {
	mat->SetElement(0, 0, 1.0);
	rhs[0] = -this->D_I;
}


/**
 * Set the non-linear vector elements.
 */
void FrozenCurrentNreCoefficient::SetVectorElements(
	real_t *F, const real_t *D_I
) {
	F[0] = D_I[0] - this->D_I;
}


/**
 * Set the jacobian block for this unknown.
 */
bool FrozenCurrentNreCoefficient::SetJacobianBlock(
	const len_t uqtyId, const len_t derivId,
	FVM::Matrix *jac, const real_t*
) {
	if (derivId == uqtyId) {
		jac->SetElement(0, 0, 1.0);
		return true;
	} else
		return false;
}



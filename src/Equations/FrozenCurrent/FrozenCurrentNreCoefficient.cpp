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
	RunawaySourceTermHandler *rsth, const real_t t_adjust
) : EquationTerm(scalarGrid), fluidGrid(fluidGrid),
	I_p_presc(I_p), REfluid(REfluid), unknowns(unknowns),
	t_adjust(t_adjust) {
	
	const len_t nr = fluidGrid->GetNCells();

	this->id_D_I   = unknowns->GetUnknownID(OptionConstants::UQTY_D_I);
	this->id_I_p   = unknowns->GetUnknownID(OptionConstants::UQTY_I_P);
	this->id_n_re  = unknowns->GetUnknownID(OptionConstants::UQTY_N_RE);
	this->id_n_tot = unknowns->GetUnknownID(OptionConstants::UQTY_N_TOT);
	this->id_n_i   = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
	this->id_j_ohm = unknowns->GetUnknownID(OptionConstants::UQTY_J_OHM);

	this->op_n_re = new FVM::Operator(fluidGrid);
	this->op_n_tot = new FVM::Operator(fluidGrid);
	this->op_n_i = new FVM::Operator(fluidGrid);

	IonHandler *ions = REfluid->GetIonHandler();
	const len_t nMultiples = ions->GetNzs();
	this->nMatrixElements = nMultiples*nr;
	this->op_mat = new FVM::Matrix(
		nr, nMultiples*nr, nMultiples*nr
	);

	// Collect RE source terms
	rsth->AddToOperators(this->op_n_re, this->op_n_tot, this->op_n_i);

	this->gamma_gen = new real_t[nr];
	this->d2ndr2 = new real_t[nr];
	this->dx_vec = new real_t[nr];
	this->dSgen_dx = new real_t[nr*nMultiples];
	this->dSloss_dx = new real_t[nr*nMultiples];

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
	delete [] this->d2ndr2;
	delete [] this->dx_vec;
	delete [] this->dSgen_dx;
	delete [] this->dSloss_dx;

	delete [] this->unit[0];
	delete [] this->unit;
}

/**
 * Rebuild the diffusion coefficient.
 */
void FrozenCurrentNreCoefficient::Rebuild(
	const real_t t, const real_t dt, FVM::UnknownQuantityHandler *unknowns
) {
	const real_t *n_re = unknowns->GetUnknownData(this->id_n_re);
	const real_t *n_tot = unknowns->GetUnknownData(this->id_n_tot);
	const real_t *n_i = unknowns->GetUnknownData(this->id_n_i);
	const real_t *j_ohm = unknowns->GetUnknownData(this->id_j_ohm);
	const real_t *j_ohm_prev = unknowns->GetUnknownDataPrevious(this->id_j_ohm);

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
	/*real_t S_gen = 0;

	FVM::RadialGrid *rgrid = this->fluidGrid->GetRadialGrid();
	const real_t *Vp = this->fluidGrid->GetVpVol();
	const real_t *dr = rgrid->GetDr();
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
	this->S_gen = S_gen * Constants::ec * Constants::c / (2*M_PI);*/

	// Transient term is defined with + sign, so
	// "net positive RE generation" = minus sign
	this->S_gen = -EvaluateCurrentIntegral(gamma_gen);

	this->op_transport->SetVectorElements(
		d2ndr2, n_re, unit,
		this->op_transport->GetDiffusionCoeff11(),
		this->op_transport->GetDiffusionCoeff12(),
		this->op_transport->GetDiffusionCoeff21(),
		this->op_transport->GetDiffusionCoeff22()
	);
	this->bc_transport->AddToVectorElements_c(d2ndr2, n_re, unit, TransportDiffusiveBC::NO_JACOBIAN);

	// Evaluate current loss rate
	/*real_t S_loss = 0;
	for (len_t ir = 0; ir < nr; ir++) {
		S_loss +=
			d2ndr2[ir] * Vp[ir] * dr[ir]
			* std::abs(
				GR0[ir] / Bmin[ir] * FSA_R02OverR2[ir]
			);
	}
	this->S_loss = S_loss * Constants::ec * Constants::c / (2*M_PI);*/
	this->S_loss = EvaluateCurrentIntegral(d2ndr2);
	this->lastDt = dt;

	/*if (t > 0.3450000000000001)
		fprintf(stderr, "%.12e,%.12e\n", this->S_gen, this->S_loss);*/

	// Set diffusion coefficient
	if (S_loss > 0 && S_gen > 0) {
		real_t Iohm = EvaluateCurrentIntegral(j_ohm, false);
		real_t Iohm_p = EvaluateCurrentIntegral(j_ohm_prev, false);
		this->dIohm_dt = (Iohm - Iohm_p) / dt;
		this->D0 = (S_gen + dIohm_dt) / S_loss;

		real_t Ip = unknowns->GetUnknownData(this->id_I_p)[0];
		this->dIp = Ip - this->I_p_presc->Eval(t)[0];

		real_t x = this->dIp / (this->t_adjust * S_gen);
		this->dD = this->D0 * tanh(x);

		/*if (this->dD > this->D0) {
			this->dD = this->D0;
			this->D_I = 2*this->D0;
		} else if (this->dD < -this->D0) {
			this->dD = -this->D0;
			this->D_I = 0;
		} else
			this->D_I = this->D0 + this->dD;*/
		this->D_I = this->D0 + this->dD;
	} else
		this->D_I = 0;
}


/**
 * Evaluate a spatial integral with the appropriate
 * weight for converting a density to a current.
 */
real_t FrozenCurrentNreCoefficient::EvaluateCurrentIntegral(
	const real_t *intg, bool densityToCurrent
) {
	FVM::RadialGrid *rgrid = this->fluidGrid->GetRadialGrid();
	const real_t *Vp = this->fluidGrid->GetVpVol();
	const real_t *dr = rgrid->GetDr();
	const real_t *GR0 = rgrid->GetBTorG();
	const real_t *Bmin = rgrid->GetBmin();
	const real_t *FSA_R02OverR2 = rgrid->GetFSA_1OverR2();
	const len_t nr = rgrid->GetNr();

	real_t S = 0;
	for (len_t ir = 0; ir < nr; ir++) {
		S +=
			intg[ir] * Vp[ir] * dr[ir]
			* std::abs(
				GR0[ir] / Bmin[ir] * FSA_R02OverR2[ir]
			);
	}

	if (densityToCurrent)
		return S * Constants::ec * Constants::c / (2*M_PI);
	else
		return S / (2*M_PI);
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
	FVM::Matrix *jac, const real_t *x
) {
	if (derivId == uqtyId) {
		jac->SetElement(0, 0, 1.0);
		return true;
	}

	if (this->S_loss <= 0 || this->S_gen <=0 || this->D_I == 0)
		return false;
	
	const real_t *n_re = this->unknowns->GetUnknownData(this->id_n_re);
	const real_t *n_tot = this->unknowns->GetUnknownData(this->id_n_tot);
	const real_t *n_i = this->unknowns->GetUnknownData(this->id_n_i);

	// Evaluate d(S_gen)/dx
	bool contributes_gen = false;
	if (op_n_re->SetJacobianBlock(uqtyId, derivId, op_mat, n_re))
		contributes_gen = true;
	if (op_n_tot->SetJacobianBlock(uqtyId, derivId, op_mat, n_tot))
		contributes_gen = true;
	if (op_n_i->SetJacobianBlock(uqtyId, derivId, op_mat, n_i))
		contributes_gen = true;
	
	const len_t Nd = this->unknowns->GetUnknown(derivId)->NumberOfElements();
	if (contributes_gen) {
		op_mat->Assemble();

		for (len_t id = 0; id < Nd; id++) {
			op_mat->GetColumn(id, dx_vec);
			this->dSgen_dx[id] = -EvaluateCurrentIntegral(dx_vec);
		}

		this->op_mat->Zero();
	} else
		for (len_t id = 0; id < Nd; id++)
			this->dSgen_dx[id] = 0;

	// Evaluate D*d(S_loss)/dx
	bool contributes_loss = false;
	if (op_transport->SetJacobianBlock(uqtyId, derivId, op_mat, x))
		contributes_loss = true;
	if (bc_transport->AddToJacobianBlock(uqtyId, derivId, op_mat, x))
		contributes_loss = true;

	if (derivId == this->id_I_p)
		contributes_loss = true;
	else if (derivId == this->id_j_ohm)
		contributes_loss = true;
	
	if (contributes_loss) {
		op_mat->Assemble();

		const real_t D = this->unknowns->GetUnknownData(this->id_D_I)[0];
		// Re-scale loss-term so that it's just
		// "S_loss" instead of "D * S_loss".
		const real_t scale = (D <= 0 ? 1 : 1 / D);
		for (len_t id = 0; id < Nd; id++) {
			op_mat->GetColumn(id, dx_vec);
			dSloss_dx[id] = EvaluateCurrentIntegral(dx_vec) * scale;
		}

		this->op_mat->Zero();
	} else
		for (len_t id = 0; id < Nd; id++)
			dSloss_dx[id] = 0;
	
	if (!contributes_gen && !contributes_loss)
		return false;
	
	// Evaluate dD/dx
	real_t dD_D0 = this->dD / this->D0;
	for (len_t id = 0; id < Nd; id++) {
		/*real_t v =
			this->dSgen_dx[id] - this->dSloss_dx[id];

		if (derivId == this->id_I_p)
			v += 1;
		v /= S_loss;
		*/
		real_t dD0_dx = (
			dSgen_dx[id] - (D0 + dIohm_dt/S_loss) * dSloss_dx[id]
		);
		if (derivId == this->id_j_ohm) {
			FVM::RadialGrid *rgrid = this->fluidGrid->GetRadialGrid();
			const real_t *Vp = this->fluidGrid->GetVpVol();
			const real_t *dr = rgrid->GetDr();
			const real_t *GR0 = rgrid->GetBTorG();
			const real_t *Bmin = rgrid->GetBmin();
			const real_t *FSA_R02OverR2 = rgrid->GetFSA_1OverR2();

			real_t dIohm_dt_dx =
				Vp[id] * dr[id]
				* std::abs(
					GR0[id] / Bmin[id] * FSA_R02OverR2[id]
				) / this->lastDt;

			dD0_dx += dIohm_dt_dx;
		}

		dD0_dx /= S_loss;
		real_t T1 = dD0_dx * (1 + dD_D0);

		real_t pf = 1/(S_loss * this->t_adjust) * (1 - dD_D0*dD_D0);
		real_t T2 = -this->dIp / S_gen * dSgen_dx[id];
		if (derivId == this->id_I_p)
			T2 += 1;
		T2 *= pf;

		jac->SetElement(0, id, -(T1+T2));
	}

	return true;
}



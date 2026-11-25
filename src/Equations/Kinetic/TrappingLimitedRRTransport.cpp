/**
 * Isotropic Rechester-Rosenbluth diffusion limited by the lowest collisional
 * eigenmode of the trapped pitch interval. The cylindrical radial loss rate
 * j01^2 D_detrap/a^2 equals the detrapping rate kappa*nuD.
 */

#include <algorithm>
#include <cmath>
#include "DREAM/Constants.hpp"
#include "DREAM/Equations/Kinetic/TrappingLimitedRRTransport.hpp"


using namespace DREAM;

static constexpr real_t J01_SQUARED = 5.783185962946785;

TrappingLimitedRRTransport::TrappingLimitedRRTransport(
	FVM::Grid *g, enum OptionConstants::momentumgrid_type mgtype,
	FVM::Interpolator1D *dBB, CollisionQuantityHandler *cqh,
	FVM::UnknownQuantityHandler *unknowns, bool withIonJacobian
) : RechesterRosenbluthTransport(g, mgtype, dBB) {

	SetName("TrappingLimitedRRTransport");

	auto mg = g->GetMomentumGrid(0);
	if (mgtype != OptionConstants::MOMENTUMGRID_TYPE_PXI || mg->GetNp2() != 1)
		throw DREAMException("The 'TrappingLimitedRRTransport' operator is only applicable in isotropic (nxi = 1) mode.");

	this->nuD = cqh->GetNuD();

	this->id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
	this->id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
	this->id_ni    = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);

	AddUnknownForJacobian(unknowns, this->id_ncold);
	AddUnknownForJacobian(unknowns, this->id_Tcold);
	if (withIonJacobian)
		AddUnknownForJacobian(unknowns, this->id_ni);

	AllocateCache();
	ComputeKappa();
}

TrappingLimitedRRTransport::~TrappingLimitedRRTransport() {
	delete [] this->kappaDetrap;
	delete [] this->jacobianPrefactor;
	if (this->drrDetrap != nullptr) {
		delete [] this->drrDetrap[0];
		delete [] this->drrDetrap;
	}
}


void TrappingLimitedRRTransport::AllocateCache() {
	const len_t nr = this->grid->GetNr();

	delete [] this->kappaDetrap;
	this->kappaDetrap = new real_t[nr+1]();
	delete [] this->jacobianPrefactor;
	if (this->drrDetrap != nullptr) {
		delete [] this->drrDetrap[0];
		delete [] this->drrDetrap;
	}

	this->np1_cache = this->grid->GetMomentumGrid(0)->GetNp1();
	this->jacobianPrefactor = new real_t[(nr+1) * this->np1_cache];
	this->drrDetrap = new real_t*[nr+1];
	this->drrDetrap[0] = new real_t[(nr+1) * this->np1_cache];
	for (len_t ir = 1; ir < nr+1; ir++)
		this->drrDetrap[ir] = this->drrDetrap[ir-1] + this->np1_cache;
}


bool TrappingLimitedRRTransport::GridRebuilt() {
	this->FVM::DiffusionTerm::GridRebuilt();

	AllocateCache();
	ComputeKappa();
	return true;
}


/** Lowest even eigenvalue of -(1/2) d_xi[(1-xi^2) d_xi] on (-xiT,xiT).
 * Evaluate the normalized even Legendre series and bisect its first root.
 * Long-double accumulation avoids cancellation near the absorbing boundary.
 */
real_t TrappingLimitedRRTransport::DetrappingEigenvalue(real_t xiT) {
	if (!std::isfinite(xiT) || xiT < 0 || xiT >= 0.98)
		throw DREAMException("DetrappingEigenvalue: xiT must be finite and in [0, 0.98).");
	if (xiT == 0)
		return 0;
	const long double z = static_cast<long double>(xiT)*xiT;
	const long double pi = std::acos(-1.0L);
	auto boundary = [z](long double kappa) {
		long double term = 1, sum = 1;
		for (unsigned int n = 0; n < 200000; n += 2) {
			const long double nn = n;
			term *= ((nn*(nn+1)-2*kappa)/((nn+1)*(nn+2)))*z;
			sum += term;
			// Exact polynomial termination or floating-point underflow also
			// terminates safely for very narrow cones (large kappa).
			if (term == 0 || (nn+2 > 2*kappa && std::abs(term) < 1e-17L*std::abs(sum)))
				return sum;
		}
		throw DREAMException("DetrappingEigenvalue: Legendre series did not converge.");
	};
	long double hi = pi*pi/(8*z), lo = (1-z)*hi;
	if (hi == lo)
		return static_cast<real_t>(hi);
	if (boundary(lo) < 0 || boundary(hi) > 0)
		throw DREAMException("DetrappingEigenvalue: failed to bracket the lowest eigenvalue.");
	for (unsigned int iteration = 0; iteration < 256; iteration++) {
		const long double mid = (lo+hi)/2;
		if (hi-lo <= 1e-13L*std::max(1.0L,std::abs(mid)))
			return static_cast<real_t>(mid);
		if (boundary(mid) > 0)
			lo = mid;
		else
			hi = mid;
	}
	throw DREAMException("DetrappingEigenvalue: bisection did not converge.");
}


void TrappingLimitedRRTransport::ComputeKappa() {
	auto *rg = this->grid->GetRadialGrid();
	for (len_t ir = 0; ir <= this->grid->GetNr(); ir++) {
		const real_t xiT = rg->GetXi0TrappedBoundary_fr(ir);
		this->kappaDetrap[ir] = xiT <= 0 ? 0 : DetrappingEigenvalue(xiT);
	}
}


/**
 * Rebuild the transport coefficient.
 */
void TrappingLimitedRRTransport::Rebuild(
	const real_t t, const real_t, FVM::UnknownQuantityHandler*
) {
	const real_t *dB_B = this->EvaluateDeltaBOverB(t);

	// All radii use the same momentum grid, as in the ordinary RR operator.

	const len_t nr = this->grid->GetNr();
	auto mg = this->grid->GetMomentumGrid(0);
	FVM::RadialGrid *rg = this->grid->GetRadialGrid();

	real_t qR0 = rg->GetR0();
	const real_t minorRadius = rg->GetMinorRadius();
	const real_t *p1 = mg->GetP1();
	const len_t np1 = this->np1_cache;
	real_t *const* nuD_f2 = this->nuD->GetValue_f2();
	const real_t *drf = this->deltaRadialFlux;
	std::fill(this->drrDetrap[0], this->drrDetrap[0] + (nr+1)*np1, 0.0);

	// Major radius is set to 'inf' in cylindrical geometry.
	if (std::isinf(qR0))
		qR0 = 1.0;

	for (len_t ir = 0; ir < nr+1; ir++) {
		const real_t xiT = rg->GetXi0TrappedBoundary_fr(ir);
		const real_t geometricFactor = 0.5 * rg->GetFSA_B_f(ir) * (1 - xiT*xiT);
		const real_t rrAmplitude = M_PI * qR0 * dB_B[ir] * dB_B[ir];
		const real_t w1 = (ir < nr) ? drf[ir] : 0.0;
		const real_t w2 = (ir > 0) ? 1.0 - drf[ir] : 0.0;
		real_t *prefactorRow = this->jacobianPrefactor + ir*np1;
		const real_t limitPerNuD = minorRadius*minorRadius * this->kappaDetrap[ir] / J01_SQUARED;

		for (len_t i = 0; i < np1; i++) {
			const real_t p = p1[i];
			const real_t gamma = std::sqrt(1 + p*p);
			const real_t v = Constants::c * p / gamma;
			const real_t DrrIso = rrAmplitude * v * geometricFactor;

			prefactorRow[i] = 0.0;
			if (minorRadius <= 0)
				continue;

			if (xiT == 0) {
				Drr(ir, i, 0) += DrrIso;
				continue;
			}

			real_t nuD_f = 0.0;
			if (ir < nr)
				nuD_f += w1 * nuD_f2[ir][i];
			if (ir > 0)
				nuD_f += w2 * nuD_f2[ir-1][i];

			if (nuD_f < 0)
				continue;

			const real_t Ddetrap = limitPerNuD * nuD_f;
			this->drrDetrap[ir][i] = Ddetrap;
			if (DrrIso <= 0)
				continue;

			// Scale the harmonic sum by its larger coefficient. This handles
			// a zero cap and retains the smaller coefficient at extreme ratios.
			const real_t larger = std::max(DrrIso, Ddetrap);
			const real_t rrScaled = DrrIso/larger;
			const real_t sumScaled = rrScaled + Ddetrap/larger;
			Drr(ir, i, 0) += std::min(DrrIso, Ddetrap) / sumScaled;

			// dD/dnuD = (D_RR/(D_RR+D_detrap))^2 * dD_detrap/dnuD.
			// At zero frequency Newton uses the derivative from positive nuD.
			const real_t rrFraction = rrScaled / sumScaled;
			prefactorRow[i] = limitPerNuD * rrFraction * rrFraction;
		}
	}
}


/**
 * Differentiate each radial-face flux with respect to its two source cells.
 * The coefficient uses an interpolated nuD, not a collision frequency
 * evaluated at interpolated plasma parameters. Keep the two local dnuD/du
 * values distinct; the generic diffusion Jacobian would interpolate them
 * first and then apply the radial interpolation weights a second time.
 * Boundary flux derivatives are handled by TransportDiffusiveBC below.
 */
bool TrappingLimitedRRTransport::SetJacobianBlock(
	const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t *f
) {
	const bool diagonal = (uqtyId == derivId);
	if (diagonal && !this->coefficientsShared)
		this->SetMatrixElements(jac, nullptr);
	len_t nMultiples;
	if (!HasJacobianContribution(derivId, &nMultiples))
		return diagonal;

	const real_t *dNuD = this->nuD->GetUnknownPartialContribution(derivId, FVM::FLUXGRIDTYPE_P2);
	if (dNuD == nullptr)
		return diagonal;
	const len_t nr = this->grid->GetNr(), np = this->np1_cache;
	const len_t stride = this->grid->GetMomentumGrid(0)->GetNCells_f2();
	const real_t *dr = this->grid->GetRadialGrid()->GetDr();
	const real_t *drf = this->grid->GetRadialGrid()->GetDr_f();
	for (len_t n = 0; n < nMultiples; n++)
		for (len_t ir = 1; ir < nr; ir++)
			for (len_t i = 0; i < np; i++) {
				const real_t prefactor = this->jacobianPrefactor[ir*np+i];
				if (prefactor == 0)
					continue;
				const real_t flux = prefactor * (f[(ir-1)*np+i]-f[ir*np+i]) *
					this->grid->GetVp_fr(ir)[i] / drf[ir-1];
				for (len_t source : {ir-1, ir}) {
					const real_t weight = (source == ir) ? this->deltaRadialFlux[ir] : 1-this->deltaRadialFlux[ir];
					const real_t value = flux * weight * dNuD[(n*nr+source)*stride+i];
					jac->SetElement((ir-1)*np+i, n*nr+source, value/(dr[ir-1]*this->grid->GetVp(ir-1)[i]));
					jac->SetElement(ir*np+i, n*nr+source, -value/(dr[ir]*this->grid->GetVp(ir)[i]));
				}
			}
	return true;
}


/**
 * Coefficient derivatives used by the radial boundary condition. At the
 * outer face the collision frequency comes entirely from the last cell.
 */
void TrappingLimitedRRTransport::SetPartialDiffusionTerm(len_t derivId, len_t nMultiples) {
	if (this->dd11 == nullptr)
		AllocateDifferentiationCoefficients();
	ResetDifferentiationCoefficients();

	if (
		derivId != this->id_ncold &&
		derivId != this->id_Tcold &&
		derivId != this->id_ni
	)
		return;

	// In isotropic PXi, p_f2 is cell-centered in momentum, so the j=0 P2 slice matches the p1[i] loop.
	const real_t *dNuD_f2 = this->nuD->GetUnknownPartialContribution(derivId, FVM::FLUXGRIDTYPE_P2);
	if (dNuD_f2 == nullptr)
		return;

	const len_t nr = this->grid->GetNr();
	auto mg = this->grid->GetMomentumGrid(0);
	const len_t np1 = this->np1_cache;
	const len_t rowStride = mg->GetNCells_f2();
	const len_t perMultiple = nr * rowStride;
	const real_t *drf = this->deltaRadialFlux;

	for (len_t n = 0; n < nMultiples; n++) {
		const real_t *dNuD_n = dNuD_f2 + n*perMultiple;
		for (len_t ir = 0; ir < nr+1; ir++) {
			const real_t w1 = (ir < nr) ? drf[ir] : 0.0;
			const real_t w2 = (ir > 0) ? 1.0 - drf[ir] : 0.0;
			const real_t *prefactorRow = this->jacobianPrefactor + ir*np1;

			for (len_t i = 0; i < np1; i++) {
				const real_t prefactor = prefactorRow[i];
				if (prefactor == 0.0)
					continue;

				real_t dNuD_f = 0.0;
				if (ir < nr)
					dNuD_f += w1 * dNuD_n[ir*rowStride + i];
				if (ir > 0)
					dNuD_f += w2 * dNuD_n[(ir-1)*rowStride + i];

				dDrr(ir, i, 0, n) = prefactor * dNuD_f;
			}
		}
	}
}

/**
 * Adjust the Newton step based on the physicalness of the proposed solution.
 */

#include "DREAM/Solver/PhysicalStepAdjuster.hpp"


using namespace DREAM;
using namespace std;


/**
 * Constructor.
 */
PhysicalStepAdjuster::PhysicalStepAdjuster(
	std::vector<len_t>& nu, FVM::UnknownQuantityHandler *uqh,
	IonHandler *ions, const len_t ms
) : NewtonStepAdjuster(nu, uqh, ms), ionHandler(ions) {

	ids_nonNegativeQuantities.push_back(unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD));
	ids_nonNegativeQuantities.push_back(unknowns->GetUnknownID(OptionConstants::UQTY_N_TOT));
	ids_nonNegativeQuantities.push_back(unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD));
	if(unknowns->HasUnknown(OptionConstants::UQTY_W_COLD))
		ids_nonNegativeQuantities.push_back(unknowns->GetUnknownID(OptionConstants::UQTY_W_COLD));
	if(unknowns->HasUnknown(OptionConstants::UQTY_WI_ENER))
		ids_nonNegativeQuantities.push_back(unknowns->GetUnknownID(OptionConstants::UQTY_WI_ENER));
	if(unknowns->HasUnknown(OptionConstants::UQTY_NI_DENS))
		ids_nonNegativeQuantities.push_back(unknowns->GetUnknownID(OptionConstants::UQTY_NI_DENS));

	this->id_ni = uqh->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
}


/**
 * Destructor.
 */
PhysicalStepAdjuster::~PhysicalStepAdjuster() { }


/**
 * Check if the Newton step needs to be updated.
 */
bool PhysicalStepAdjuster::AdjustmentNeeded(
	const len_t, Vec&, FVM::BlockMatrix*
) {
	return false;
}


/**
 * Adjust solution (not used by the physical step adjuster).
 */
void PhysicalStepAdjuster::AdjustSolution(const len_t iteration, real_t *sol) {
	this->UpdateSolution(sol, this->maximalPhysicalStepLength);
}


/**
 * First evaluation of the next Newton solution.
 */
void PhysicalStepAdjuster::Reset(
	Vec&, FVM::BlockMatrix*
) {
	this->maximalPhysicalStepLength = 1;
}

/**
 * Set the initial solution and Newton step.
 */
void PhysicalStepAdjuster::SetX0(const len_t iteration, const real_t *x0, const real_t *dx) {
	this->NewtonStepAdjuster::SetX0(iteration, x0, dx);

	this->maximalPhysicalStepLength = MaximalPhysicalStepLength(
		x0, dx, iteration
	);
}


/**
 * Helper function to damping factor calculation; given a quantity
 * X0 and its change dX in the iteration, returns the maximal step
 * length (damping) such that X>threshold*X0 after the iteration
 */
real_t PhysicalStepAdjuster::MaximalStepLengthAtGridPoint(
	real_t X0, real_t dX, real_t threshold
){
	real_t maxStepAtI = 1.0;
	if(dX>X0*std::numeric_limits<real_t>::min()) // dX positive, with check to avoid overflow
		maxStepAtI = (1-threshold) * X0 / dX;
	return maxStepAtI;
}


/**
 * Returns a dampingFactor such that x1 = x0 - dampingFactor*dx satisfies
 * physically-motivated constraints, such as positivity of temperature.
 * If initial guess dx from Newton step satisfies all constraints, returns 1.
 */
real_t PhysicalStepAdjuster::MaximalPhysicalStepLength(
	const real_t *x0, const real_t *dx, len_t iteration
) {
	real_t maxStepLength = 1.0;
	real_t threshold = 0.1;

	// add those quantities which we expect to be non-negative
	// T_cold and n_cold will crash the simulation if negative, so they should always be added
	bool nonNegativeZeff = true;

	const len_t N = nontrivial_unknowns.size();
	const len_t N_nn = ids_nonNegativeQuantities.size();
	len_t offset = 0;
	// sum over non-trivial unknowns
	for (len_t it=0; it<N; it++) {
		const len_t id = nontrivial_unknowns[it];
		FVM::UnknownQuantity *uq = unknowns->GetUnknown(id);
		len_t NCells = uq->NumberOfElements();

		// check whether unknown it is a non-negative quantity
		bool isNonNegativeQuantity = false;
		for (len_t it_nn = 0; it_nn < N_nn; it_nn++)
			if(id==ids_nonNegativeQuantities[it_nn])
				isNonNegativeQuantity = true;

		// Quantities which physically cannot be negative, require that they cannot
		// be reduced by more than some threshold in each iteration.
		if(isNonNegativeQuantity)
			for(len_t i=0; i<NCells; i++){
				// require x1 > threshold*x0
				real_t maxStepAtI = MaximalStepLengthAtGridPoint(x0[offset+i], dx[offset+i], threshold);
				// if this is a stronger constaint than current maxlength, override
				if(maxStepAtI < maxStepLength && maxStepAtI>0){
					maxStepLength = maxStepAtI;
					this->limitingUnknown = id;
				}
			}

		if(nonNegativeZeff && id==id_ni){
			len_t nZ = ionHandler->GetNZ();
			const len_t *Zs = ionHandler->GetZs();
			len_t nr = NCells/uq->NumberOfMultiples();
			for(len_t ir=0; ir<nr; ir++){
				real_t nZ0Z0=0;
				real_t dnZ0Z0=0;
				for(len_t iz=0; iz<nZ; iz++)
					for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
						len_t ind = ionHandler->GetIndex(iz,Z0);
						nZ0Z0 += Z0*Z0*x0[offset+ind*nr+ir];
						dnZ0Z0 += Z0*Z0*dx[offset+ind*nr+ir];
					}
				real_t maxStepAtI = MaximalStepLengthAtGridPoint(nZ0Z0, dnZ0Z0, threshold);
				if(maxStepAtI < maxStepLength && maxStepAtI>0){
					maxStepLength = maxStepAtI;
					this->limitingUnknown = id;
				}
			}
		}
		offset += NCells;
	}

	// Add automatic damping for abnormally high number of iterations to force convergence
	bool automaticDampingWithItertion = false; // skip the below for now; the method did not seem to stabilize ill-posed cases
	if(automaticDampingWithItertion){
		real_t minDamping = 0.1;
		len_t itMax = 100;
		len_t itThresh = 30;
		if(iteration>itThresh)
			maxStepLength *= std::max(minDamping,
				1.0 - ((1.0-minDamping)*(iteration-itThresh))/(itMax-itThresh));
	}

	return maxStepLength;
}


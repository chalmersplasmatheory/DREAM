/**
 * Implementation of a prescribed kinetic parameter.
 */

#include "DREAM/Equations/Kinetic/PrescribedKineticParameter.hpp"
#include "FVM/Interpolator1D.hpp"

using namespace DREAM;


/**
 * Constructor.
 */
PrescribedKineticParameter::PrescribedKineticParameter(
	FVM::Grid *g
) : PrescribedParameter(g, FVM::Interpolator1D::INTERP_NEAREST) { }

PrescribedKineticParameter::PrescribedKineticParameter(
	FVM::Grid *g, struct dream_4d_data *data
) : PrescribedParameter(g, FVM::Interpolator1D::INTERP_NEAREST), data(data) { }


/**
 * Destructor
 */
PrescribedKineticParameter::~PrescribedKineticParameter() {
	if (data != nullptr)
		DeallocateData();
}
void PrescribedKineticParameter::DeallocateData() {
	delete [] this->data->x[0];
	delete [] this->data->x;
	delete [] this->data->t;
	delete [] this->data->r;
	delete [] this->data->p1;
	delete [] this->data->p2;
	delete this->data;
}

/**
 * Rebuild the data.
 */
void PrescribedKineticParameter::Rebuild(
	const real_t t, const real_t, FVM::UnknownQuantityHandler*
) {
	len_t it = 0;

	if (this->data->nt > 1) {
		while (this->data->t[it] < t && (it+1) < this->data->nt)
			it++;
		
		if (it > 0) {
			if ((it+1)<this->data->nt) {
				real_t dt1 = std::abs(this->data->t[it]-t);
				real_t dt2 = std::abs(this->data->t[it+1]-t);

				if (dt1 > dt2)
					it++;
			}
		}
	}

	// Have we moved to a new time step in the data?
	if (this->evalData != nullptr && this->prevTimestep == it)
		return;
	
	const len_t N = this->grid->GetNCells();
	if (this->evalData == nullptr)
		this->evalData = new real_t[N];
	
	DREAM::FVM::Interpolator3D intp3(
		data->nr, data->np2, data->np1,
		data->r, data->p2, data->p1, data->x[it],
		data->gridtype, data->ps_interp, false
	);

	intp3.Eval(
		this->grid, this->data->gridtype, FVM::FLUXGRIDTYPE_DISTRIBUTION,
		this->evalData
	);

	this->interpolatedData = this->evalData;
}


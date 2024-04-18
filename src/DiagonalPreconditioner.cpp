/**
 * Implementation of a diagonal preconditioner which transforms the equation
 * system
 *
 *   Ax = b
 *
 * into
 *
 *   PAQ^-1 Qx = Pb
 *
 * where 'P' and 'Q' are diagonal matrices. Here, 'P' can be used to rescale
 * equations, while 'Q' is used to normalize the values of unknowns.
 */

#include <string>
#include "DREAM/DiagonalPreconditioner.hpp"
#include "DREAM/DREAMException.hpp"
#include "DREAM/IO.hpp"
#include "DREAM/Settings/OptionConstants.hpp"


using namespace DREAM;
using namespace std;


/**
 * Constructor.
 *
 * unknowns:    Unknown quantity handler.
 * nontrivials: List of IDs of non-trivial unknowns (i.e. those showing up in
 *              the equation system matrix).
 */
DiagonalPreconditioner::DiagonalPreconditioner(
    FVM::UnknownQuantityHandler *unknowns, const std::vector<len_t> &nontrivials
) : uqh(unknowns), nontrivials(nontrivials) {
    
    const len_t N = unknowns->GetLongVectorSize(nontrivials);

    VecCreateSeq(MPI_COMM_WORLD, N, &this->iuqn);
    VecCreateSeq(MPI_COMM_WORLD, N, &this->eqn);

    this->SetDefaultScalings();
}

/**
 * Destructor.
 */
DiagonalPreconditioner::~DiagonalPreconditioner() {
	// These seem to (somehow) be deleted anyway...
    //VecDestroy(&this->eqn);
    //VecDestroy(&this->iuqn);
}


/**
 * Build the preconditioner vectors (this method
 * need only be called once per run).
 */
void DiagonalPreconditioner::Build() {
    real_t *p, *iq;
    len_t offset = 0;

    VecGetArray(this->iuqn, &iq);
    VecGetArray(this->eqn, &p);

    const real_t REL_ERROR = 1e-6;
    for (len_t id : this->nontrivials) {
        const len_t N = uqh->GetUnknown(id)->NumberOfElements();

        for (len_t i = 0; i < N; i++) {
            p[offset+i] = 1/this->eqn_scales[id];
            iq[offset+i]  = REL_ERROR * this->uqn_scales[id];
        }

        offset += N;
    }

    VecRestoreArray(this->iuqn, &iq);
    VecRestoreArray(this->eqn, &p);
}

/**
 * Set preconditioner equation scalings for the specified unknown.
 *
 * uqty:  ID of unknown quantity to set equation scaling for.
 * scale: Scale to normalize equation for unknown to.
 */
void DiagonalPreconditioner::SetEquationScale(
    const len_t uqty, const real_t scale
) {
    auto nt = this->nontrivials;
    if (find(nt.begin(), nt.end(), uqty) == nt.end())
        throw DREAMException(
            "DiagonalPreconditioner: Cannot set equation scale for unknown "
            "quantity '%s' as it is not a non-trivial quantity.",
            this->uqh->GetUnknown(uqty)->GetName().c_str()
        );

    this->eqn_scales[uqty] = scale;
}

/**
 * Set precondition unknown scalings for the specified unknown.
 *
 * uqty:  ID of unknown quantity to set equation scaling for.
 * scale: Scale to normalize equation for unknown to.
 */
void DiagonalPreconditioner::SetUnknownScale(
    const len_t uqty, const real_t scale
) {
    auto nt = this->nontrivials;
    if (find(nt.begin(), nt.end(), uqty) == nt.end())
        throw DREAMException(
            "DiagonalPreconditioner: Cannot set equation scale for unknown "
            "quantity '%s' as it is not a non-trivial quantity.",
            this->uqh->GetUnknown(uqty)->GetName().c_str()
        );
    
    this->uqn_scales[uqty] = scale;
}


// Default characteristic sizes of various types of quantities
const real_t DiagonalPreconditioner::DENSITY_SCALE = 1e20;
const real_t DiagonalPreconditioner::CURRENT_SCALE = 1e6;
const real_t DiagonalPreconditioner::ENERGY_SCALE = 1e6;
const real_t DiagonalPreconditioner::FLUX_SCALE = 1;
const real_t DiagonalPreconditioner::TEMPERATURE_SCALE = 100;
const real_t DiagonalPreconditioner::RUNAWAY_FRACTION = 1e-10;

/**
 * Set default scalings.
 */
void DiagonalPreconditioner::SetDefaultScalings() {
    for (len_t id : this->nontrivials) {
        const string &name = uqh->GetUnknown(id)->GetName();

        if (name == OptionConstants::UQTY_E_FIELD) {
            uqn_scales[id] = eqn_scales[id] = FLUX_SCALE;
        } else if (name == OptionConstants::UQTY_F_HOT) {
            uqn_scales[id] = eqn_scales[id] = DENSITY_SCALE;
        } else if (name == OptionConstants::UQTY_F_RE) {
            uqn_scales[id] = eqn_scales[id] = RUNAWAY_FRACTION*DENSITY_SCALE;
        } else if (name == OptionConstants::UQTY_ION_SPECIES) {
            uqn_scales[id] = DENSITY_SCALE;
            eqn_scales[id] = 1e6*DENSITY_SCALE;  // dn_i/dt has characteristic times of 1e-6 seconds
        } else if (name == OptionConstants::UQTY_I_WALL) {
            uqn_scales[id] = eqn_scales[id] = CURRENT_SCALE;  // 1 MA
        } else if (name == OptionConstants::UQTY_I_P) {
            uqn_scales[id] = eqn_scales[id] = CURRENT_SCALE;  // 1 MA
        } else if (name == OptionConstants::UQTY_J_HOT) {
            uqn_scales[id] = eqn_scales[id] = CURRENT_SCALE;  // 1 MA/m^2
        } else if (name == OptionConstants::UQTY_J_OHM) {
            uqn_scales[id] = eqn_scales[id] = CURRENT_SCALE;  // 1 MA/m^2
        } else if (name == OptionConstants::UQTY_J_RE) {
            uqn_scales[id] = eqn_scales[id] = CURRENT_SCALE;  // 1 MA/m^2
        } else if (name == OptionConstants::UQTY_J_TOT) {
            uqn_scales[id] = eqn_scales[id] = CURRENT_SCALE;  // 1 MA/m^2
        } else if (name == OptionConstants::UQTY_N_COLD) {
            uqn_scales[id] = eqn_scales[id] = DENSITY_SCALE;
        } else if (name == OptionConstants::UQTY_N_HOT) {
            uqn_scales[id] = eqn_scales[id] = DENSITY_SCALE;
        } else if (name == OptionConstants::UQTY_N_RE) {
            uqn_scales[id] = eqn_scales[id] = RUNAWAY_FRACTION*DENSITY_SCALE;
        } else if (name == OptionConstants::UQTY_N_RE_NEG) {
            uqn_scales[id] = eqn_scales[id] = RUNAWAY_FRACTION*DENSITY_SCALE;
        } else if (name == OptionConstants::UQTY_N_TOT) {
            uqn_scales[id] = eqn_scales[id] = DENSITY_SCALE;
        } else if (name == OptionConstants::UQTY_NI_DENS) {
            uqn_scales[id] = eqn_scales[id] = DENSITY_SCALE;
        } else if (name == OptionConstants::UQTY_POL_FLUX) {
            uqn_scales[id] = eqn_scales[id] = FLUX_SCALE;
        } else if (name == OptionConstants::UQTY_PSI_WALL) {
            uqn_scales[id] = eqn_scales[id] = FLUX_SCALE;
        } else if (name == OptionConstants::UQTY_PSI_EDGE) {
            uqn_scales[id] = eqn_scales[id] = FLUX_SCALE;
        } else if (name == OptionConstants::UQTY_S_PARTICLE) {
            uqn_scales[id] = eqn_scales[id] = 1e-10*DENSITY_SCALE;
        } else if (name == OptionConstants::UQTY_TAU_COLL) {
            uqn_scales[id] = eqn_scales[id] = 1e-3;
        } else if (name == OptionConstants::UQTY_T_COLD) {
            uqn_scales[id] = TEMPERATURE_SCALE; // a medium temperature
            eqn_scales[id] = ENERGY_SCALE;   // 1 MJ/m^3
        } else if (name == OptionConstants::UQTY_V_LOOP_WALL) {
            uqn_scales[id] = eqn_scales[id] = FLUX_SCALE;
        } else if (name == OptionConstants::UQTY_W_COLD) {
            uqn_scales[id] = eqn_scales[id] = ENERGY_SCALE; // 1 MJ/m^3
        } else if (name == OptionConstants::UQTY_WI_ENER) {
            uqn_scales[id] = eqn_scales[id] = ENERGY_SCALE;
        } else {
            DREAM::IO::PrintWarning(
                "DiagonalPreconditioner: Unrecognized unknown '%s'. Unknown "
                "will not be preconditioned.", name.c_str()
            );

            uqn_scales[id] = eqn_scales[id] = 1;
        }
    }
}

/**
 * Rescale the given matrix according to the transformation
 *
 *   B = PAQ^-1
 *
 * where 'P' is the equation rescaling matrix and 'Q' is the
 * unknown rescaling matrix.
 */
void DiagonalPreconditioner::RescaleMatrix(FVM::Matrix *mat) {
    mat->DiagonalScale(this->eqn, this->iuqn);
}

/**
 * Rescale the given RHS vector according to the transformation
 *
 *   c = Pb
 *
 * where 'P' is the equation rescaling matrix.
 */
void DiagonalPreconditioner::RescaleRHSVector(Vec b) {
    VecPointwiseMult(b, b, this->eqn);
}

/**
 * Rescale the given unknown vector according to the transformation
 *
 *   y = Q^-1x
 *
 * where 'Q' is the unknown rescaling matrix.
 */
void DiagonalPreconditioner::UnscaleUnknownVector(Vec Qx) {
    VecPointwiseMult(Qx, Qx, this->iuqn);
}


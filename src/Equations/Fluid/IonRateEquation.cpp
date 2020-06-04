/**
 * Implementation of the ion rate equation
 *
 *   d n_i^(j) / dt =
 *      [I_i^(j-1) n_cold + II_i^(j-1)] n_i^(j-1) -
 *      [I_i^(j) n_cold  II_i^(j)] n_i^(j-1) +
 *      R_i^(j+1) n_i^(j+1) n_cold - R_i^(j) n_i^(j) n_cold
 *
 * where
 *
 *   I_i^(j)  = ionization rate coefficient for charge state 'j'
 *              of ion species 'i'.
 *   II_i^(j) = fast-electron impact ionization coefficient for
 *              charge state 'j' of ion species 'i'.
 *   R_i^(j)  = radiative recombination rate for charge state 'j'
 *              of ion species 'i'.
 *
 * Note that this equation is applied to a single _ion species_,
 * (and to all its charge states).
 */

#include "DREAM/ADAS.hpp"
#include "DREAM/Equations/Fluid/IonRateEquation.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/NotImplementedException.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
IonRateEquation::IonRateEquation(
    FVM::Grid *g, IonHandler *ihdl, const len_t iIon,
    ADAS *adas, FVM::UnknownQuantityHandler *unknowns
) : IonEquationTerm<FVM::EquationTerm>(g, ihdl, iIon), adas(adas) {
    
    this->unknowns  = unknowns;
    this->id_ions   = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    this->id_n_cold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    this->id_n_hot  = unknowns->GetUnknownID(OptionConstants::UQTY_N_HOT);
    this->id_n_tot  = unknowns->GetUnknownID(OptionConstants::UQTY_N_TOT);
    this->id_T_cold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);

    AllocateRateCoefficients();
}

/**
 * Destructor.
 */
IonRateEquation::~IonRateEquation() {
    DeallocateRateCoefficients();
}


/**
 * Allocate memory for storing the ionization rate coefficients.
 */
void IonRateEquation::AllocateRateCoefficients() {
    const len_t Nr  = this->grid->GetNr();

    this->Rec = new real_t*[(Zion+1)];
    this->Ion = new real_t*[(Zion+1)];
    this->Imp = new real_t*[(Zion+1)];

    this->Rec[0] = new real_t[Nr*(Zion+1)];
    this->Ion[0] = new real_t[Nr*(Zion+1)];
    this->Imp[0] = new real_t[Nr*(Zion+1)];

    for (len_t i = 1; i <= Zion; i++) {
        this->Rec[i] = this->Rec[i-1] + Nr;
        this->Ion[i] = this->Ion[i-1] + Nr;
        this->Imp[i] = this->Imp[i-1] + Nr;
    }
}

/**
 * Deallocate memory for the ionization rate coefficients.
 */
void IonRateEquation::DeallocateRateCoefficients() {
    delete [] this->Imp[0];
    delete [] this->Ion[0];
    delete [] this->Rec[0];

    delete [] this->Imp;
    delete [] this->Ion;
    delete [] this->Rec;
}

/**
 * Method called whenever the grid is rebuilt.
 */
bool IonRateEquation::GridRebuilt() {
    this->IonEquationTerm<FVM::EquationTerm>::GridRebuilt();

    DeallocateRateCoefficients();
    AllocateRateCoefficients();

    return true;
}

/**
 * Rebuild rate coefficients.
 */
void IonRateEquation::Rebuild(
    const real_t, const real_t, FVM::UnknownQuantityHandler *unknowns
) {
    const len_t Nr = this->grid->GetNr();

    real_t *T = unknowns->GetUnknownData(id_T_cold);
    real_t *n = unknowns->GetUnknownData(id_n_cold);

    ADASRateInterpolator *acd = adas->GetACD(Zion);
    ADASRateInterpolator *scd = adas->GetSCD(Zion);

    // Iterate over charge state (0 ... Z)
    for (len_t Z0 = 0; Z0 <= Zion; Z0++) {
        for (len_t i = 0; i < Nr; i++) {
            Rec[Z0][i] = acd->Eval(Z0, n[i], T[i]);
            Ion[Z0][i] = scd->Eval(Z0, n[i], T[i]);
        }
    }

    // ///////
    // Construct fast-electron ionization rate
    real_t *n_hot = unknowns->GetUnknownData(id_n_hot);
    real_t *n_tot = unknowns->GetUnknownData(id_n_tot);

    // Evaluate 'Imp_i(r) = nhot(r) * I_i(r, n_tot, T_cold)'
    // Iterate over charge states (0 ... Z)
    for (len_t Z0 = 0; Z0 <= Zion; Z0++)
        for (len_t i = 0; i < Nr; i++)
            Imp[Z0][i] = scd->Eval(Z0, n_tot[i], T[i]) * n_hot[i];
}

/**
 * Build block of Jacobian matrix for the given charge state.
 *
 * derivId: ID of unknown quantity with respect to which differentiation
 *          should be carried out.
 * uqtyId:  ID of unknown quantity to differentiate.
 * jac:     Jacobian matrix to build.
 * x:       Current value of the unknown quantity.
 * iIon:    Index of ion to build jacobian for.
 * Z0:      Ion charge state.
 * rOffset: Offset in matrix block to set elements of.
 */
void IonRateEquation::SetCSJacobianBlock(
    const len_t /*derivId*/, const len_t /*uqtyId*/, FVM::Matrix * /*jac*/,
    const real_t* /*x*/,
    const len_t /*iIon*/, const len_t /*Z0*/, const len_t /*rOffset*/
) {
    throw NotImplementedException("Jacobian for ion rate equation not implemented.");
}

/**
 * Build linear operator matrix for this equation.
 *
 * mat:     Linear operator matrix.
 * rhs:     Vector representing the equation right-hand-side.
 * iIon:    Index of ion to build matrix for.
 * Z0:      Ion charge state.
 * rOffset: Offset in matrix block to set elements of.
 */
void IonRateEquation::SetCSMatrixElements(
    FVM::Matrix *mat, real_t*, const len_t iIon, const len_t Z0, const len_t rOffset
) {
    #define NI(J,V) \
        mat->SetElement(\
            rOffset+ir, rOffset+ir+(J)*Nr, \
            (V) \
        )
    #   include "IonRateEquation.set.cpp"
    #undef NI
}

/**
 * Build function vector.
 *
 * vec:     Function vector to set elements of.
 * nions:   Ion densities.
 * iIon:    Index of ion to build matrix for.
 * Z0:      Ion charge state.
 * rOffset: Offset in matrix block to set elements of.
 */
void IonRateEquation::SetCSVectorElements(
    real_t *vec, const real_t *nions,
    const len_t iIon, const len_t Z0, const len_t rOffset
) {
    #define NI(J,V) \
        vec[rOffset+ir] += (V) * nions[rOffset+ir+(J)*Nr]
    #   include "IonRateEquation.set.cpp"
    #undef NI
}


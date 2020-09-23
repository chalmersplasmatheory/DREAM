/**
 * Implementation of the ion rate equation
 *
 *   d n_i^(j) / dt =
 *      I_i^(j-1) n_i^(j-1) n_cold - I_i^(j) n_i^(j) n_cold  +
 *      R_i^(j+1) n_i^(j+1) n_cold - R_i^(j) n_i^(j) n_cold
 *
 * where
 *
 *   I_i^(j)  = ionization rate coefficient for charge state 'j'
 *              of ion species 'i'.
 *   R_i^(j)  = radiative recombination rate for charge state 'j'
 *              of ion species 'i'.
 *
 * Note that this equation is applied to a single _ion species_,
 * (and to all its charge states).
 * If using collfreq_mode FULL and kinetic ionization is used,
 * ionization rates will here be set to 0 to avoid double counting.
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
    ADAS *adas, FVM::UnknownQuantityHandler *unknowns,
    bool addFluidIonization
) : IonEquationTerm<FVM::EquationTerm>(g, ihdl, iIon), adas(adas), addFluidIonization(addFluidIonization) {
    
    this->unknowns  = unknowns;
    this->id_ions   = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    this->id_n_cold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
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
    this->PartialNRec = new real_t*[(Zion+1)];
    this->PartialTRec = new real_t*[(Zion+1)];
    this->Ion = new real_t*[(Zion+1)];
    this->PartialNIon = new real_t*[(Zion+1)];
    this->PartialTIon = new real_t*[(Zion+1)];

    this->Rec[0]         = new real_t[Nr*(Zion+1)];
    this->PartialNRec[0] = new real_t[Nr*(Zion+1)];
    this->PartialTRec[0] = new real_t[Nr*(Zion+1)];
    this->Ion[0]         = new real_t[Nr*(Zion+1)];
    this->PartialNIon[0] = new real_t[Nr*(Zion+1)];
    this->PartialTIon[0] = new real_t[Nr*(Zion+1)];

    for (len_t i = 1; i <= Zion; i++) {
        this->Rec[i]         = this->Rec[i-1] + Nr;
        this->PartialNRec[i] = this->PartialNRec[i-1] + Nr;
        this->PartialTRec[i] = this->PartialTRec[i-1] + Nr;
        this->Ion[i]         = this->Ion[i-1] + Nr;
        this->PartialNIon[i] = this->PartialNIon[i-1] + Nr;
        this->PartialTIon[i] = this->PartialTIon[i-1] + Nr;
    }
}

/**
 * Deallocate memory for the ionization rate coefficients.
 */
void IonRateEquation::DeallocateRateCoefficients() {
    delete [] this->Ion[0];
    delete [] this->PartialNIon[0];
    delete [] this->PartialTIon[0];
    delete [] this->Rec[0];
    delete [] this->PartialNRec[0];
    delete [] this->PartialTRec[0];

    delete [] this->Ion;
    delete [] this->PartialNIon;
    delete [] this->PartialTIon;
    delete [] this->Rec;
    delete [] this->PartialNRec;
    delete [] this->PartialTRec;
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
    for (len_t Z0 = 0; Z0 <= Zion; Z0++)
        for (len_t i = 0; i < Nr; i++) {
            Rec[Z0][i]         = acd->Eval(Z0, n[i], T[i]);
            PartialNRec[Z0][i] = acd->Eval_deriv_n(Z0, n[i], T[i]);
            PartialTRec[Z0][i] = acd->Eval_deriv_T(Z0, n[i], T[i]);
            Ion[Z0][i]         = 0;
            PartialNIon[Z0][i] = 0;
            PartialTIon[Z0][i] = 0;
        }
    // if not using kinetic ionization model, set fluid ionization rates
    if(addFluidIonization)
        for (len_t Z0 = 0; Z0 <= Zion; Z0++)
            for (len_t i = 0; i < Nr; i++) {
                Ion[Z0][i]         = scd->Eval(Z0, n[i], T[i]);
                PartialNIon[Z0][i] = scd->Eval_deriv_n(Z0, n[i], T[i]);
                PartialTIon[Z0][i] = scd->Eval_deriv_T(Z0, n[i], T[i]);
            }
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
    const len_t uqtyId, const len_t derivId, FVM::Matrix *jac,
    const real_t* nions,
    const len_t iIon, const len_t Z0, const len_t rOffset
) {
    if (derivId == uqtyId) {
        this->SetCSMatrixElements(jac, nullptr, iIon, Z0, rOffset);
    }

    #define NI(J,V) \
        jac->SetElement(\
            rOffset+ir, ir, \
            (V) * nions[rOffset+ir+(J)*Nr] \
        )

    if(derivId == id_T_cold) {
        #include "IonRateEquation.setDT.cpp"
    }

    if(derivId == id_n_cold){
        #include "IonRateEquation.setDN.cpp"        
    }
    #undef NI
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


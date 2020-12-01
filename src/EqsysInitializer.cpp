/**
 * Equation system initializer.
 */

#include <gsl/gsl_interp.h>
#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include "DREAM/EqsysInitializer.hpp"
#include "DREAM/IO.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Interpolator3D.hpp"


using namespace DREAM;
using namespace std;

const int_t
    EqsysInitializer::COLLQTYHDL_HOTTAIL=-1,
    EqsysInitializer::COLLQTYHDL_RUNAWAY=-2,
    EqsysInitializer::RUNAWAY_FLUID=-3;


/**
 * Constructor.
 */
EqsysInitializer::EqsysInitializer(
    FVM::UnknownQuantityHandler *unknowns,
    vector<UnknownQuantityEquation*> *unknown_equations,
    FVM::Grid *fluidGrid, FVM::Grid *hottailGrid, FVM::Grid *runawayGrid,
    enum OptionConstants::momentumgrid_type hottail_type,
    enum OptionConstants::momentumgrid_type runaway_type

) : unknowns(unknowns), unknown_equations(unknown_equations),
    fluidGrid(fluidGrid), hottailGrid(hottailGrid), runawayGrid(runawayGrid),
    hottail_type(hottail_type), runaway_type(runaway_type) { }

/**
 * Destructor.
 */
EqsysInitializer::~EqsysInitializer() { }

/**
 * Add a new initialization rule.
 *
 * uqtyId: ID of the unknown quantity to which the rule applies.
 * type:   Rule type (see 'enum initrule_t').
 */
void EqsysInitializer::AddRule(const string& uqtyName, const enum initrule_t type, initfunc_t fnc) {
    const int_t uqtyId = (int_t)this->unknowns->GetUnknownID(uqtyName);
    this->AddRule(uqtyId, type, fnc);
}
void EqsysInitializer::AddRule(const int_t uqtyId, const enum initrule_t type, initfunc_t fnc) {
    if (this->HasRuleFor(uqtyId))
        throw EqsysInitializerException(
            "An initialization rule for '%s' has already been added.",
            this->unknowns->GetUnknown(uqtyId)->GetName().c_str()
        );

    if (type == INITRULE_STEADY_STATE_SOLVE)
        throw NotImplementedException(
            "%s: No support for solving for initial steady state implemented yet.",
            this->unknowns->GetUnknown(uqtyId)->GetName().c_str()
        );

    this->rules[uqtyId] = new struct initrule(
        uqtyId, type, vector<int_t>(0), fnc
    );
}

/**
 * Execute all initialization rules and initialize the unknown
 * quantities of the equation system.
 *
 * t0: Time for which the system should be initialized.
 */
void EqsysInitializer::Execute(const real_t t0) {
    // List specifying order of execution
    vector<int_t> order;

    // Add default rules for special objects...
    int_t
        id_n_cold  = (int_t)this->unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD),
        id_n_tot   = (int_t)this->unknowns->GetUnknownID(OptionConstants::UQTY_N_TOT),
        id_T_cold  = (int_t)this->unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD),
        id_ions    = (int_t)this->unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES),
        id_E_field = (int_t)this->unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD);

    if (this->cqhHottail != nullptr)
        this->AddRule(COLLQTYHDL_HOTTAIL, INITRULE_EVAL_EQUATION, nullptr, id_n_cold, id_T_cold, id_ions);
    if (this->cqhRunaway != nullptr)
        this->AddRule(COLLQTYHDL_RUNAWAY, INITRULE_EVAL_EQUATION, nullptr, id_n_cold, id_T_cold, id_ions);
    this->AddRule(RUNAWAY_FLUID, INITRULE_EVAL_EQUATION, nullptr, id_n_cold, id_n_tot, id_T_cold, id_ions, id_E_field);

    // Build the 'order' list (sort the rules in order of
    // "least" to "most" dependent)
    for (auto rule : this->rules) {
        // Check if the rule has already been executed
        if (rule.second->executed)
            continue;

        vector<int_t> tmp = ConstructExecutionOrder(rule.second);
        order.insert(order.end(), tmp.begin(), tmp.end());
    }

    // Execute the rules in the calculated order
    for (int_t uqtyId : order) {
        struct initrule *rule = this->rules[uqtyId];

        // Verify that all dependences are initialized...
        /*for (int_t dep : rule->dependencies) {
            if (dep < 0) continue;

            if (!this->unknowns->HasInitialValue(dep))
                throw EqsysInitializerException(
                    "Bug in EqsysInitializer: dependency '%s' for '%s' not satisfied.",
                    this->unknowns->GetUnknown(dep)->GetName().c_str(),
                    this->unknowns->GetUnknown(uqtyId)->GetName().c_str()
                );
        }*/

        // Special object?
        if (uqtyId < 0) {
            switch (uqtyId) {
                case COLLQTYHDL_HOTTAIL:
                    this->cqhHottail->Rebuild();
                    break;

                case COLLQTYHDL_RUNAWAY:
                    this->cqhRunaway->Rebuild();
                    break;

                case RUNAWAY_FLUID:
                    this->runawayFluid->Rebuild();
                    break;

                default:
                    throw EqsysInitializerException(
                        "Unrecognizd ID of special object to initialize: " INT_T_PRINTF_FMT ".",
                        uqtyId
                    );
            }

        // Regular unknown quantities...
        } else {
            switch (rule->type) {
                case INITRULE_EVAL_EQUATION:
                    this->EvaluateEquation(t0, uqtyId);
                    break;

                case INITRULE_EVAL_FUNCTION:
                    this->EvaluateFunction(t0, uqtyId);
                    break;

                case INITRULE_STEADY_STATE_SOLVE:
                    throw NotImplementedException("The 'STEADY_STATE_SOLVE' initialization rule has not been implemented yet.");
                    break;

                default:
                    throw EqsysInitializerException(
                        "Unrecognized initialization rule type: %d.",
                        rule->type
                    );
            }
            // when ions have been initialized, initialize the ionHandler
            if(uqtyId==id_ions) 
                ionHandler->Rebuild();
        }
    }

    // Verfiy that all unknown quantities have now been initialized
    this->VerifyAllInitialized();
}

/**
 * Calculates the order in which quantities on which the specified
 * rule depends should be initialized. The ID of the specified rule
 * is added last to the execution order.
 *
 * rule: Rule to resolve dependency execution order for.
 */
vector<int_t> EqsysInitializer::ConstructExecutionOrder(struct initrule *rule) {
    vector<int_t> order;

    // Traverse dependencies and ensure
    for (int_t dep : rule->dependencies) {
        if (dep >= 0 && this->unknowns->HasInitialValue(dep)) {
            continue;
        } else if (this->HasRuleFor(dep)) {
            struct initrule *deprule = this->rules[dep];

            // If the rule has already been executed, we can
            // safely proceed to the next dependency
            if (deprule->executed)
                continue;
            // Seeing this rule a second time means we have
            // discovered a circular dependency...
            else if (deprule->marked) {
                throw EqsysInitializerException(
                    "Unable to resolve circular dependency: %s -> %s.",
                    this->unknowns->GetUnknown(rule->uqtyId)->GetName().c_str(),
                    this->unknowns->GetUnknown(deprule->uqtyId)->GetName().c_str()
                );
            } else {
                // This rule has no un-initialized dependencies so we add
                // it to the list of rules.
                deprule->marked = true;
                vector<int_t> tmp = ConstructExecutionOrder(deprule);
                order.insert(order.end(), tmp.begin(), tmp.end());
            }
        } else
            throw EqsysInitializerException(
                "Unable to resolve initialization dependencies. No rule to initialize '%s'.",
                this->unknowns->GetUnknown(dep)->GetName().c_str()
            );
    }

    rule->executed = true;
    order.push_back(rule->uqtyId);
    return order;
}

/**
 * Check an initialization rule for the specified unknown
 * quantity exists.
 */
bool EqsysInitializer::HasRuleFor(const int_t uqtyId) const {
    return (this->rules.find(uqtyId) != this->rules.end());
}

/**
 * Initialize unknowns from DREAM output stored in the given
 * output file.
 *
 * This method will remove initialization rules for unknown
 * quantities which can be initialized from the output. Unknown
 * quantities which are not available among the output will be
 * initialized according to their rules.
 *
 * filename:   Name of file to load quantities from.
 * sf:         SFile object to use for loading output quantities.
 * t0:         Time for which to initialize the system (does not
 *             necessarily correspond to the time from which to load
 *             the data from; t0 denotes the time in the current simulation,
 *             which can have a completely different offset).
 * tidx:       Index of time slice to initialize simulation from
 *             (if negative, it is transformed to "nt+tidx", so
 *             that '-1' accesses the very last time step).
 * ignoreList: List of unknown quantities to ignore and NOT initialize
 *             from the given output file.
 */
void EqsysInitializer::InitializeFromOutput(
    const string& filename, const real_t t0, int_t tidx, IonHandler *ions,
    vector<string>& ignoreList
) {
    SFile *sf = SFile::Create(filename, SFILE_MODE_READ);

    this->InitializeFromOutput(sf, t0, tidx, ions, ignoreList);

    sf->Close();
    delete sf;
}
void EqsysInitializer::InitializeFromOutput(
    SFile *sf, const real_t t0, int_t tidx, IonHandler *ions,
    vector<string>& ignoreList
) {
    sfilesize_t nr, nt, np1_hot, np2_hot, np1_re, np2_re;
    enum OptionConstants::momentumgrid_type
        momtype_hot, momtype_re;

    // Load grids
    real_t *r = sf->GetList("grid/r", &nr);
    real_t *t = sf->GetList("grid/t", &nt);
    real_t *hot_p1 = nullptr, *hot_p2 = nullptr;
    real_t *re_p1  = nullptr, *re_p2  = nullptr;

    if (sf->HasVariable("grid/hottail/p1")) {
        hot_p1      = sf->GetList("grid/hottail/p1", &np1_hot);
        hot_p2      = sf->GetList("grid/hottail/p2", &np2_hot);
        momtype_hot = (enum OptionConstants::momentumgrid_type)sf->GetInt("grid/hottail/type");
    }
    if (sf->HasVariable("grid/runaway/p1")) {
        re_p1      = sf->GetList("grid/runaway/p1", &np1_re);
        re_p2      = sf->GetList("grid/runaway/p2", &np2_re);
        momtype_re = (enum OptionConstants::momentumgrid_type)sf->GetInt("grid/hottail/type");
    }

    // Shift time index if negative
    if (tidx < 0)
        tidx = nt+tidx;

    // Iterate over unknown quantities
    const int_t nUnknowns = (int_t)this->unknowns->GetNUnknowns();
    for (int_t i = 0; i < nUnknowns; i++) {
        FVM::UnknownQuantity *uqn = this->unknowns->GetUnknown(i);
        string name = "eqsys/" + uqn->GetName();

        // If the unknown quantity is in the ignore list, quietly
        // skip it...
        if (find(ignoreList.begin(), ignoreList.end(), uqn->GetName()) != ignoreList.end())
            continue;

        // Check if unknown quantity is available in output...
        if (!sf->HasVariable(name))
            continue;
        else
            DREAM::IO::PrintInfo(
                "Initializing '%s' from '%s'...",
                uqn->GetName().c_str(), sf->filename.c_str()
            );

        // Get data
        sfilesize_t dims[4], ndims;
        real_t *data = sf->GetMultiArray_linear(name, 4, ndims, dims);

        // Time + radius
        if (ndims == 2) {
            this->__InitTR(uqn, t0, tidx, nr, r, data, dims);

        // Time + multiples + radius
        } else if (ndims == 3) {
            // If ions, verify that the ion density structure has not changed...
            if (uqn->GetName() == OptionConstants::UQTY_ION_SPECIES) {
                sfilesize_t n;
                int64_t *Z = sf->GetIntList("ionmeta/Z", &n);

                if (ions->GetNZ() != n)
                    throw EqsysInitializerException(
                        "%s: ion density structure has changed from previous "
                        "simulation. More ion species have been added.",
                        name.c_str()
                    );

                for (int_t j = 0; j < (int_t)n; j++) {
                    if (ions->GetZ(j) != (len_t)Z[j])
                        throw EqsysInitializerException(
                            "%s: ion density structure has changed from previous "
                            "simulation. Ion at index " INT_T_PRINTF_FMT " has changed.",
                            name.c_str(), j
                        );
                }

                delete [] Z;
            }

            this->__InitTRmult(uqn, t0, tidx, nr, r, data, dims);

        // Time + radius + momentum
        } else if (ndims == 4) {
            // Hot-tail
            if (uqn->GetGrid() == this->hottailGrid) {
                if (nr != dims[1])
                    throw EqsysInitializerException("Initializing from output '%s': invalid size of dimension 1. Size was " LEN_T_PRINTF_FMT ", expected " LEN_T_PRINTF_FMT ".", name.c_str(), dims[1], nr);
                else if (np2_hot != dims[2])
                    throw EqsysInitializerException("Initializing from output '%s': invalid size of dimension 2. Size was " LEN_T_PRINTF_FMT ", expected " LEN_T_PRINTF_FMT ".", name.c_str(), dims[2], np2_hot);
                else if (np1_hot != dims[3])
                    throw EqsysInitializerException("Initializing from output '%s': invalid size of dimension 3. Size was " LEN_T_PRINTF_FMT ", expected " LEN_T_PRINTF_FMT ".", name.c_str(), dims[3], np1_hot);

                this->__InitTR2P(
                    uqn, t0, tidx, r, hot_p1, hot_p2,
                    data, dims, momtype_hot, this->hottail_type
                );
            // Runaway
            } else if (uqn->GetGrid() == this->runawayGrid) {
                if (nr != dims[1])
                    throw EqsysInitializerException("Initializing from output '%s': invalid size of dimension 1. Size was " LEN_T_PRINTF_FMT ", expected " LEN_T_PRINTF_FMT ".", name.c_str(), dims[1], nr);
                else if (np2_re != dims[2])
                    throw EqsysInitializerException("Initializing from output '%s': invalid size of dimension 2. Size was " LEN_T_PRINTF_FMT ", expected " LEN_T_PRINTF_FMT ".", name.c_str(), dims[2], np2_re);
                else if (np1_re != dims[3])
                    throw EqsysInitializerException("Initializing from output '%s': invalid size of dimension 3. Size was " LEN_T_PRINTF_FMT ", expected " LEN_T_PRINTF_FMT ".", name.c_str(), dims[3], np1_re);

                this->__InitTR2P(
                    uqn, t0, tidx, r, re_p1, re_p2,
                    data, dims, momtype_re, this->runaway_type
                );
            } else
                throw EqsysInitializerException(
                    "Initializing from output: '%s': unrecognized momentum grid for quantity.",
                    name.c_str()
                );
        } else
            throw EqsysInitializerException(
                "Initializing from output: '%s': unrecognized dimensions of quantity: "
                LEN_T_PRINTF_FMT,
                name.c_str(), ndims
            );

        // Remove initialization rule
        this->RemoveRule(i);

        delete [] data;
    }

    delete [] r;
    delete [] t;

    if (hot_p1 != nullptr) delete [] hot_p1;
    if (hot_p2 != nullptr) delete [] hot_p2;
    if (re_p1  != nullptr) delete [] re_p1;
    if (re_p2  != nullptr) delete [] re_p2;
}

/**
 * Initialize the given unknown quantity from the given
 * spatiotemporal (time+radius) data.
 *
 * uqn:  Unknown quantity to set initial value of.
 * t0:   Time for which the quantity should be initialized.
 * tidx: Index of time point in data to initialize from.
 * nr:   Number of radial grid points.
 * r:    Radial grid for given data.
 * data: Data to initialize from.
 * dims: Dimensions of the given array.
 */
void EqsysInitializer::__InitTR(
    FVM::UnknownQuantity *uqn, const real_t t0, const int_t tidx,
    const len_t nr, const real_t *r, const real_t *data,
    const sfilesize_t *dims
) {
    const len_t /*nt = dims[0],*/ nrel = dims[1];

    // Handle scalar quantities correctly
    const len_t nmult = (nrel==1 ? 1 : nrel/nr);
    const len_t NR = (nrel==1 ? 1 : uqn->GetGrid()->GetNr());

    // Get requested time step...
    const real_t *d = data + tidx*nrel;

    real_t *intpdata;

    if (nr % nrel != 0)
        throw EqsysInitializerException(
            "Initializing from output: '%s': dimensions mismatch. nrel is not a multiple of nr.",
            uqn->GetName().c_str()
        );

    // Scalar quantities and un-interpolatables
    if (nr == 1 || nrel == 1) {
        intpdata = new real_t[nmult*NR];

        // More than one value per radius
        // (as in the case of ion densities)?
        if (nmult == 1) {
            for (len_t i = 0; i < NR; i++)
                intpdata[i] = d[0];
        } else {
            for (len_t j = 0; j < nmult; j++)
                for (len_t i = 0; i < NR; i++)
                    intpdata[j*NR + i] = d[j];
        }

    // Interpolate in radius
    } else {
        // Regular fluid quantities...
        if (nmult == 1) {
            intpdata = SimulationGenerator::InterpolateR(
                nr, r, d, uqn->GetGrid()->GetRadialGrid(),
                gsl_interp_linear
            );

        // Ions...
        } else {
            intpdata = SimulationGenerator::InterpolateIonR(
                uqn->GetGrid()->GetRadialGrid(), nr, nmult,
                r, d, gsl_interp_linear
            );
        }
    }

    uqn->SetInitialValue(intpdata, t0);

    delete [] intpdata;
}

/**
 * Initialize the given unknown quantity from the given
 * spatiotemporal (time+radius) data. The data is assumed
 * to consist of several multiples of the time+radius grid
 * (e.g. ion species/charge state data).
 *
 * uqn:        Unknown quantity to set initial value of.
 * t0:         Time for which the quantity should be initialized.
 * tidx:       Index of time point in data to initialize from.
 * nMultiples: Number of multiples in the data, e.g. the number
 *             of ion charge states.
 * nr:         Number of radial grid points.
 * r:          Radial grid for given data.
 * data:       Data to initialize from.
 * dims:       Dimensions of the given array.
 */
void EqsysInitializer::__InitTRmult(
    FVM::UnknownQuantity *uqn, const real_t t0, const int_t tidx,
    const len_t nr, const real_t *r,
    const real_t *data, const sfilesize_t *dims
) {
    const len_t /*nt = dims[0],*/ nmult = dims[1], nrel = dims[2];
    const len_t NR = uqn->GetGrid()->GetNr();
    // Get requested time step...
    const real_t *d = data + tidx*nrel*nmult;
    //const real_t *d = data + tidx*NR*nmult;

    real_t *intpdata;

    if (nr != dims[2])
        throw EqsysInitializerException(
            "Initializing from output: '%s': dimensions mismatch. dims[2] != nr "
            "(" LEN_T_PRINTF_FMT " != " LEN_T_PRINTF_FMT ").",
            uqn->GetName().c_str()
        );

    // Scalar quantities and un-interpolatables
    if (nr == 1) {
        intpdata = new real_t[nmult*NR];
        for (len_t j = 0; j < nmult; j++)
            for (len_t i = 0; i < NR; i++)
                intpdata[j*NR + i] = d[j];

    // Interpolate in radius
    } else {
        intpdata = SimulationGenerator::InterpolateIonR(
            uqn->GetGrid()->GetRadialGrid(), nr, nmult,
            r, d, gsl_interp_linear
        );
    }

    uqn->SetInitialValue(intpdata, t0);

    delete [] intpdata;
}

/**
 * Initialize the given unknown quantity from the given
 * phase space (time, radius and momentum) data.
 *
 * uqn:         Unknown quantity to set initial value of.
 * t0:          Time for which the quantity should be initialized.
 * tidx:        Index of time point in data to initialize from.
 * r:           Radial grid for given data.
 * p1:          Grid for first momentum parameter on which 'data' is defined.
 * p2:          Grid for second momentum parameter on which 'data' is defined.
 * data:        Data to initialize from.
 * dims:        Dimensions of the given array.
 * momtype:     Type of momentum coordinates in 'p1' and 'p2'.
 * uqn_momtype: Type of momentum coordinates for grid used by 'uqn'.
 */
void EqsysInitializer::__InitTR2P(
    FVM::UnknownQuantity *uqn, const real_t t0, const int_t tidx,
    const real_t *r, const real_t *p1, const real_t *p2, const real_t *data,
    const sfilesize_t *dims,
    enum OptionConstants::momentumgrid_type momtype,
    enum OptionConstants::momentumgrid_type uqn_momtype
) {
    const len_t
        /*nt  = dims[0],*/
        nr  = dims[1],
        np2 = dims[2],
        np1 = dims[3];
    const real_t
        *d = data + tidx*(nr*np1*np2);

    enum FVM::Interpolator3D::interp_method interp_meth =
        FVM::Interpolator3D::INTERP_LINEAR;

    enum FVM::Interpolator3D::momentumgrid_type _momtype = (
        momtype == OptionConstants::MOMENTUMGRID_TYPE_PXI ?
            FVM::Interpolator3D::GRID_PXI :
            FVM::Interpolator3D::GRID_PPARPPERP
    );
    enum FVM::Interpolator3D::momentumgrid_type _uqn_momtype = (
        uqn_momtype == OptionConstants::MOMENTUMGRID_TYPE_PXI ?
            FVM::Interpolator3D::GRID_PXI :
            FVM::Interpolator3D::GRID_PPARPPERP
    );

    // Interpolate onto the grid of 'uqn'...
    FVM::Interpolator3D *interp = new FVM::Interpolator3D(
        nr, np2, np1, r, p2, p1, d, _momtype, interp_meth, false
    );

    const real_t *intp = interp->Eval(uqn->GetGrid(), _uqn_momtype);
    uqn->SetInitialValue(intp, t0);

    delete [] intp;
    delete interp;
}

/**
 * Remove any initialization rule for the specified quantity
 * if it exists.
 *
 * uqtyId: ID of unknown quantity to remove rule for.
 */
void EqsysInitializer::RemoveRule(const int_t uqtyId) {
    const auto &it = this->rules.find(uqtyId);

    // If the quantity has not init rule, just return...
    if (it == this->rules.end())
        return;

    this->rules.erase(it);
}

/**
 * Verify that all unknown quantities in the UnknownQuantityHandler
 * have been successfully initialized.
 */
void EqsysInitializer::VerifyAllInitialized() const {
    bool initialized = true;
    for (auto it : this->unknowns->GetUnknowns()) {
        if (!it->HasInitialValue()) {
            initialized = false;
            DREAM::IO::PrintError(
                "%s: No initial value defined for unknown quantity.",
                it->GetName().c_str()
            );
        }
    }

    if (!initialized)
        throw EqsysInitializerException(
            "Some quantities have no defined initial values."
        );
}


/***********************
 * EVALUATION ROUTINES *
 ***********************/
/**
 * Initialize the specified unknown quantity by evaluating
 * its equation.
 *
 * t0:     Time for which the system should be initialized.
 * uqtyId: ID of unknown quantity to evaluate equation for.
 */
void EqsysInitializer::EvaluateEquation(const real_t t0, const int_t uqtyId) {
    UnknownQuantityEquation *eqn = this->unknown_equations->at(uqtyId);

    if (!eqn->IsEvaluable())
        throw EqsysInitializerException(
            "Unable to initialize '%s': equation is not evaluable.",
            this->unknowns->GetUnknown(uqtyId)->GetName().c_str()
        );
    
    // Construct vector for temporary storage of quantity value
    const len_t N = eqn->NumberOfElements();
    real_t *vec = new real_t[N];
    for (len_t i = 0; i < N; i++)
        vec[i] = 0.0;

    // Evaluate equation
    eqn->RebuildEquations(t0, 0, unknowns);
    eqn->Evaluate(uqtyId, vec, unknowns);

    // Store initial value
    this->unknowns->SetInitialValue(uqtyId, vec, t0);

    delete [] vec;
}

/**
 * Initialize the specified unknown quantity by evaluating
 * its initialization function.
 *
 * t0:     Time for which the system should be initialized.
 * uqtyId: ID of the unknown quantity to initialize.
 */
void EqsysInitializer::EvaluateFunction(const real_t t0, const int_t uqtyId) {
    struct initrule *rule = this->rules[uqtyId];
    UnknownQuantityEquation *eqn = this->unknown_equations->at(uqtyId);

    if (!rule->init)
        throw EqsysInitializerException(
            "Unable to initialize '%s': no initialization function given.",
            this->unknowns->GetUnknown(uqtyId)->GetName().c_str()
        );

    // Construct vector for temporary storage of quantity value
    const len_t N = eqn->NumberOfElements();
    real_t *vec = new real_t[N];
    for (len_t i = 0; i < N; i++)
        vec[i] = 0.0;
    
    // Evaluate the unknown quantity
    rule->init(this->unknowns, vec);

    // Store initial value
    this->unknowns->SetInitialValue(uqtyId, vec, t0);

    delete [] vec;
}



#include <vector>
#include <string>
#include "DREAM/EquationSystem.hpp"
#include "DREAM/OtherQuantityHandler.hpp"
#include "DREAM/PostProcessor.hpp"
#include "DREAM/Settings/Settings.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"


using namespace DREAM;
using namespace std;


#define EQUATIONSYSTEM "eqsys"
#define INITIALIZATION "init"

/**
 * Define the options which can be set for things
 * related to the equation system.
 *
 * s: Settings object to define the options for.
 */
void SimulationGenerator::DefineOptions_EquationSystem(Settings *s) {
    s->DefineSetting(EQUATIONSYSTEM "/n_cold/type", "Type of equation to use for determining the cold electron density", (int_t)OptionConstants::UQTY_N_COLD_EQN_PRESCRIBED);
    DefineDataRT(EQUATIONSYSTEM "/n_cold", s);
    
//    s->DefineSetting(EQUATIONSYSTEM "/T_cold/type", "Type of equation to use for determining the electron temperature evolution", (int_t)OptionConstants::UQTY_T_COLD_EQN_PRESCRIBED);
//    DefineDataRT(EQUATIONSYSTEM "/T_cold", s);
}

/**
 * Define options for initialization.
 */
void SimulationGenerator::DefineOptions_Initializer(Settings *s) {
    s->DefineSetting(INITIALIZATION "/eqsysignore", "List of unknown quantities to NOT initialize from output file.", (const string)"");
    s->DefineSetting(INITIALIZATION "/filetimeindex", "Time index to take initialization data for from output file.", (int_t)-1);
    s->DefineSetting(INITIALIZATION "/fromfile", "Name of DREAM output file from which simulation should be initialized.", (const string)"");
    s->DefineSetting(INITIALIZATION "/t0", "Simulation at which to initialize the simulation.", (real_t)0.0);

	s->DefineSetting(INITIALIZATION "/solver_maxiter", "Maximum number of iterations for non-linear steady-state solver.", (int_t)100);
	s->DefineSetting(INITIALIZATION "/solver_reltol", "Relative tolerance used for non-linear steady-state solver.", (real_t)1e-6);
	s->DefineSetting(INITIALIZATION "/solver_verbose", "Whether or not to print convergence information for non-linear steady-state solver.", (bool)false);
	s->DefineSetting(INITIALIZATION "/solver_linear", "Primary linear solver to use.", (int_t)OptionConstants::LINEAR_SOLVER_LU);
	s->DefineSetting(INITIALIZATION "/solver_backup", "Secondary linear solver to use.", (int_t)OptionConstants::LINEAR_SOLVER_NONE);
}

/**
 * Construct an equation system, based on the specification
 * in the given 'Settings' object.
 *
 * s:           Settings object specifying how to construct
 *              the equation system.
 * fluidGrid:   Radial grid for the computation.
 * ht_type:     Exact type of the hot-tail momentum grid.
 * hottailGrid: Grid on which the hot-tail electron population
 *              is computed.
 * re_type:     Exact type of the runaway momentum grid.
 * runawayGrid: Grid on which the runaway electron population
 *              is computed.
 *
 * NOTE: The 'hottailGrid' and 'runawayGrid' will be 'nullptr'
 *       if disabled.
 */
EquationSystem *SimulationGenerator::ConstructEquationSystem(
    Settings *s, FVM::Grid *scalarGrid, FVM::Grid *fluidGrid,
    enum OptionConstants::momentumgrid_type ht_type, FVM::Grid *hottailGrid,
    enum OptionConstants::momentumgrid_type re_type, FVM::Grid *runawayGrid,
    ADAS *adas, NIST *nist, AMJUEL *amjuel
) {
    EquationSystem *eqsys = new EquationSystem(scalarGrid, fluidGrid, ht_type, hottailGrid, re_type, runawayGrid, s);
    struct OtherQuantityHandler::eqn_terms *oqty_terms = new OtherQuantityHandler::eqn_terms;

    // Timing information
    eqsys->SetTiming(s->GetBool("/output/timingstdout"), s->GetBool("/output/timingfile"));

    // Initialize from previous simulation output?
    const real_t t0 = ConstructInitializer(eqsys, s);

    // Construct unknowns
    ConstructUnknowns(eqsys, s, scalarGrid, fluidGrid, hottailGrid, runawayGrid);


    // Construct equations according to settings
    ConstructEquations(eqsys, s, adas, nist, amjuel, oqty_terms);

    // Construct the "other" quantity handler
    ConstructOtherQuantityHandler(eqsys, s, oqty_terms);

    // Figure out which unknowns must be part of the matrix,
    // and set initial values for those quantities which don't
    // yet have an initial value.
    eqsys->ProcessSystem(t0);

    // (these must be initialized AFTER calling 'ProcessSystem()' on
    // the equation system, since we need to which unknowns are
    // "non-trivial", i.e. need to show up in the solver matrices,
    // in order to build them)
    
    // Construct the time stepper
    ConstructTimeStepper(eqsys, s);

    // Construct solver (must be done after processing equation system,
    // since we need to know which unknowns are "non-trivial",
    // i.e. need to show up in the solver matrices)
    ConstructSolver(eqsys, s);

    return eqsys;
}

/**
 * Set the equations of the equation system.
 *
 * eqsys: Equation system to define quantities in.
 * s:     Settings object specifying how to construct
 *        the equation system.
 * adas:  ADAS database object.
 * nist:  NIST database object.
 * amjuel: AMJUEL database object.
 *
 * NOTE: The 'hottailGrid' and 'runawayGrid' will be 'nullptr'
 *       if disabled.
 */
void SimulationGenerator::ConstructEquations(
    EquationSystem *eqsys, Settings *s, ADAS *adas, NIST *nist, AMJUEL *amjuel,
    struct OtherQuantityHandler::eqn_terms *oqty_terms
) {
    FVM::Grid *hottailGrid = eqsys->GetHotTailGrid();
    FVM::Grid *runawayGrid = eqsys->GetRunawayGrid();
    FVM::Grid *fluidGrid   = eqsys->GetFluidGrid();
    FVM::UnknownQuantityHandler *unknowns = eqsys->GetUnknownHandler();
    enum OptionConstants::momentumgrid_type ht_type = eqsys->GetHotTailGridType();
    enum OptionConstants::momentumgrid_type re_type = eqsys->GetRunawayGridType();
    SPIHandler* SPI;
    enum OptionConstants::eqterm_spi_ablation_mode spi_ablation_mode = (enum OptionConstants::eqterm_spi_ablation_mode)s->GetInteger("eqsys/spi/ablation");
    if(spi_ablation_mode!=OptionConstants::EQTERM_SPI_ABLATION_MODE_NEGLECT){
        SPI = ConstructSPIHandler(fluidGrid, unknowns, s);
        eqsys->SetSPIHandler(SPI);
    }
    
    // Fluid equations
    ConstructEquation_Ions(eqsys, s, adas, amjuel, oqty_terms);


    IonHandler *ionHandler = eqsys->GetIonHandler();
    // Construct collision quantity handlers
    if (hottailGrid != nullptr) {
        CollisionQuantityHandler *cqh = ConstructCollisionQuantityHandler(ht_type, hottailGrid, unknowns, ionHandler, s);
        eqsys->SetHotTailCollisionHandler(cqh);
    }
	if (runawayGrid != nullptr) {
        CollisionQuantityHandler *cqh = ConstructCollisionQuantityHandler(re_type, runawayGrid, unknowns, ionHandler, s);
        eqsys->SetRunawayCollisionHandler(cqh);
    }
    ConstructRunawayFluid(fluidGrid,unknowns,ionHandler,re_type,eqsys,s);
    eqsys->SetREFluid(eqsys->GetREFluid());
    // Post processing handler
    FVM::MomentQuantity::pThresholdMode pMode = FVM::MomentQuantity::P_THRESHOLD_MODE_MIN_THERMAL;
    real_t pThreshold = 0.0;
    enum OptionConstants::collqty_collfreq_mode collfreq_mode =
        (enum OptionConstants::collqty_collfreq_mode)s->GetInteger("collisions/collfreq_mode");
    if(eqsys->HasHotTailGrid() && collfreq_mode == OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL){
        pThreshold = (real_t)s->GetReal("eqsys/f_hot/pThreshold");
        pMode = (FVM::MomentQuantity::pThresholdMode)s->GetInteger("eqsys/f_hot/pThresholdMode");
    }
    PostProcessor *postProcessor = new PostProcessor(fluidGrid, unknowns, pThreshold, pMode);
    eqsys->SetPostProcessor(postProcessor);

    // Hot-tail quantities
    if (eqsys->HasHotTailGrid()) {
        ConstructEquation_f_hot(eqsys, s, oqty_terms);
    }

    // Runaway quantities
    FVM::Operator *transport_fre = nullptr;
    if (eqsys->HasRunawayGrid()) {
        ConstructEquation_f_re(eqsys, s, oqty_terms, &transport_fre);
    }
    ConstructEquation_E_field(eqsys, s, oqty_terms);
    ConstructEquation_j_hot(eqsys, s);
    ConstructEquation_j_tot(eqsys, s);
    ConstructEquation_j_ohm(eqsys, s);
    ConstructEquation_j_re(eqsys, s);
    ConstructEquation_n_cold(eqsys, s);
    ConstructEquation_n_hot(eqsys, s);
    ConstructEquation_T_cold(eqsys, s, adas, nist, amjuel, oqty_terms);
    
    if(spi_ablation_mode==OptionConstants::EQTERM_SPI_ABLATION_MODE_NGPS){
		ConstructEquation_Ions_abl(eqsys, s, adas, amjuel);
		ConstructEquation_n_abl(eqsys, s);
		ConstructEquation_T_abl(eqsys, s, adas, nist, amjuel, oqty_terms);
	}

    if(eqsys->GetSPIHandler()!=nullptr){
        ConstructEquation_SPI(eqsys,s);
        if(hottailGrid != nullptr){
        	ConstructEquation_W_hot(eqsys,s);
        	ConstructEquation_q_hot(eqsys,s);
        }
    }

    // Add equations for net ion density of each species and its energy density
    // only if including the cross-species collisional energy transfer
    OptionConstants::uqty_T_i_eqn typeTi = (OptionConstants::uqty_T_i_eqn) s->GetInteger("eqsys/n_i/typeTi");
    if(typeTi == OptionConstants::UQTY_T_I_INCLUDE /* && typeTcold == OptionConstants::UQTY_T_COLD_SELF_CONSISTENT */){
        ConstructEquation_Ion_Ni(eqsys,s);
        ConstructEquation_T_i(eqsys,s);
    }
    // NOTE: The runaway number may depend explicitly on
    // either f_hot or f_re and must therefore be constructed
    // AFTER the calls to 'ConstructEquation_f_hot()' and
    // 'ConstructEquation_f_re()'.
    ConstructEquation_n_re(eqsys, s, oqty_terms, transport_fre);

    ConstructEquation_psi_p(eqsys, s);
    ConstructEquation_psi_edge(eqsys, s);
    
    // Helper quantities
    ConstructEquation_n_tot(eqsys, s);
    OptionConstants::eqterm_hottail_mode hottail_mode = (enum OptionConstants::eqterm_hottail_mode)s->GetInteger("eqsys/n_re/hottail");
    OptionConstants::uqty_f_hot_dist_mode ht_dist_mode = (enum OptionConstants::uqty_f_hot_dist_mode)s->GetInteger("eqsys/f_hot/dist_mode");    
    if(hottail_mode != OptionConstants::EQTERM_HOTTAIL_MODE_DISABLED && ht_dist_mode == OptionConstants::UQTY_F_HOT_DIST_MODE_NONREL){
        ConstructEquation_tau_coll(eqsys);
    }

}

/**
 * Load initialization settings for the EquationSystem.
 *
 * eqsys:       Equation system to define quantities in.
 * s:           Settings object specifying how to construct
 *              the equation system.
 */
real_t SimulationGenerator::ConstructInitializer(
    EquationSystem *eqsys, Settings *s
) {
    real_t t0 = s->GetReal(INITIALIZATION "/t0");
    const string& filename = s->GetString(INITIALIZATION "/fromfile");
    int_t timeIndex = s->GetInteger(INITIALIZATION "/filetimeindex");

    // Initialize from previous output
    if (filename != "") {
        vector<string> ignoreList = s->GetStringList(INITIALIZATION "/eqsysignore");
        eqsys->SetInitializerFile(filename, ignoreList, timeIndex);
    }

	len_t maxiter = (len_t)s->GetInteger(INITIALIZATION "/solver_maxiter");
	real_t reltol = s->GetReal(INITIALIZATION "/solver_reltol");
	bool verbose  = s->GetBool(INITIALIZATION "/solver_verbose");
	enum OptionConstants::linear_solver linear_solver =
		(enum OptionConstants::linear_solver)s->GetInteger(INITIALIZATION "/solver_linear");
	enum OptionConstants::linear_solver backup_solver =
		(enum OptionConstants::linear_solver)s->GetInteger(INITIALIZATION "/solver_backup");

	eqsys->SetInitializerSolver(maxiter, reltol, linear_solver, backup_solver, verbose);

    return t0;
}

/**
 * Construct the unknowns of the equation system.
 *
 * eqsys:       Equation system to define quantities in.
 * s:           Settings object specifying how to construct
 *              the equation system.
 * fluidGrid:   Radial grid for the computation.
 * hottailGrid: Grid on which the hot-tail electron population
 *              is computed.
 * runawayGrid: Grid on which the runaway electron population
 *              is computed.
 *
 * NOTE: The 'hottailGrid' and 'runawayGrid' will be 'nullptr'
 *       if disabled.
 */
void SimulationGenerator::ConstructUnknowns(
    EquationSystem *eqsys, Settings *s, FVM::Grid *scalarGrid, FVM::Grid *fluidGrid,
    FVM::Grid *hottailGrid, FVM::Grid *runawayGrid
) {
    #define DEFU_HOT(NAME) eqsys->SetUnknown( \
        OptionConstants::UQTY_ ## NAME, \
        OptionConstants::UQTY_ ## NAME ## _DESC, \
        hottailGrid)
    #define DEFU_RE(NAME) eqsys->SetUnknown( \
        OptionConstants::UQTY_ ## NAME, \
        OptionConstants::UQTY_ ## NAME ## _DESC, \
        runawayGrid)
    #define DEFU_FLD(NAME) eqsys->SetUnknown( \
        OptionConstants::UQTY_ ## NAME, \
        OptionConstants::UQTY_ ## NAME ## _DESC, \
        fluidGrid)
    #define DEFU_FLD_N(NAME,NMULT) eqsys->SetUnknown( \
        OptionConstants::UQTY_ ## NAME, \
        OptionConstants::UQTY_ ## NAME ## _DESC, \
        fluidGrid, (NMULT))
    #define DEFU_SCL(NAME) eqsys->SetUnknown( \
        OptionConstants::UQTY_ ## NAME, \
        OptionConstants::UQTY_ ## NAME ## _DESC, \
        scalarGrid)
    #define DEFU_SCL_N(NAME,NMULT) eqsys->SetUnknown( \
        OptionConstants::UQTY_ ## NAME, \
        OptionConstants::UQTY_ ## NAME ## _DESC, \
        scalarGrid,(NMULT))

    // Hot-tail quantities
    if (hottailGrid != nullptr) {
        DEFU_HOT(F_HOT);
    }

    // Runaway quantities
    if (runawayGrid != nullptr) {
        DEFU_RE(F_RE);
    }

    // Fluid quantities
    len_t nIonChargeStates = GetNumberOfIonChargeStates(s);
    DEFU_FLD_N(ION_SPECIES, nIonChargeStates);
    DEFU_FLD(N_HOT);
    DEFU_FLD(N_COLD);
    DEFU_FLD(N_RE);
    DEFU_FLD(J_OHM);
    DEFU_FLD(J_HOT);
    DEFU_FLD(J_RE);
    DEFU_FLD(J_TOT);
    DEFU_FLD(T_COLD);
    DEFU_FLD(W_COLD);
    DEFU_FLD(E_FIELD);
    DEFU_FLD(POL_FLUX);
    DEFU_SCL(PSI_EDGE);
    DEFU_SCL(PSI_WALL);
    DEFU_SCL(I_P);

    if (s->GetBool("eqsys/n_re/negative_re"))
        DEFU_FLD(N_RE_NEG);

    enum OptionConstants::eqterm_spi_ablation_mode spi_ablation_mode = (enum OptionConstants::eqterm_spi_ablation_mode)s->GetInteger("eqsys/spi/ablation");
    if(spi_ablation_mode!=OptionConstants::EQTERM_SPI_ABLATION_MODE_NEGLECT){
        len_t nShard;
        s->GetRealArray("eqsys/spi/init/rp", 1, &nShard);
        DEFU_SCL_N(Y_P,nShard);
        DEFU_SCL_N(X_P,3*nShard);
        DEFU_SCL_N(V_P,3*nShard);
        
        if (hottailGrid != nullptr){
        	DEFU_FLD(Q_HOT);
        	DEFU_FLD(W_HOT);
    	}
    	if(spi_ablation_mode==OptionConstants::EQTERM_SPI_ABLATION_MODE_NGPS){
    		DEFU_FLD_N(ION_SPECIES_ABL, nIonChargeStates);
    		DEFU_FLD(N_ABL);
    		DEFU_FLD(T_ABL);
    		DEFU_FLD(W_ABL);
		}
    }

    if( (OptionConstants::uqty_T_i_eqn)s->GetInteger("eqsys/n_i/typeTi") == OptionConstants::UQTY_T_I_INCLUDE ){
        len_t nIonSpecies = GetNumberOfIonSpecies(s);
        DEFU_FLD_N(WI_ENER, nIonSpecies);
        DEFU_FLD_N(NI_DENS, nIonSpecies);
    }
    
    // Fluid helper quantities
    DEFU_FLD(N_TOT);
    if (hottailGrid != nullptr){
        DEFU_FLD(S_PARTICLE);
    }
    OptionConstants::eqterm_hottail_mode hottail_mode = (enum OptionConstants::eqterm_hottail_mode)s->GetInteger("eqsys/n_re/hottail");
    OptionConstants::uqty_f_hot_dist_mode ht_dist_mode = (enum OptionConstants::uqty_f_hot_dist_mode)s->GetInteger("eqsys/f_hot/dist_mode");    
    if(hottail_mode != OptionConstants::EQTERM_HOTTAIL_MODE_DISABLED && ht_dist_mode == OptionConstants::UQTY_F_HOT_DIST_MODE_NONREL){
        DEFU_FLD(TAU_COLL);
    }

}

/**
 * Implementation the EquationSystem class.
 */

#include <iostream>
#include <string>
#include <softlib/Timer.h>
#include "DREAM/EquationSystem.hpp"
#include "DREAM/IO.hpp"
#include "DREAM/QuitException.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Solver/SolverLinearlyImplicit.hpp"
#include "FVM/QuantityData.hpp"


using namespace DREAM;
using namespace std;

/**
 * Constructor.
 */
EquationSystem::EquationSystem(
    FVM::Grid *emptygrid, FVM::Grid *rgrid,
    enum OptionConstants::momentumgrid_type ht_type, FVM::Grid *hottailGrid,
    enum OptionConstants::momentumgrid_type re_type, FVM::Grid *runawayGrid,
    Settings *s
) : scalarGrid(emptygrid), fluidGrid(rgrid),
    hottailGrid(hottailGrid), runawayGrid(runawayGrid),
    hottailGrid_type(ht_type), runawayGrid_type(re_type),
    settings(s) {

    this->initializer = new EqsysInitializer(
        &this->unknowns, &this->unknown_equations,
        this, fluidGrid, hottailGrid, runawayGrid,
        ht_type, re_type
    );
}

/**
 * Destructor.
 */
EquationSystem::~EquationSystem() {
    if (this->ionHandler != nullptr)
        delete this->ionHandler;

    if (this->solver != nullptr)
        delete this->solver;
    
    if (this->timestepper != nullptr)
        delete this->timestepper;

    if (this->cqh_hottail != nullptr)
        delete this->cqh_hottail;

    if (this->cqh_runaway != nullptr)
        delete this->cqh_runaway;

    if (this->REFluid != nullptr)
        delete this->REFluid;

    if (this->distRE != nullptr)
        delete this->distRE;

    if (this->distHT != nullptr)
        delete this->distHT;
    
    if (this->postProcessor != nullptr)
        delete this->postProcessor;

    if (this->initializer != nullptr)
        delete this->initializer;

    if (this->settings != nullptr)
        delete this->settings;
	
	for (auto rsth : this->rsths)
		delete rsth;

    if (this->otherQuantityHandler != nullptr)
		delete this->otherQuantityHandler;

    if (this->SPI != nullptr)
		delete this->SPI;

	for (auto eqn : this->unknown_equations)
		delete eqn;
    
	FVM::RadialGrid *rgrid=nullptr;
    if (this->fluidGrid != nullptr){
        rgrid = this->fluidGrid->GetRadialGrid();
		delete this->fluidGrid;
    }

    if (this->hottailGrid != nullptr) {
		if (rgrid == nullptr)
			rgrid = this->hottailGrid->GetRadialGrid();
		delete this->hottailGrid;
	}

    if (this->runawayGrid != nullptr) {
		if (rgrid == nullptr)
			rgrid = this->runawayGrid->GetRadialGrid();
		delete this->runawayGrid;
	}

    if (this->scalarGrid != nullptr) {
		FVM::RadialGrid *rgs = this->scalarGrid->GetRadialGrid();
		delete this->scalarGrid;
		delete rgs;
	}

	if (rgrid != nullptr)
		delete rgrid;
}

/**
 * Processes the system after it has been initialized and
 * prepares it for being solved. This includes setting
 * initial values for those unknown quantities which do
 * not yet have initial values.
 *
 * t0: Time at which to set initial values.
 */
void EquationSystem::ProcessSystem(const real_t t0) {
    // Construct a list of the unknowns that will appear
    // in any matrices later on
    len_t totsize = 0;
    const len_t N = unknowns.GetNUnknowns();
    bool unknownMissing = false;
    for (len_t i = 0; i < N; i++) {
        if (unknown_equations[i] == nullptr) {
            DREAM::IO::PrintError("No equation has been declared for unknown '%s'", unknowns.GetUnknown(i)->GetName().c_str());
            unknownMissing = true;
        } else {
			if (unknown_equations[i]->IsSolvedExternally()) {
				external_unknowns.push_back(i);
            } else if (!unknown_equations[i]->IsPredetermined()) {
                nontrivial_unknowns.push_back(i);
                totsize += unknowns[i]->NumberOfElements();
            }
        }
    }

    // Initialize from output...
    if (this->initializerFile != "") 
        this->initializer->InitializeFromOutput(
            this->initializerFile, this->currentTime, this->initializerFileIndex,
            this->ionHandler, this->initializerFileIgnore
        );
	
	// Set external iterator
	this->extiter = new ExternalIterator(
		&this->unknowns, &this->unknown_equations
	);
	this->extiter->Initialize(this->external_unknowns);
    
    // Set initial values
    this->initializer->Execute(t0);

    if (unknownMissing)
        throw EquationSystemException("While processing equation system: Equations not declared for some unknowns.");

    this->matrix_size = totsize;
}

/**
 * Set one equation of the specified unknown.
 */
void EquationSystem::SetOperator(
	const len_t blockrow, const len_t blockcol, FVM::Operator *op,
	const std::string& desc, const bool solvedExternally
) {
    // Verify that the list is sufficiently large
    if (unknown_equations.size() < blockrow+1)
        unknown_equations.resize(unknowns.Size(), nullptr);

    // Does the unknown have any equations yet? If not, create
    // first the equation container
    if (unknown_equations[blockrow] == nullptr)
        unknown_equations[blockrow] = new UnknownQuantityEquation(blockrow, GetUnknown(blockrow), desc);

    unknown_equations[blockrow]->SetOperator(blockcol, op);

    if (desc != "") {
        unknown_equations[blockrow]->SetDescription(desc);
        unknown_equations[blockrow]->GetUnknown()->SetEquationDescription(desc);
    }

	if (solvedExternally)
		unknown_equations[blockrow]->SetExternallySolved(true);
}

/**
 * Same as 'SetEquation(len_t, len_t, Equation*)', but specifies
 * the unknowns by name rather than by index.
 */
void EquationSystem::SetOperator(len_t blockrow, const std::string& blockcol, FVM::Operator *op, const std::string& desc, const bool solvedExternally) {
    SetOperator(blockrow, GetUnknownID(blockcol), op, desc, solvedExternally);
}
void EquationSystem::SetOperator(const std::string& blockrow, len_t blockcol, FVM::Operator *op, const std::string& desc, const bool solvedExternally) {
    SetOperator(GetUnknownID(blockrow), blockcol, op, desc, solvedExternally);
}
void EquationSystem::SetOperator(const std::string& blockrow, const std::string& blockcol, FVM::Operator *op, const std::string& desc, const bool solvedExternally) {
    SetOperator(GetUnknownID(blockrow), GetUnknownID(blockcol), op, desc, solvedExternally);
}

/**
 * Set the initial value of the specified unknown quantity. If
 * the initial value has previously been specified, it is overwritten.
 *
 * id:  ID of unknown quantity.
 * val: Initial value of the quantity.
 * t0:  Initial time.
 */
void EquationSystem::SetInitialValue(const std::string& name, const real_t *val, const real_t t0) {
    this->SetInitialValue(this->unknowns.GetUnknownID(name), val, t0);
}
void EquationSystem::SetInitialValue(const len_t id, const real_t *val, const real_t t0) {
    this->unknowns.SetInitialValue(id, val, t0);
}

/**
 * Set the equation solver for this system.
 * Note that this method must be called only AFTER 'ProcessSystem()'
 * has been called.
 *
 * solver: Solver to assign
 */
void EquationSystem::SetSolver(Solver *solver) {
    if (this->matrix_size == 0)
        throw EquationSystemException("It appears that either 'EquationSystem::ProcessSystem()' has not been called prior to 'SetSolver()', or the matrix system is trivial (matrix_size = 0).");

    this->solver = solver;
    this->solver->Initialize(this->matrix_size, this->nontrivial_unknowns);

	if (this->extiter != nullptr)
		this->solver->SetExternalIterator(this->extiter);
}

/**
 * Solve this equation system.
 */
void EquationSystem::Solve() {
    this->currentTime = 0;
    this->times.push_back(this->currentTime);
    this->timestepper->SetSolver(solver);

    this->PrintNonTrivialUnknowns();
	this->PrintExternallyIteratedUnknowns();
    this->PrintTrivialUnknowns();

    // Set initial guess in solver
    const real_t *guess = unknowns.GetLongVector(this->nontrivial_unknowns);
    solver->SetInitialGuess(guess);
    delete [] guess;
    
    cout << "Beginning time advance..." << endl;

    Timer tim;
    len_t istep = 0;    // Number of times 'solver->Solve()' has been called...
    while (!timestepper->IsFinished()) {
        // Take step
        real_t tNext = timestepper->NextTime();
        this->currentTime = timestepper->CurrentTime();
        real_t dt = tNext - this->currentTime;

        this->fluidGrid->Rebuild(tNext);

        try {
            istep++;
            solver->Solve(tNext, dt);

            timestepper->ValidateStep();

            // Post-process solution (should be done before saving any
            // time step)
            this->postProcessor->Process(tNext);

            if (timestepper->IsSaveStep()) {
                this->TimestepFinished();

                // true = Really save the step (if it's false, we just
                // indicate that we have taken another timestep). This
                // should only be true for time steps which we want to
                // push to the output file.
                unknowns.SaveStep(tNext, true);
                this->times.push_back(tNext);

                otherQuantityHandler->StoreAll(tNext);
            } else
                unknowns.SaveStep(tNext, false);
            
            timestepper->PrintProgress();
        } catch (DREAM::QuitException& ex) {
            // Rethrow quit exception
            throw ex;
        } catch (FVM::FVMException& ex) {
            timestepper->HandleException(ex);
        }
    }

    cout << endl;

    this->simulationTime = tim.GetMicroseconds();
    string duration = tim.ToString();

    DREAM::IO::PrintInfo("Solved equation system in %s.", duration.c_str());

    if (this->timingStdout) {
        this->solver->PrintTimings();
        this->REFluid->PrintTimings();
    }
}

/**
 * Call all functions in 'callbacks_timestepFinished'.
 */
void EquationSystem::TimestepFinished() {
    for (auto f : this->callbacks_timestepFinished)
        (*f)(this->simulation);
}

/**
 * Register a function to call whenever a time step
 * has been taken.
 */
void EquationSystem::RegisterCallback_TimestepFinished(
    timestep_finished_func_t f
) {
    this->callbacks_timestepFinished.push_back(f);
}

/**
 * Register a function to call whenever a solver iteration
 * has been finished.
 */
void EquationSystem::RegisterCallback_IterationFinished(
    Solver::iteration_finished_func_t f
) {
    this->solver->RegisterCallback_IterationFinished(f);
}

/**
 * Returns the maximum simulation time for this simulation.
 */
real_t EquationSystem::GetMaxTime() const {
    return this->timestepper->MaxTime();
}

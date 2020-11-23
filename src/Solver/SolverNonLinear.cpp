/**
 * Implementation of a custom Newton solver which only utilizes
 * the linear solvers of PETSc.
 */

#include <string>
#include <vector>
#include "DREAM/IO.hpp"
#include "DREAM/OutputGeneratorSFile.hpp"
#include "DREAM/Solver/SolverNonLinear.hpp"
#include "FVM/Solvers/MILU.hpp"
#include "FVM/Solvers/MIKSP.hpp"
#include "FVM/Solvers/MIMUMPS.hpp"


using namespace DREAM;
using namespace std;


/**
 * Constructor.
 */
SolverNonLinear::SolverNonLinear(
	FVM::UnknownQuantityHandler *unknowns,
	vector<UnknownQuantityEquation*> *unknown_equations,
    EquationSystem *eqsys,
    enum OptionConstants::linear_solver ls,
	const int_t maxiter, const real_t reltol,
	bool verbose
) : Solver(unknowns, unknown_equations), eqsys(eqsys), linearSolver(ls),
	maxiter(maxiter), reltol(reltol), verbose(verbose) {

    this->timeKeeper = new FVM::TimeKeeper("Solver non-linear");
    this->timerTot = this->timeKeeper->AddTimer("total", "Total time");
    this->timerRebuild = this->timeKeeper->AddTimer("rebuildtot", "Rebuild coefficients");
    this->timerResidual = this->timeKeeper->AddTimer("residual", "Construct residual");
    this->timerJacobian = this->timeKeeper->AddTimer("jacobian", "Construct jacobian");
    this->timerInvert = this->timeKeeper->AddTimer("invert", "Invert jacobian");
}

/**
 * Destructor.
 */
SolverNonLinear::~SolverNonLinear() {
	Deallocate();

    delete this->timeKeeper;
}


/**
 * "Accept" the current solution and prepare for taking
 * another Newton step.
 */
void SolverNonLinear::AcceptSolution() {
	real_t *x = this->x1;
	this->x1  = this->x0;
	this->x0  = x;

	this->StoreSolution(x);
}

/**
 * Allocate memory for all objects used by this solver.
 */
void SolverNonLinear::Allocate() {
	jacobian = new FVM::BlockMatrix();

	for (len_t i = 0; i < nontrivial_unknowns.size(); i++) {
		len_t id = nontrivial_unknowns[i];
		UnknownQuantityEquation *eqn = this->unknown_equations->at(id);

		unknownToMatrixMapping[id] =
			jacobian->CreateSubEquation(eqn->NumberOfElements(), eqn->NumberOfNonZeros_jac(), id);
	}

	jacobian->ConstructSystem();

	const len_t N = jacobian->GetNRows();

    // Select linear solver
    if (this->linearSolver == OptionConstants::LINEAR_SOLVER_LU)
        this->inverter = new FVM::MILU(N);
    else if (this->linearSolver == OptionConstants::LINEAR_SOLVER_MUMPS)
        this->inverter = new FVM::MIMUMPS(N);
    else
        throw SolverException(
            "Unrecognized linear solver specified: %d.", this->linearSolver
        );

	VecCreateSeq(PETSC_COMM_WORLD, N, &this->petsc_F);
	VecCreateSeq(PETSC_COMM_WORLD, N, &this->petsc_dx);

	this->x0 = new real_t[N];
	this->x1 = new real_t[N];
	this->dx = new real_t[N];

	this->x_2norm  = new real_t[this->unknown_equations->size()];
	this->dx_2norm = new real_t[this->unknown_equations->size()];
}

/**
 * Deallocate memory used by this solver.
 */
void SolverNonLinear::Deallocate() {
	delete inverter;
	delete jacobian;

	delete [] this->x_2norm;
	delete [] this->dx_2norm;

	delete [] this->x0;
	delete [] this->x1;
	delete [] this->dx;

	VecDestroy(&this->petsc_F);
	VecDestroy(&this->petsc_dx);
}

/**
 * Returns the name of the specified non-trivial unknown quantity.
 *
 * idx: Index into 'this->nontrivial_unknowns' of the non-trivial unknown
 *      to return the name of.
 */
const string& SolverNonLinear::GetNonTrivialName(const len_t idx) {
    return this->unknowns->GetUnknown(this->nontrivial_unknowns[idx])->GetName();
}

/**
 * Initialize the solver.
 */
void SolverNonLinear::initialize_internal(
	const len_t, vector<len_t>&
) {
	this->Allocate();

    if (this->convChecker == nullptr)
        this->SetConvergenceChecker(
            new ConvergenceChecker(unknowns, this->nontrivial_unknowns, this->reltol)
        );
}

/**
 * Check if the solver has converged.
 */
bool SolverNonLinear::IsConverged(const real_t *x, const real_t *dx) {
	if (this->GetIteration() >= (len_t)this->MaxIter()){
		throw SolverException(
			"Non-linear solver reached the maximum number of allowed "
			"iterations: " LEN_T_PRINTF_FMT ".",
			this->MaxIter()
		);
	}

    if (this->Verbose())
        DREAM::IO::PrintInfo("ITERATION %d", this->GetIteration());

    return convChecker->IsConverged(x, dx, this->Verbose());
}

/**
 * Set the initial guess for the solver.
 *
 * guess: Vector containing values of initial guess.
 */
void SolverNonLinear::SetInitialGuess(const real_t *guess) {
	if (guess != nullptr) {
		for (len_t i = 0; i < this->matrix_size; i++)
			this->x0[i] = guess[i];
	} else {
		for (len_t i = 0; i < this->matrix_size; i++)
			this->x0[i] = 0;
	}
}

/**
 * Solve the equation system (advance the system in time
 * by one step).
 *
 * t:  Time at which the current solution is given.
 * dt: Time step to take.
 *
 * (the obtained solution will correspond to time t'=t+dt)
 */
void SolverNonLinear::Solve(const real_t t, const real_t dt) {
    this->nTimeStep++;

	this->t  = t;
	this->dt = dt;

    this->timeKeeper->StartTimer(timerTot);

	// Take Newton steps
	len_t iter = 0;
	const real_t *x, *dx;
	do {
		iter++;
		this->SetIteration(iter);

		dx = this->TakeNewtonStep();
		x  = UpdateSolution(dx);

		// TODO backtracking...
		
		AcceptSolution();
	} while (!IsConverged(x, dx));
	

    this->timeKeeper->StopTimer(timerTot);
}

/**
 * Debugging routine for saving both the "analytically" computed
 * Jacobian, as well as the Jacobian evaluated numerically using
 * finite differences, to file. When this method is called, the
 * 'jacobian' variable is assumed to contain the "analytical"
 * Jacobian matrix for the current time/iteration. This routine
 * will then save that matrix, compute the corresponding numerical
 * Jacobian, and save that.
 *
 * name: Base name to use for files.
 */
void SolverNonLinear::SaveNumericalJacobian(const std::string& name) {
    this->_EvaluateJacobianNumerically(this->jacobian);
    this->jacobian->View(FVM::Matrix::BINARY_MATLAB, name + "_num");
}

void SolverNonLinear::SaveJacobian() {
    this->jacobian->View(FVM::Matrix::BINARY_MATLAB, "petsc_jacobian");
}
void SolverNonLinear::SaveJacobian(const std::string& name) {
    this->jacobian->View(FVM::Matrix::BINARY_MATLAB, name);
}

/**
 * Store the current solution to the UnknownQuantityHandler.
 */
void SolverNonLinear::StoreSolution(const real_t *x) {
	this->unknowns->Store(this->nontrivial_unknowns, x);
}

/**
 * Calculate the next Newton step to take.
 */
const real_t *SolverNonLinear::TakeNewtonStep() {
    this->timeKeeper->StartTimer(timerRebuild);
	this->RebuildTerms(this->t, this->dt);
    this->timeKeeper->StopTimer(timerRebuild);

	// Evaluate function vector
    this->timeKeeper->StartTimer(timerResidual);
	real_t *fvec;
	VecGetArray(this->petsc_F, &fvec);
	this->BuildVector(this->t, this->dt, fvec, this->jacobian);
	VecRestoreArray(this->petsc_F, &fvec);
    this->timeKeeper->StopTimer(timerResidual);


	// Evaluate jacobian
    this->timeKeeper->StartTimer(timerJacobian);
	this->BuildJacobian(this->t, this->dt, this->jacobian);
    this->timeKeeper->StopTimer(timerJacobian);

    // Print/save debug info (if requested)
    this->SaveDebugInfo(this->nTimeStep, this->iteration);

    // Apply preconditioner (if enabled)
    this->Precondition(this->jacobian, this->petsc_F);

	// Solve J*dx = F
    this->timeKeeper->StartTimer(timerInvert);
	inverter->Invert(this->jacobian, &this->petsc_F, &this->petsc_dx);
    this->timeKeeper->StopTimer(timerInvert);

    // Undo preconditioner (if enabled)
    this->UnPrecondition(this->petsc_dx);

	// Copy dx
	VecGetArray(this->petsc_dx, &fvec);
	for (len_t i = 0; i < this->matrix_size; i++)
		this->dx[i] = fvec[i];
	VecRestoreArray(this->petsc_dx, &fvec);
	
	return this->dx;
}



/**
 * Helper function to damping factor calculation; given a quantity
 * X0 and its change dX in the iteration, returns the maximal step 
 * length (damping) such that X>threshold*X0 after the iteration
 */
const real_t MaximalStepLengthAtGridPoint(
	real_t X0, real_t dX, real_t threshold
){
	real_t maxStepAtI = 1.0;
	if(dX)
		maxStepAtI = (1-threshold) * X0 / dX;
	return maxStepAtI;
}

/**
 * Returns a dampingFactor such that x1 = x0 - dampingFactor*dx satisfies 
 * physically-motivated constraints, such as positivity of temperature.
 * If initial guess dx from Newton step satisfies all constraints, returns 1.
 */
const real_t MaximalPhysicalStepLength(real_t *x0, const real_t *dx,len_t iteration, std::vector<len_t> nontrivial_unknowns, FVM::UnknownQuantityHandler *unknowns, IonHandler *ionHandler, len_t &id_uqn){
	real_t maxStepLength = 1.0;
	real_t threshold = 0.1;

	std::vector<len_t> ids_nonNegativeQuantities;
	// add those quantities which we expect to be non-negative
	// T_cold and n_cold will crash the simulation if negative, so they should always be added
	ids_nonNegativeQuantities.push_back(unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD));
	ids_nonNegativeQuantities.push_back(unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD));

	bool nonNegativeZeff = true;
	const len_t id_ni = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);

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
					id_uqn = id;
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
					id_uqn = id;
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

/**
 * Update the current solution with the Newton step 'dx'.
 *
 * dx: Newton step to take.
 */
const real_t *SolverNonLinear::UpdateSolution(const real_t *dx) {
    len_t id_uqn;
	real_t dampingFactor = MaximalPhysicalStepLength(x0,dx,iteration,nontrivial_unknowns,unknowns,ionHandler,id_uqn);
	
	if(dampingFactor < 1 && this->Verbose()) {
        DREAM::IO::PrintInfo();
		DREAM::IO::PrintInfo("Newton iteration dynamically damped for unknown quantity: %s",unknowns->GetUnknown(id_uqn)->GetName().c_str());
		DREAM::IO::PrintInfo("to conserve positivity, by a factor: %e", dampingFactor);
        DREAM::IO::PrintInfo();
	}
	for (len_t i = 0; i < this->matrix_size; i++)
		this->x1[i] = this->x0[i] - dampingFactor*dx[i];
	
	return this->x1;
}

/**
 * Print timing information after the solve.
 */
void SolverNonLinear::PrintTimings() {
    this->timeKeeper->PrintTimings(true, 0);
    this->Solver::PrintTimings_rebuild();
}

/**
 * Save timing information to the given SFile object.
 *
 * sf:   SFile object to save timing information to.
 * path: Path in file to save timing information to.
 */
void SolverNonLinear::SaveTimings(SFile *sf, const string& path) {
    this->timeKeeper->SaveTimings(sf, path);

    sf->CreateStruct(path+"/rebuild");
    this->Solver::SaveTimings_rebuild(sf, path+"/rebuild");
}

/**
 * Save debugging information for the current iteration.
 *
 * iTimeStep:  Current time step index.
 * iIteration: Current iteration index.
 */
void SolverNonLinear::SaveDebugInfo(
    len_t iTimeStep, len_t iIteration
) {
    if ((this->savetimestep == iTimeStep &&
        (this->saveiteration == iIteration || this->saveiteration == 0)) ||
         this->savetimestep == 0) {

        string suffix = "_" + to_string(iTimeStep) + "_" + to_string(iIteration);
        
        if (this->savejacobian) {
            string jacname;
            if (this->savetimestep == 0 || this->saveiteration == 0)
                jacname = "petsc_jac" + suffix;
            else
                jacname = "petsc_jac";

            SaveJacobian(jacname);
        }

        if (this->savevector) {
            string resname;
            if (this->savetimestep == 0 || this->saveiteration == 0)
                resname = "residual" + suffix + ".mat";
            else
                resname = "residual.mat";

            real_t *fvec;
            VecGetArray(this->petsc_F, &fvec);

            SFile *sf = SFile::Create(resname, SFILE_MODE_WRITE);
            sf->WriteList("F", fvec, this->jacobian->GetNRows());
            sf->Close();

            VecRestoreArray(this->petsc_F, &fvec);
        }

        if (this->savenumjac) {
            string jacname;
            if (this->savetimestep == 0 || this->saveiteration == 0)
                jacname = "petsc_jac" + suffix;
            else
                jacname = "petsc_jac";

            SaveNumericalJacobian(jacname);
        }

        // Save full output?
        if (this->savesystem) {
            string outname = "debugout";
            if (this->savetimestep == 0 || this->saveiteration == 0)
                outname += suffix;
            outname += ".h5";

            OutputGeneratorSFile *outgen = new OutputGeneratorSFile(this->eqsys, outname);
            outgen->SaveCurrent();
            delete outgen;
        }

        if (this->printjacobianinfo)
            this->jacobian->PrintInfo();
    }
}

/**
 * Enable or disable debug mode.
 *
 * printjacobianinfo: If true, prints detailed debug information about the
 *                    PETSc matrix for the jacobian in every iteration.
 * savejacobian:      If true, saves the jacobian using a PETSc viewer in
 *                    the specified time step(s).
 * savevector:        If true, saves the residual vector in the specified
 *                    time step(s).
 * savenumjac:        If true, calculates the jacobian numerically and saves
 *                    it using a PETSc viewer in the specified time step(s).
 * timestep:          Time step index to save debug info for. If 0, saves
 *                    the information in every iteration of every time step.
 * iteration:         Iteration of specified time step to save debug info for.
 * savesystem:        If true, saves the full equation system, including grid information,
 *                    to a proper DREAMOutput file. However, only the most recently obtained
 *                    solution is saved.
 */
void SolverNonLinear::SetDebugMode(
    bool printjacobianinfo, bool savejacobian, bool savevector,
    bool savenumjac, int_t timestep, int_t iteration, bool savesystem
) {
    this->printjacobianinfo = printjacobianinfo;
    this->savejacobian      = savejacobian;
    this->savevector        = savevector;
    this->savenumjac        = savenumjac;
    this->savetimestep      = timestep;
    this->saveiteration     = iteration;
    this->savesystem        = savesystem;
}


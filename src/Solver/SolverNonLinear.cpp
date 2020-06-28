/**
 * Implementation of a custom Newton solver which only utilizes
 * the linear solvers of PETSc.
 */

#include <string>
#include <vector>
#include "DREAM/IO.hpp"
#include "DREAM/Solver/SolverNonLinear.hpp"
#include "FVM/Solvers/MILU.hpp"
#include "FVM/Solvers/MIKSP.hpp"


using namespace DREAM;
using namespace std;


/**
 * Constructor.
 */
SolverNonLinear::SolverNonLinear(
	FVM::UnknownQuantityHandler *unknowns,
	vector<UnknownQuantityEquation*> *unknown_equations,
    enum OptionConstants::linear_solver ls,
	const int_t maxiter, const real_t reltol,
	bool verbose
) : Solver(unknowns, unknown_equations), linearSolver(ls),
	maxiter(maxiter), reltol(reltol), verbose(verbose) {
}

/**
 * Destructor.
 */
SolverNonLinear::~SolverNonLinear() {
	Deallocate();
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
			jacobian->CreateSubEquation(eqn->NumberOfElements(), eqn->NumberOfNonZeros_jac());
	}

	jacobian->ConstructSystem();

	const len_t N = jacobian->GetNRows();

    // Select linear solver
    if (this->linearSolver == OptionConstants::LINEAR_SOLVER_LU)
        this->inverter = new FVM::MILU(N);
    else if (this->linearSolver == OptionConstants::LINEAR_SOLVER_GMRES)
        this->inverter = new FVM::MIKSP(N);
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
}

/**
 * Check if the solver has converged.
 */
bool SolverNonLinear::IsConverged(const real_t *x, const real_t *dx) {
	if (this->GetIteration() >= (len_t)this->MaxIter())
		throw SolverException(
			"Non-linear solver reached the maximum number of allowed "
			"iterations: " LEN_T_PRINTF_FMT ".",
			this->MaxIter()
		);

    this->CalculateNonTrivial2Norm(x, this->x_2norm);
    this->CalculateNonTrivial2Norm(dx, this->dx_2norm);

    // Iterate over norms and ensure that all are small
    const len_t N = this->nontrivial_unknowns.size();
    bool converged = true;

    if (this->Verbose())
        DREAM::IO::PrintInfo("ITERATION %d", this->GetIteration());

    for (len_t i = 0; i < N; i++) {
        bool conv = (dx_2norm[i] < this->reltol*x_2norm[i]) && (x_2norm[i]>0);

        if (this->Verbose()) {
#ifdef COLOR_TERMINAL
            if (conv)
                DREAM::IO::PrintInfo(
                    "   \x1B[32m%10s  |x| = %e, |dx| = %e\x1B[0m",
                    this->GetNonTrivialName(i).c_str(),
                    x_2norm[i], dx_2norm[i]
                );
             else
                DREAM::IO::PrintInfo(
                    "   \x1B[1;31m%10s  |x| = %e, |dx| = %e\x1B[0m",
                    this->GetNonTrivialName(i).c_str(),
                    x_2norm[i], dx_2norm[i]
                );
#else
            DREAM::IO::PrintInfo(
                "   %10s  |x| = %e, |dx| = %e",
                this->GetNonTrivialName(i).c_str(),
                x_2norm[i], dx_2norm[i]
            );
#endif
        }

        converged = converged && conv;
    }

//	this->jacobian->PrintInfo();

	return converged;
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
	this->t  = t;
	this->dt = dt;

	// Take Newton steps
	len_t iter = 0;
	const real_t *x, *dx;
	do {
		iter++;
		this->SetIteration(iter);

		dx = this->TakeNewtonStep();
		x  = UpdateSolution(dx);


/*
		if (iter==1) {
			SaveJacobian();
//            SaveNumericalJacobian();
            throw SolverException("Stopping now.");
        }
*/		
		// TODO backtracking...
		
		AcceptSolution();

	} while (!IsConverged(x, dx));
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
	this->RebuildTerms(this->t, this->dt);

	// Evaluate function vector
	real_t *fvec;
	VecGetArray(this->petsc_F, &fvec);
	this->BuildVector(this->t, this->dt, fvec, this->jacobian);
	VecRestoreArray(this->petsc_F, &fvec);


	// Evaluate jacobian
	this->BuildJacobian(this->t, this->dt, this->jacobian);

	

	// Solve J*dx = F
	inverter->Invert(this->jacobian, &this->petsc_F, &this->petsc_dx);

	// Copy dx
	VecGetArray(this->petsc_dx, &fvec);
	for (len_t i = 0; i < this->matrix_size; i++)
		this->dx[i] = fvec[i];
	VecRestoreArray(this->petsc_dx, &fvec);
	
	return this->dx;
}


/**
 * Returns a dampingFactor such that x1 = x0 - dampingFactor*dx satisfies 
 * physically-motivated constraints, such as positivity of temperature.
 * If initial guess dx from Newton step satisfies all constraints, returns 1.
 */
const real_t MaximalPhysicalStepLength(real_t *x0, const real_t *dx, std::vector<len_t> nontrivial_unknowns, FVM::UnknownQuantityHandler *unknowns ){
	real_t maxStepLength = 1;
	real_t threshold = 0.1;

	std::vector<len_t> ids_nonNegativeQuantities;
	// add those quantities which we expect to be non-negative
	ids_nonNegativeQuantities.push_back(unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD));
	ids_nonNegativeQuantities.push_back(unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD));
	//ids_nonNegativeQuantities.push_back(unknowns->GetUnknownID(OptionConstants::UQTY_N_RE));

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
		if(isNonNegativeQuantity){
			for(len_t i=0; i<NCells; i++){
				// require x1 > threshold*x0
				real_t maxStepAtI = (1-threshold) * x0[offset + i] / abs(dx[offset + i]);
				// if this is a stronger constaint than current maxlength, override
				if(maxStepAtI < maxStepLength)
					maxStepLength = maxStepAtI;
			}
		}
		offset += NCells;
	}
	return maxStepLength;
}

/**
 * Update the current solution with the Newton step 'dx'.
 *
 * dx: Newton step to take.
 */
#include <iostream>
const real_t *SolverNonLinear::UpdateSolution(const real_t *dx) {

	real_t dampingFactor = MaximalPhysicalStepLength(x0,dx,nontrivial_unknowns,unknowns);
	if(dampingFactor < 1){
		std::cout << std::endl;
		std::cout << "Newton iteration dynamically damped" << std::endl;
		std::cout << "to conserve positivity, by a factor: " << dampingFactor << std::endl;
		std::cout << std::endl;
	}
	for (len_t i = 0; i < this->matrix_size; i++)
		this->x1[i] = this->x0[i] - dampingFactor*dx[i];
	
	return this->x1;
}


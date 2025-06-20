/**
 * Definition of equations relating to the poloidal flux.
 * In DREAM, the poloidal flux is normalised to the major
 * radius of the magnetic axis,
 * psi_p = poloidal flux / R0
 */

#include <iostream>
#include <string>
#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Equation/Operator.hpp"
#include "FVM/Interpolator3D.hpp"
#include "FVM/Equation/PrescribedParameter.hpp"
#include "FVM/Equation/ConstantParameter.hpp"
#include "FVM/Equation/DiagonalLinearTerm.hpp"
#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Equation/TransientTerm.hpp"
#include "DREAM/Equations/Fluid/AmperesLawBoundaryAtRMax.hpp"
#include "DREAM/Equations/Scalar/WallCurrentTerms.hpp"

using namespace DREAM;
using namespace std;

// Implementation of equation terms in Ampere's law
namespace DREAM {
    // Term representing the j_||/(B/Bmin) term
    class AmperesLawJTotTerm : public FVM::DiagonalLinearTerm {
    public:
        AmperesLawJTotTerm(FVM::Grid* g) : FVM::DiagonalLinearTerm(g){}

        virtual void SetWeights() override {
            len_t offset = 0;
            for (len_t ir = 0; ir < nr; ir++){
                real_t w = 2*M_PI*Constants::mu0 * grid->GetRadialGrid()->GetFSA_1OverR2(ir) 
                    * grid->GetRadialGrid()->GetBTorG(ir) / grid->GetRadialGrid()->GetBmin(ir);
                for(len_t i = 0; i < n1[ir]*n2[ir]; i++)
                    weights[offset + i] = w;
                offset += n1[ir]*n2[ir];
            }
        }
    };
    // Term representing the diffusion term on psi_p
    class AmperesLawDiffusionTerm : public FVM::DiffusionTerm {
    public:
        AmperesLawDiffusionTerm(FVM::Grid *g ) : FVM::DiffusionTerm(g) {}        
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override {
            for (len_t ir = 0; ir <= nr; ir++) {
                real_t drr = grid->GetRadialGrid()->GetFSA_NablaR2OverR2_f(ir);
                for (len_t j = 0; j < n2[0]; j++) 
                    for (len_t i = 0; i < n1[0]; i++) 
                        Drr(ir, i, j) += drr;
            }
        }
    };
}

#define MODULENAME "eqsys/E_field/bc"

/**
 * Construct the equation for the poloidal flux j_|| ~ Laplace(psi)
 * as well as impose boundary conditions from a 0D circuit equation
 * including a resistive wall.
 * 
 * eqsys: Equation system to put the equation in.
 * s:     Settings object describing how to construct the equation.
 */
void SimulationGenerator::ConstructEquation_psi_p(
    EquationSystem *eqsys, Settings*
) {
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    FVM::Grid *scalarGrid = eqsys->GetScalarGrid();
    
    const len_t id_I_p = eqsys->GetUnknownID(OptionConstants::UQTY_I_P);    
    const len_t id_psi_p = eqsys->GetUnknownID(OptionConstants::UQTY_POL_FLUX);
    const len_t id_psi_edge = eqsys->GetUnknownID(OptionConstants::UQTY_PSI_EDGE);    
    const len_t id_j_tot = eqsys->GetUnknownID(OptionConstants::UQTY_J_TOT);

    // Set equation j_tot ~ d_r^2(psi_p)
    FVM::Operator *eqn_j1 = new FVM::Operator(fluidGrid);
    FVM::Operator *eqn_j2 = new FVM::Operator(fluidGrid);
    FVM::Operator *eqn_j3 = new FVM::Operator(fluidGrid);

    eqn_j1->AddTerm(new AmperesLawJTotTerm(fluidGrid));
    eqn_j2->AddTerm(new AmperesLawDiffusionTerm(fluidGrid));
    eqsys->SetOperator(id_psi_p, id_j_tot, eqn_j1, "Poloidal flux Ampere's law");

	/*enum OptionConstants::uqty_E_field_eqn Etype =
		(enum OptionConstants::uqty_E_field_eqn)s->GetInteger("eqsys/E_field/type");*/

	ConstructEquation_psi_init_nl(
		eqsys, id_psi_p, id_j_tot, id_I_p, id_psi_edge
	);
	/*ConstructEquation_psi_init_integral(
		eqsys, fluidGrid, id_psi_p, id_j_tot, id_I_p
	);*/

    // Set outgoing flux from diffusion term due to dpsi/dr at r=a,
    // obtained from psi_edge = psi(a)
    eqn_j2->AddBoundaryCondition(new FVM::BC::AmperesLawBoundaryAtRMax(fluidGrid,fluidGrid,eqn_j2,-1.0));
    eqn_j3->AddBoundaryCondition(new FVM::BC::AmperesLawBoundaryAtRMax(fluidGrid,scalarGrid,eqn_j2,+1.0));
    eqsys->SetOperator(id_psi_p, id_psi_edge, eqn_j3);
    eqsys->SetOperator(id_psi_p, id_psi_p, eqn_j2);
    
}

/**
 * Construct initialization rule for the poloidal flux which solves
 * for the psi_p directly from Ampere's law.
 */
void SimulationGenerator::ConstructEquation_psi_init_nl(
	EquationSystem *eqsys,
	const len_t id_psi_p, const len_t id_j_tot, const len_t id_I_p,
	const len_t id_psi_edge
) {
    eqsys->initializer->AddRule(
        id_psi_p,
        EqsysInitializer::INITRULE_STEADY_STATE_SOLVE,
		nullptr,
        // Dependencies
        id_j_tot,
		id_psi_edge,
        id_I_p,
		EqsysInitializer::RUNAWAY_FLUID
    );
}

/**
 * Construct initialization rule for the poloidal flux, utilizing
 * the integral relationship between the poloidal flux and total
 * current density.
 */
void SimulationGenerator::ConstructEquation_psi_init_integral(
	EquationSystem *eqsys, Settings *s, FVM::Grid *fluidGrid,
	const len_t id_psi_p, const len_t id_j_tot, const len_t id_I_p
) {
    /**
     * Initialization: define the function which integrates j_tot.
     * In principle, this expression is a direct inversion of the 
     * equation given previously with boundary condition psi(r_max) = 0 
     * It would be nicer to just solve the equation. (may be tricky since 
     * the boundary condittion may differ from the one we wish to use later)
     */
    FVM::RadialGrid *rGrid = fluidGrid->GetRadialGrid();
    real_t a = fluidGrid->GetRadialGrid()->GetMinorRadius();
    real_t b = (real_t)s->GetReal("radialgrid/wall_radius");
    real_t M_inductance = PlasmaEdgeToWallInductanceTerm::GetInductance(a,b);
    std::function<void(FVM::UnknownQuantityHandler*, real_t*)> initfunc_PsiPFromJtot 
        = [rGrid,M_inductance](FVM::UnknownQuantityHandler*u, real_t *psi_p_init)
    {
        len_t id_j_tot = u->GetUnknownID(OptionConstants::UQTY_J_TOT);
        len_t id_I_p = u->GetUnknownID(OptionConstants::UQTY_I_P);
        
        len_t nr = rGrid->GetNr();
        real_t *Itot = new real_t[nr];

        real_t *j_tot_init = u->GetUnknownData(id_j_tot);
        real_t *I_p_init = u->GetUnknownData(id_I_p);

        Itot[0] = TotalPlasmaCurrentFromJTot::GetIpIntegrand(0,rGrid) * j_tot_init[0];
        for(len_t ir=1; ir<nr; ir++)
            Itot[ir] = Itot[ir-1] + TotalPlasmaCurrentFromJTot::GetIpIntegrand(ir,rGrid) * j_tot_init[ir];

        // we use the convention that the initial poloidal flux at the wall is 0
        real_t psi_edge_init = -M_inductance*Itot[nr-1]; 

        const real_t *r = rGrid->GetR();
        const real_t *dr = rGrid->GetDr();
        const real_t a = rGrid->GetR_f(nr);
        #define integrand(I, Ip) 2*M_PI*Constants::mu0*Ip/(rGrid->GetVpVol(I)*rGrid->GetFSA_NablaR2OverR2_f(I))
        psi_p_init[nr-1] = psi_edge_init - (a-r[nr-1])*integrand(nr-1, I_p_init[0]);
        if(nr>1)
            for(len_t ir = nr-2; true; ir--){
                psi_p_init[ir] = psi_p_init[ir+1] - dr[ir]*integrand(ir, Itot[ir]);
                if(ir==0)
                    break;
            }

        #undef integrand
        delete [] Itot;
    };
            
    eqsys->initializer->AddRule(
        id_psi_p,
        EqsysInitializer::INITRULE_EVAL_FUNCTION,
        initfunc_PsiPFromJtot,
        // Dependencies
        id_j_tot,
        id_I_p
    );
}


/**
 * Sets the equation for psi_edge = psi(r=a)
 */
void SimulationGenerator::ConstructEquation_psi_edge(
    EquationSystem *eqsys, Settings *s
) {
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    FVM::Grid *scalarGrid = eqsys->GetScalarGrid();
    FVM::UnknownQuantityHandler *unknowns = eqsys->GetUnknownHandler();

    const len_t id_psi_edge    = unknowns->GetUnknownID(OptionConstants::UQTY_PSI_EDGE);
    const len_t id_psi_wall    = unknowns->GetUnknownID(OptionConstants::UQTY_PSI_WALL);
    const len_t id_I_p         = unknowns->GetUnknownID(OptionConstants::UQTY_I_P);
    
    real_t a = fluidGrid->GetRadialGrid()->GetMinorRadius();
    real_t b = (real_t)s->GetReal("radialgrid/wall_radius");

    // Set equation "psi_edge = psi_w - I_p*M"
    FVM::Operator *Op_psi_edge_1 = new FVM::Operator(scalarGrid);
    FVM::Operator *Op_psi_edge_2 = new FVM::Operator(scalarGrid);
    Op_psi_edge_1->AddTerm(new FVM::IdentityTerm(scalarGrid,-1.0));
    Op_psi_edge_2->AddTerm(new FVM::IdentityTerm(scalarGrid));
    eqsys->SetOperator(id_psi_edge, id_psi_wall, Op_psi_edge_2);      
    std::string desc = "psi_edge = psi_wall";  
    if(b>a){
        FVM::Operator *Op_psi_edge_3 = new FVM::Operator(scalarGrid);
        Op_psi_edge_3->AddTerm(new PlasmaEdgeToWallInductanceTerm(scalarGrid,a,b));
        eqsys->SetOperator(id_psi_edge, id_I_p, Op_psi_edge_3);
        desc += " - M_pw*I_p";
    }
    eqsys->SetOperator(id_psi_edge, id_psi_edge, Op_psi_edge_1,desc);

    eqsys->initializer->AddRule(
        id_psi_edge,
        EqsysInitializer::INITRULE_EVAL_EQUATION,
        nullptr,
        id_psi_wall,
        id_I_p
    );
}


/*
 * Sets the equation for the poloidal flux at the wall psi_w = psi(r=b),
 * the loop voltage at the wall V_loop_wall and the wall current I_w (if applicable).
 * Either the loop voltage is prescribed, or it is modelled with a circuit model,
 * documented in doc/notes/theory
 */
void SimulationGenerator::ConstructEquation_psi_wall_selfconsistent(
    EquationSystem *eqsys, Settings *s
) {
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    FVM::Grid *scalarGrid = eqsys->GetScalarGrid();
    FVM::UnknownQuantityHandler *unknowns = eqsys->GetUnknownHandler();

    eqsys->SetUnknown(OptionConstants::UQTY_V_LOOP_WALL, OptionConstants::UQTY_V_LOOP_WALL_DESC, scalarGrid);

    const len_t id_V_loop_wall = unknowns->GetUnknownID(OptionConstants::UQTY_V_LOOP_WALL);
    const len_t id_psi_wall    = unknowns->GetUnknownID(OptionConstants::UQTY_PSI_WALL);
    const len_t id_I_p         = unknowns->GetUnknownID(OptionConstants::UQTY_I_P);

    // Set V_loop_wall equation
    enum OptionConstants::uqty_V_loop_wall_eqn type = (enum OptionConstants::uqty_V_loop_wall_eqn)s->GetInteger(MODULENAME "/type");
    if(type == OptionConstants::UQTY_V_LOOP_WALL_EQN_PRESCRIBED){
        // Set V_loop_wall to prescribed
        FVM::Operator *Op_psi_wall_1 = new FVM::Operator(scalarGrid);
        FVM::Interpolator1D *interp = LoadDataT(MODULENAME, s, "V_loop_wall");
        Op_psi_wall_1->AddTerm(new FVM::PrescribedParameter(scalarGrid, interp));
        eqsys->SetOperator(id_V_loop_wall, OptionConstants::UQTY_V_LOOP_WALL, Op_psi_wall_1, "Prescribed");

    } else if (
        type == OptionConstants::UQTY_V_LOOP_WALL_EQN_SELFCONSISTENT ||
        type == OptionConstants::UQTY_V_LOOP_WALL_EQN_TRANSFORMER
    ) {
        // Inverse wall time in 1/s
        real_t wall_freq = (real_t)s->GetReal(MODULENAME "/inverse_wall_time");
        if(wall_freq == 0){
            // Prescribe V_loop_wall = 0.
            // Same as type PRESCRIBED with V=0.
            FVM::Operator *Op_V_loop_wall_1 = new FVM::Operator(scalarGrid);
            
            Op_V_loop_wall_1->AddTerm(new FVM::ConstantParameter(scalarGrid,0.0));
            eqsys->SetOperator(id_V_loop_wall, id_V_loop_wall, Op_V_loop_wall_1, "zero");
        } else {
            // Introduce I_w and set 
            //      V_loop_wall = R_W * I_w
            //      dpsi_w/dt = -L_w*(dI_p/dt+dI_w/dt)
            eqsys->SetUnknown(OptionConstants::UQTY_I_WALL, OptionConstants::UQTY_I_WALL_DESC, scalarGrid);
            const len_t id_I_w = unknowns->GetUnknownID(OptionConstants::UQTY_I_WALL);

            real_t R0 = (real_t)s->GetReal(MODULENAME "/R0");
            if(R0==0)
                R0 = fluidGrid->GetRadialGrid()->GetR0();
            if(isinf(R0))
                throw FVM::FVMException("Invalid major radius: Cannot be inf (cylindrical plasma) "
                                      "with finite wall time due to divergent external inductance.");

            real_t b = (real_t)s->GetReal("radialgrid/wall_radius");
            // External inductance normalized to R0
            real_t L_ext = Constants::mu0 * log(R0/b);
            // Wall resistivity
            real_t R_W = L_ext * wall_freq;

            /** 
             * Typically L_ext ~ mu0*log(R0/b), and the wall resistivity R_W is such that the wall
             * time tau_wall = L_W/R_W takes a value which is known experimentally/from simulation
             * for various devices. For example, it has been estimated to ~10 ms in DIII-D and in
             * JET, ~1 ms in FTU and ~500 ms in ITER. 
             */
            FVM::Operator *Op_V_loop_wall_1 = new FVM::Operator(scalarGrid);
            FVM::Operator *Op_V_loop_wall_2 = new FVM::Operator(scalarGrid);

            Op_V_loop_wall_1->AddTerm(new FVM::IdentityTerm(scalarGrid,-1.0));
            Op_V_loop_wall_2->AddTerm(new FVM::IdentityTerm(scalarGrid,R_W));
            eqsys->SetOperator(id_V_loop_wall, id_V_loop_wall, Op_V_loop_wall_1, "V_loop_wall = R_w*I_w");
            eqsys->SetOperator(id_V_loop_wall, id_I_w, Op_V_loop_wall_2);

            // Set psi_w equation
            FVM::Operator *Op_I_w_1 = new FVM::Operator(scalarGrid);
            FVM::Operator *Op_I_w_2 = new FVM::Operator(scalarGrid);
            FVM::Operator *Op_I_w_3 = new FVM::Operator(scalarGrid);
            
            Op_I_w_1->AddTerm(new FVM::TransientTerm(scalarGrid,id_psi_wall));
            Op_I_w_2->AddTerm(new FVM::TransientTerm(scalarGrid,id_I_w, L_ext));
            Op_I_w_3->AddTerm(new FVM::TransientTerm(scalarGrid,id_I_p, L_ext));

            string psiw_desc = "psi_w = ";

            // If loop voltage is applied via a transformer...
            if (type == OptionConstants::UQTY_V_LOOP_WALL_EQN_TRANSFORMER) {
                // Add unknown for psi_trans...
                eqsys->SetUnknown(OptionConstants::UQTY_PSI_TRANS, OptionConstants::UQTY_PSI_TRANS_DESC, scalarGrid);
                eqsys->SetUnknown(OptionConstants::UQTY_V_LOOP_TRANS, OptionConstants::UQTY_V_LOOP_TRANS_DESC, scalarGrid);
                const len_t id_psi_trans = unknowns->GetUnknownID(OptionConstants::UQTY_PSI_TRANS);
                const len_t id_V_loop_trans = unknowns->GetUnknownID(OptionConstants::UQTY_V_LOOP_TRANS);

                FVM::Operator *Op_psi_w_psi_trans = new FVM::Operator(scalarGrid);

                // Add psi_trans to psi_wall equation
                Op_psi_w_psi_trans->AddTerm(new FVM::TransientTerm(scalarGrid, id_psi_trans, -1.0));
                eqsys->SetOperator(id_I_w, id_psi_trans, Op_psi_w_psi_trans);
                psiw_desc += "psi_trans ";

                // Get prescribed loop voltage...
                FVM::Interpolator1D *V_loop_trans = LoadDataT(MODULENAME, s, "V_loop_wall");

                // Add equation for psi_trans...
                FVM::Operator *Op_psi_trans_1 = new FVM::Operator(scalarGrid);
                FVM::Operator *Op_psi_trans_2 = new FVM::Operator(scalarGrid);

                Op_psi_trans_1->AddTerm(new FVM::TransientTerm(scalarGrid, id_psi_trans, -1.0));
                Op_psi_trans_2->AddTerm(new FVM::IdentityTerm(scalarGrid, 1.0));

                eqsys->SetOperator(id_psi_trans, id_psi_trans, Op_psi_trans_1, "d psi_trans / dt = V_loop_trans");
                eqsys->SetOperator(id_psi_trans, id_V_loop_trans, Op_psi_trans_2);

                // Add equation for V_loop_trans...
                FVM::Operator *Op_V_loop_trans = new FVM::Operator(scalarGrid);
                Op_V_loop_trans->AddTerm(new FVM::PrescribedParameter(scalarGrid, V_loop_trans));
                eqsys->SetOperator(id_V_loop_trans, id_V_loop_trans, Op_V_loop_trans, "Prescribed");

                // Initialize psi_trans to zero...
                eqsys->SetInitialValue(id_psi_trans, nullptr);
                // ...and V_loop_trans to whatever its prescribed value is...
                eqsys->initializer->AddRule(
                    id_V_loop_trans,
                    EqsysInitializer::INITRULE_EVAL_EQUATION
                );
            }

            psiw_desc += "- L_ext*(I_p+I_w)";

            eqsys->SetOperator(id_I_w,id_psi_wall,Op_I_w_1, psiw_desc);
            eqsys->SetOperator(id_I_w,id_I_w,Op_I_w_2);
            eqsys->SetOperator(id_I_w,id_I_p,Op_I_w_3);

            // Initialize I_w to 0
			real_t I_wall_0 = s->GetReal(MODULENAME "/I_wall_0");
            eqsys->SetInitialValue(id_I_w, &I_wall_0);
        }
    } else
        FVM::FVMException("Unrecognized equation type for '%s': %d.",
                OptionConstants::UQTY_V_LOOP_WALL, type);
    // Set equation dpsi_w/dt = V_loop_wall
    FVM::Operator *Op_psi_wall_1 = new FVM::Operator(scalarGrid);
    FVM::Operator *Op_psi_wall_2 = new FVM::Operator(scalarGrid);

    Op_psi_wall_1->AddTerm(new FVM::TransientTerm(scalarGrid,id_psi_wall,-1.0));
    Op_psi_wall_2->AddTerm(new FVM::IdentityTerm(scalarGrid));
    eqsys->SetOperator(id_psi_wall, id_V_loop_wall, Op_psi_wall_2, "dpsi_w/dt = V_loop_wall");
    eqsys->SetOperator(id_psi_wall, id_psi_wall, Op_psi_wall_1);

    // Initialize psi_w to 0
    eqsys->SetInitialValue(id_psi_wall, nullptr);

    // Regardless of setting, V_loop_wall is initialized from its equation.
    eqsys->initializer->AddRule(
        id_V_loop_wall,
        EqsysInitializer::INITRULE_EVAL_EQUATION
    );
}


/**
 * Sets the equation for the poloidal flux the wall to psi_wall = 0
 */
void SimulationGenerator::ConstructEquation_psi_wall_zero(
    EquationSystem *eqsys, Settings* /*s*/
) {
    FVM::Grid *scalarGrid   = eqsys->GetScalarGrid();
    const len_t id_psi_wall = eqsys->GetUnknownID(OptionConstants::UQTY_PSI_WALL);

    FVM::Operator *Op = new FVM::Operator(scalarGrid);
    Op->AddTerm(new FVM::ConstantParameter(scalarGrid,0.0));
    eqsys->SetOperator(id_psi_wall, id_psi_wall, Op, "zero");

    eqsys->initializer->AddRule(
        id_psi_wall,
        EqsysInitializer::INITRULE_EVAL_EQUATION
    );
}

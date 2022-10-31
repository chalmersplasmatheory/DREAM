/**
 * Definition of equations relating to the runaway electron
 * distribution function.
 */

#include <string>
#include "DREAM/DREAMException.hpp"
#include "DREAM/EquationSystem.hpp"
#include "DREAM/Equations/Kinetic/BCIsotropicSourcePXi.hpp"
#include "DREAM/Equations/Kinetic/ElectricFieldTerm.hpp"
#include "DREAM/Equations/Kinetic/ElectricFieldDiffusionTerm.hpp"
#include "DREAM/Equations/Kinetic/EnergyDiffusionTerm.hpp"
#include "DREAM/Equations/Kinetic/PitchScatterTerm.hpp"
#include "DREAM/Equations/Kinetic/SlowingDownTerm.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Equation/BoundaryConditions/PXiExternalKineticKinetic.hpp"
#include "FVM/Equation/BoundaryConditions/PXiExternalLoss.hpp"
#include "FVM/Equation/BoundaryConditions/PInternalBoundaryCondition.hpp"
#include "FVM/Equation/BoundaryConditions/XiInternalBoundaryCondition.hpp"
#include "FVM/Equation/Operator.hpp"
#include "FVM/Equation/TransientTerm.hpp"
#include "FVM/Interpolator3D.hpp"


using namespace DREAM;
using namespace std;


#define MODULENAME "eqsys/f_re"


/**
 * Define settings for the runaway electron distribution function.
 *
 * s: Settings object to define options in.
 */
void SimulationGenerator::DefineOptions_f_re(Settings *s) {
    DefineOptions_f_general(s, MODULENAME);

    s->DefineSetting(MODULENAME "/inittype", "Specifies how to initialize f_re from n_re.", (int_t)OptionConstants::UQTY_F_RE_INIT_FORWARD);
	//s->DefineSetting(MODULENAME "/initavag0", "Gamma0 parameter in analytical avalanche distribution when initializing using this distribution.", (real_t)20.0);
}

/**
 * Construct the equation for the runaway electron distribution function.
 * This method is only called if the runaway grid is enabled.
 *
 * eqsys: Equation system to put the equation in.
 * s:     Settings object describing how to construct the equations.
 */
void SimulationGenerator::ConstructEquation_f_re(
    EquationSystem *eqsys, Settings *s,
    struct OtherQuantityHandler::eqn_terms *oqty_terms,
    FVM::Operator **transport
) {
    len_t id_f_re = eqsys->GetUnknownID(OptionConstants::UQTY_F_RE);
    FVM::Grid *runawayGrid = eqsys->GetRunawayGrid();

    // EXTERNAL BOUNDARY CONDITIONS
    // Lose particles to runaway region
    bool addExternalBC = true;
    bool addInternalBC = false;

    FVM::Operator *eqn = ConstructEquation_f_general(
        s, MODULENAME, eqsys, id_f_re, runawayGrid, eqsys->GetRunawayGridType(),
        eqsys->GetRunawayCollisionHandler(), addExternalBC, addInternalBC,
        transport, &oqty_terms->f_re_advective_bc, &oqty_terms->f_re_diffusive_bc,
        &oqty_terms->f_re_ripple_Dxx, &oqty_terms->f_re_timevaryingb
    );

    // Add fluid source terms (and kinetic avalanche, if enabled)
    RunawaySourceTermHandler *rsth = ConstructRunawaySourceTermHandler(
        runawayGrid, eqsys->GetHotTailGrid(), runawayGrid, eqsys->GetFluidGrid(),
        eqsys->GetUnknownHandler(), eqsys->GetREFluid(),
        eqsys->GetIonHandler(), eqsys->GetAnalyticHottailDistribution(),oqty_terms, s
    );

    len_t id_n_re  = eqsys->GetUnknownHandler()->GetUnknownID(OptionConstants::UQTY_N_RE);
    len_t id_n_tot = eqsys->GetUnknownHandler()->GetUnknownID(OptionConstants::UQTY_N_TOT);
    len_t id_n_i   = eqsys->GetUnknownHandler()->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    len_t id_n_re_neg = 0;

    FVM::Operator *Op_nRE  = new FVM::Operator(runawayGrid);
    FVM::Operator *Op_nTot = new FVM::Operator(runawayGrid);
    FVM::Operator *Op_ni   = new FVM::Operator(runawayGrid);
    FVM::Operator *Op_nREn = nullptr;

    if (s->GetBool("eqsys/n_re/negative_re")) {
        id_n_re_neg = eqsys->GetUnknownHandler()->GetUnknownID(OptionConstants::UQTY_N_RE_NEG);
        Op_nREn = new FVM::Operator(runawayGrid);
    }

    rsth->AddToOperators(Op_nRE, Op_nTot, Op_ni, Op_nREn);

    if (!Op_nRE->IsEmpty())
        eqsys->SetOperator(id_f_re, id_n_re, Op_nRE);
    if (!Op_nTot->IsEmpty())
        eqsys->SetOperator(id_f_re, id_n_tot, Op_nTot);
    if (!Op_ni->IsEmpty())
        eqsys->SetOperator(id_f_re, id_n_i, Op_nTot);
    if (Op_nREn != nullptr && !Op_nREn->IsEmpty())
        eqsys->SetOperator(id_f_re, id_n_re_neg, Op_nREn);

    // Add kinetic-kinetic boundary condition if necessary...
    if (eqsys->HasHotTailGrid()) {
		len_t id_f_hot = eqsys->GetUnknownID(OptionConstants::UQTY_F_HOT);
        const FVM::Operator *eqn_f_hot = eqsys->GetOperator(id_f_hot, id_f_hot);
		eqn->AddBoundaryCondition(new FVM::BC::PXiExternalKineticKinetic(
			runawayGrid, eqsys->GetHotTailGrid(), runawayGrid, eqn_f_hot,
			id_f_hot, id_f_re, FVM::BC::PXiExternalKineticKinetic::TYPE_UPPER
		));
    } // else {
      //    The fluid runaway source terms are the "boundary condition" at p = pMin
      // }
    // INITIALIZATION
    // This is generally handled in 'ConstructEquation_f_general()'
    // but when the hot-tail grid is disabled and an initial n_re profile
    // is prescribed, we would like to make sure that f_re integrates
    // properly to n_re.
    //if (!eqsys->HasHotTailGrid()) {
        //const len_t id_n_re    = eqsys->GetUnknownID(OptionConstants::UQTY_N_RE);
        const len_t id_E_field = eqsys->GetUnknownID(OptionConstants::UQTY_E_FIELD);
		//const len_t id_n_i     = eqsys->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
		RunawayFluid *REFluid = eqsys->GetREFluid();

        enum OptionConstants::uqty_f_re_inittype inittype =
            (enum OptionConstants::uqty_f_re_inittype)s->GetInteger(MODULENAME "/inittype");

        eqsys->initializer->AddRule(
            id_f_re, EqsysInitializer::INITRULE_EVAL_FUNCTION,
            [id_f_re,inittype,REFluid](FVM::UnknownQuantityHandler *unknowns, real_t *finit) {
                const real_t *n_re = unknowns->GetUnknownData(OptionConstants::UQTY_N_RE);
                const real_t *E    = unknowns->GetUnknownData(OptionConstants::UQTY_E_FIELD);

                FVM::Grid *runawayGrid = unknowns->GetUnknown(id_f_re)->GetGrid();
                const len_t nr = runawayGrid->GetNr();
                for (len_t ir = 0, offset = 0; ir < nr; ir++) {
                    FVM::MomentumGrid *mg = runawayGrid->GetMomentumGrid(ir);
                    const len_t np1 = mg->GetNp1();
                    const len_t np2 = mg->GetNp2();
                    real_t dp = mg->GetDp1(0);

                    real_t VpVol = runawayGrid->GetVpVol(ir);
                    if (inittype == OptionConstants::UQTY_F_RE_INIT_ISOTROPIC) {

                        // Add an equal number of particles in every cell
                        for (len_t j = 0; j < np2; j++) {
                            real_t Vp = runawayGrid->GetVp(ir, 0, j);
                            finit[offset + j*np1] = n_re[ir]*VpVol / (2.0*dp*Vp);   // 2 = integral over xi from -1 to 1
                        }
					} else if (inittype == OptionConstants::UQTY_F_RE_INIT_AVALANCHE) {
						const real_t *p = runawayGrid->GetMomentumGrid(ir)->GetP1();
						const real_t *xi = runawayGrid->GetMomentumGrid(ir)->GetP2();
						AnalyticDistributionRE *fRE = REFluid->GetAnalyticDistributionRE();

						// NOTE: This distribution function is normalized such
						// that
						//
						//   <n_re> = integral(f_re * V', 0, inf),
						//
						// i.e. unless pmax is set sufficiently high, the
						// density moment of the discretized f_re
						// will in general be different from the quantity n_re
						for (len_t j = 0; j < np2; j++)
							for (len_t i = 0; i < np1; i++)
								finit[offset + j*np1 + i] = fRE->evaluateFullDistribution(ir, xi[j], p[i]);
                    } else {
                        len_t xiIndex = 0;
                        // Select either xi=+1 or xi=-1, depending on the sign of E
                        if (inittype == OptionConstants::UQTY_F_RE_INIT_FORWARD)
                            xiIndex = (E[ir]>=0 ? np2-1 : 0);
                        else if (inittype == OptionConstants::UQTY_F_RE_INIT_XI_NEGATIVE)
                            xiIndex = 0;
                        else if (inittype == OptionConstants::UQTY_F_RE_INIT_XI_POSITIVE)
                            xiIndex = np2-1;

                        real_t Vp  = runawayGrid->GetVp(ir, 0, xiIndex);
                        real_t dxi = mg->GetDp2(xiIndex);

                        finit[offset + xiIndex*np1] = n_re[ir]*VpVol / (dxi*dp*Vp);
                    }

                    offset += np1*np2;
                }
            },
            // Dependencies
            id_n_re, id_E_field, id_n_i
        );
    //}
}


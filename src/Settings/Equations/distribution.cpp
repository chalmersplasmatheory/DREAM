/**
 * General initialization of distribution functions.
 */

#include <string>
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/EquationSystem.hpp"
#include "DREAM/Equations/Kinetic/BCIsotropicSourcePXi.hpp"
#include "DREAM/Equations/Kinetic/ElectricFieldTerm.hpp"
#include "DREAM/Equations/Kinetic/ElectricFieldDiffusionTerm.hpp"
#include "DREAM/Equations/Kinetic/EnergyDiffusionTerm.hpp"
#include "DREAM/Equations/Kinetic/PitchScatterTerm.hpp"
#include "DREAM/Equations/Kinetic/RipplePitchScattering.hpp"
#include "DREAM/Equations/Kinetic/SlowingDownTerm.hpp"
#include "DREAM/Equations/Kinetic/SynchrotronTerm.hpp"
#include "DREAM/IO.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Equation/BoundaryConditions/PXiExternalKineticKinetic.hpp"
#include "FVM/Equation/BoundaryConditions/PXiExternalLoss.hpp"
#include "FVM/Equation/BoundaryConditions/PXiInternalTrapping.hpp"
#include "FVM/Equation/BoundaryConditions/PInternalBoundaryCondition.hpp"
#include "FVM/Equation/BoundaryConditions/XiInternalBoundaryCondition.hpp"
#include "FVM/Equation/Operator.hpp"
#include "FVM/Equation/TransientTerm.hpp"
#include "FVM/Interpolator3D.hpp"


using namespace DREAM;
using namespace std;


/**
 * Define options which are present for any distribution function.
 *
 * s:   Settings object to define settings in.
 * mod: Name of module to define settings for.
 */
void SimulationGenerator::DefineOptions_f_general(Settings *s, const string& mod) {
    // External boundary condition
    s->DefineSetting(mod + "/boundarycondition", "Type of boundary condition to use at pMax.", (int_t)FVM::BC::PXiExternalLoss::BC_PHI_CONST);

    // Flux limiter settings
    s->DefineSetting(mod + "/adv_interp/r", "Type of interpolation method to use in r-component of advection term of kinetic equation.", (int_t)FVM::AdvectionInterpolationCoefficient::AD_INTERP_CENTRED);
    s->DefineSetting(mod + "/adv_jac_mode/r", "Type of interpolation method to use in the jacobian of the r-component of advection term of kinetic equation.", (int_t)OptionConstants::AD_INTERP_JACOBIAN_LINEAR);
    s->DefineSetting(mod + "/adv_interp/p1", "Type of interpolation method to use in p1-component of advection term of kinetic equation.", (int_t)FVM::AdvectionInterpolationCoefficient::AD_INTERP_CENTRED);
    s->DefineSetting(mod + "/adv_jac_mode/p1", "Type of interpolation method to use in the jacobian of the p1-component of advection term of kinetic equation.", (int_t)OptionConstants::AD_INTERP_JACOBIAN_LINEAR);
    s->DefineSetting(mod + "/adv_interp/p2", "Type of interpolation method to use in p2-component of advection term of kinetic equation.", (int_t)FVM::AdvectionInterpolationCoefficient::AD_INTERP_CENTRED);
    s->DefineSetting(mod + "/adv_jac_mode/p2", "Type of interpolation method to use in the jacobian of the p2-component of advection term of kinetic equation.", (int_t)OptionConstants::AD_INTERP_JACOBIAN_LINEAR);
    s->DefineSetting(mod + "/adv_interp/fluxlimiterdamping", "Underrelaxation parameter that may be needed to achieve convergence with flux limiter methods", (real_t) 1.0);

    s->DefineSetting(mod + "/ripplemode", "Enables/disables pitch scattering due to the magnetic ripple", (int_t)OptionConstants::EQTERM_RIPPLE_MODE_NEGLECT);
    s->DefineSetting(mod + "/synchrotronmode", "Enables/disables synchrotron losses on the distribution function", (int_t)OptionConstants::EQTERM_SYNCHROTRON_MODE_NEGLECT);

    // Initial distribution
    DefineDataR(mod, s, "n0");
    DefineDataR(mod, s, "T0");
    DefineDataR2P(mod, s, "init");

    // Kinetic transport model
    DefineOptions_Transport(mod, s, true);

    // Approximate jacobian settings
    s->DefineSetting(mod + "/fullIonJacobian", "Enables/disables the ion jacobian.", (bool) true);
}

/**
 * Construct the equation for a general distribution function.
 *
 * s:   Object to load settings from.
 * mod: Name of module to load settings from.
 */
FVM::Operator *SimulationGenerator::ConstructEquation_f_general(
    Settings *s, const string& mod, EquationSystem *eqsys,
    len_t id_f, FVM::Grid *grid, enum OptionConstants::momentumgrid_type gridtype,
    CollisionQuantityHandler *cqty, bool addExternalBC, bool addInternalBC,
    TransportAdvectiveBC **advective_bc, TransportDiffusiveBC **diffusive_bc,
    RipplePitchScattering **ripple_Dxx, bool rescaleMaxwellian
) {
    FVM::Operator *eqn = new FVM::Operator(grid);

    // Add transient term
    eqn->AddTerm(new FVM::TransientTerm(grid, id_f));

    bool withFullIonJacobian = (bool) s->GetBool(mod + "/fullIonJacobian");

    string desc;
    bool isReducedEquation = 
        (gridtype == OptionConstants::MOMENTUMGRID_TYPE_PXI &&
        grid->GetMomentumGrid(0)->GetNp2() == 1);
    // Determine whether electric field acceleration should be
    // modelled with an advection or a diffusion term
    //
    // XXX Here we assume that all momentum grids have
    // the same grid points
    if (isReducedEquation) {        
        desc = "Reduced kinetic equation";

        eqn->AddTerm(new ElectricFieldDiffusionTerm(
            grid, cqty, eqsys->GetUnknownHandler(), withFullIonJacobian
        ));

    } else {
        desc = "3D kinetic equation";

        // Electric field term
        eqn->AddTerm(new ElectricFieldTerm(
            grid, eqsys->GetUnknownHandler(), gridtype
        ));

        // Pitch scattering term
        eqn->AddTerm(new PitchScatterTerm(
            grid, cqty, gridtype,
            eqsys->GetUnknownHandler(),
            withFullIonJacobian
        ));
        
        // Add ripple effects?
        if ((*ripple_Dxx = ConstructEquation_f_ripple(s, mod, grid, gridtype)) != nullptr)
            eqn->AddTerm(*ripple_Dxx);

        // Add trapping boundary condition which mirrors the solution in the trapping region.
        // Only affects the dynamics in inhomogeneous magnetic fields.
        if(grid->HasTrapped())
            eqn->AddBoundaryCondition(new FVM::BC::PXiInternalTrapping(grid, eqn));
    }

    // ALWAYS PRESENT
    // Slowing down term
    eqn->AddTerm(new SlowingDownTerm(
        grid, cqty, gridtype, 
        eqsys->GetUnknownHandler(),
        withFullIonJacobian
    ));

    // Energy diffusion
    eqn->AddTerm(new EnergyDiffusionTerm(
        grid, cqty, gridtype,
        eqsys->GetUnknownHandler(),
        withFullIonJacobian
    ));

    // Synchrotron losses
    enum OptionConstants::eqterm_synchrotron_mode synchmode =
        (enum OptionConstants::eqterm_synchrotron_mode)s->GetInteger(mod + "/synchrotronmode");
    if (synchmode == OptionConstants::EQTERM_SYNCHROTRON_MODE_INCLUDE)
        eqn->AddTerm(new SynchrotronTerm(
            grid, gridtype
        ));

    // Add transport term
    ConstructTransportTerm(
        eqn, mod, grid,
        gridtype, eqsys->GetUnknownHandler(),
        s, true, false, advective_bc, diffusive_bc
    );

    // EXTERNAL BOUNDARY CONDITIONS
    // Lose particles to n_re?
	if (addExternalBC) {
		enum FVM::BC::PXiExternalLoss::bc_type bc =
			(enum FVM::BC::PXiExternalLoss::bc_type)s->GetInteger(mod + "/boundarycondition");

		eqn->AddBoundaryCondition(new FVM::BC::PXiExternalLoss(
			grid, eqn, id_f, nullptr,
			FVM::BC::PXiExternalLoss::BOUNDARY_KINETIC, bc
		));
	}

    // Set interpolation scheme for advection term
    enum FVM::AdvectionInterpolationCoefficient::adv_interpolation adv_interp_r =
			(enum FVM::AdvectionInterpolationCoefficient::adv_interpolation)s->GetInteger(mod + "/adv_interp/r");
    OptionConstants::adv_jacobian_mode adv_jac_mode_r = 
            (OptionConstants::adv_jacobian_mode)s->GetInteger(mod + "/adv_jac_mode/r");
    enum FVM::AdvectionInterpolationCoefficient::adv_interpolation adv_interp_p1 =
			(enum FVM::AdvectionInterpolationCoefficient::adv_interpolation)s->GetInteger(mod + "/adv_interp/p1");
    OptionConstants::adv_jacobian_mode adv_jac_mode_p1 = 
            (OptionConstants::adv_jacobian_mode)s->GetInteger(mod + "/adv_jac_mode/p1");
    enum FVM::AdvectionInterpolationCoefficient::adv_interpolation adv_interp_p2 =
			(enum FVM::AdvectionInterpolationCoefficient::adv_interpolation)s->GetInteger(mod + "/adv_interp/p2");
    OptionConstants::adv_jacobian_mode adv_jac_mode_p2 = 
            (OptionConstants::adv_jacobian_mode)s->GetInteger(mod + "/adv_jac_mode/p2");
    real_t fluxLimiterDamping = (real_t)s->GetReal(mod + "/adv_interp/fluxlimiterdamping");
    eqn->SetAdvectionInterpolationMethod(
        adv_interp_r,  adv_jac_mode_r,  FVM::FLUXGRIDTYPE_RADIAL, 
        id_f, fluxLimiterDamping
    );
    eqn->SetAdvectionInterpolationMethod(
        adv_interp_p1, adv_jac_mode_p1, FVM::FLUXGRIDTYPE_P1, 
        id_f, fluxLimiterDamping
    );
    eqn->SetAdvectionInterpolationMethod(
        adv_interp_p2, adv_jac_mode_p2, FVM::FLUXGRIDTYPE_P2, 
        id_f, fluxLimiterDamping
    );

    // Set lower boundary condition to 'mirrored' so that interpolation coefficients can set
    // boundary condition at p=0
    if (addInternalBC)
        eqn->SetAdvectionBoundaryConditions(FVM::FLUXGRIDTYPE_P1, FVM::AdvectionInterpolationCoefficient::AD_BC_MIRRORED, FVM::AdvectionInterpolationCoefficient::AD_BC_DIRICHLET);

    eqsys->SetOperator(id_f, id_f, eqn, desc);


    // Set initial value of distribution
    //   First, we check whether the distribution has been specified numerically.
    //   If it hasn't, we prescribe a Maxwellian with the correct temperature.
    len_t nx[3];
    if (s->GetRealArray(mod + "/init/x", 3, nx, false) != nullptr) {
        FVM::Interpolator3D *interp = LoadDataR2P(mod, s, "init");
        enum FVM::Interpolator3D::momentumgrid_type momtype = GetInterp3DMomentumGridType(gridtype);
        const real_t *init = interp->Eval(grid, momtype);

        eqsys->SetInitialValue(id_f, init);

        delete [] init;
        delete interp;
    } else {
        real_t *n0 = LoadDataR(mod, grid->GetRadialGrid(), s, "n0");
        real_t *T0 = LoadDataR(mod, grid->GetRadialGrid(), s, "T0");

        ConstructEquation_f_maxwellian(id_f, eqsys, grid, n0, T0, rescaleMaxwellian);

        delete [] T0;
        delete [] n0;
    }

    return eqn;
}

/**
 * Construct and add the magnetic ripple pitch scattering
 * term (if enabled).
 */
RipplePitchScattering *SimulationGenerator::ConstructEquation_f_ripple(
    Settings *s, const std::string& mod, FVM::Grid *grid,
    enum OptionConstants::momentumgrid_type mgtype
) {
    enum OptionConstants::eqterm_ripple_mode rmode =
        (enum OptionConstants::eqterm_ripple_mode)s->GetInteger(mod + "/ripplemode");

    if (rmode == OptionConstants::EQTERM_RIPPLE_MODE_NEGLECT)
        return nullptr;

    len_t ncoils = s->GetInteger("radialgrid/ripple/ncoils");
    real_t deltaCoils = s->GetReal("radialgrid/ripple/deltacoils");

    if (ncoils == 0 && deltaCoils == 0)
        return nullptr;

    // Load in ripple
    len_t nModes_m, nModes_n;
    const int_t *m    = s->GetIntegerArray("radialgrid/ripple/m", 1, &nModes_m);
    const int_t *n    = s->GetIntegerArray("radialgrid/ripple/n", 1, &nModes_n);

    if (m == nullptr || n == nullptr)
        throw SettingsException("%s: Both 'm' and 'n' must be set.", mod.c_str());

    MultiInterpolator1D *dB_B = LoadDataIonRT("radialgrid", grid->GetRadialGrid(), s, nModes_m, "ripple");

    RipplePitchScattering *rps;
    if (ncoils > 0)
        rps = new RipplePitchScattering(
            grid, rmode, mgtype, (len_t)ncoils, nModes_m, m, n, dB_B
        );
    else
        rps = new RipplePitchScattering(
            grid, rmode, mgtype, (real_t)deltaCoils, nModes_m, m, n, dB_B
        );

    return rps;
}

/**
 * Initializes the electron distribution function as a
 * Maxwellian with the specified temperature and density
 * profiles.
 *
 * n0: Initial density profile of electrons.
 * T0: Initial temperature profile of electrons.
 */
void SimulationGenerator::ConstructEquation_f_maxwellian(
    len_t id_f, EquationSystem *eqsys, FVM::Grid *grid,
    const real_t *n0, const real_t *T0, bool rescaleMaxwellian
) {
    const len_t nr = grid->GetNr();
    real_t *init = new real_t[grid->GetNCells()];

    // Construct Maxwellian
    len_t offset = 0;
    for (len_t ir = 0; ir < nr; ir++) {
        const len_t np1 = grid->GetMomentumGrid(ir)->GetNp1();
        const len_t np2 = grid->GetMomentumGrid(ir)->GetNp2();
        const real_t *pvec = grid->GetMomentumGrid(ir)->GetP();

        // Define distribution offset vector
        real_t *f = init + offset;
        for (len_t j = 0; j < np2; j++) 
            for (len_t i = 0; i < np1; i++) {
                const real_t p = pvec[j*np1+i];
                f[j*np1 + i] = Constants::RelativisticMaxwellian(p, n0[ir], T0[ir]);
            }
    
        // Rescale initial distribution so that it integrates numerically exactly to n0
        if(rescaleMaxwellian){
            real_t normalizationFactor = n0[ir]/grid->IntegralMomentumAtRadius(ir,f);
            for (len_t i = 0; i < np1*np2; i++)
                f[i] *=  normalizationFactor;
        }
        offset += np1*np2;
    }
    eqsys->SetInitialValue(id_f, init);

    delete [] init;
}


/**
 * Common routines for adding transport terms to equations.
 */

#include "DREAM/Equations/Fluid/HeatTransportDiffusion.hpp"
#include "DREAM/Equations/Fluid/HeatTransportRechesterRosenbluth.hpp"
#include "DREAM/Equations/Kinetic/RechesterRosenbluthTransport.hpp"
#include "DREAM/Equations/TransportPrescribed.hpp"
#include "DREAM/Equations/Fluid/SvenssonTransport.hpp"
#include "DREAM/Equations/TransportBC.hpp"
#include "DREAM/Equations/FrozenCurrentCoefficient.hpp"
#include "DREAM/Equations/FrozenCurrentTransport.hpp"
#include "DREAM/IO.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Equation/Operator.hpp"


using namespace DREAM;
using namespace std;


/**
 * Define options for the transport model.
 */
void SimulationGenerator::DefineOptions_Transport(
    const string& mod, Settings *s, bool kinetic, const string& subname
) {
    // Advection
    if (kinetic)
        DefineDataTR2P(mod + "/" + subname, s, "ar");
    else
        DefineDataRT(mod + "/" + subname, s, "ar");

    // Diffusion
    if (kinetic)
        DefineDataTR2P(mod + "/" + subname, s, "drr");
    else
        DefineDataRT(mod + "/" + subname, s, "drr");

    DefineDataTR2P(mod + "/" + subname, s, "s_ar");
    DefineDataTR2P(mod + "/" + subname, s, "s_drr");
    s->DefineSetting(mod + "/" + subname + "/pstar",
        "The lower momentum bound for the (source-free) runaway radial-transport region.",
        (real_t)0.0 );    
    s->DefineSetting(mod + "/" + subname + "/interp1d_param",
        "Which parameter (time or plasma current) to use in the 1D interpolation.",
                     (int_t) OptionConstants::svensson_interp1d_param::SVENSSON_INTERP1D_TIME);
    
    // Boundary condition
    s->DefineSetting(mod + "/" + subname + "/boundarycondition", "Boundary condition to use for radial transport.", (int_t)OptionConstants::EQTERM_TRANSPORT_BC_F_0);

    // Rechester-Rosenbluth diffusion
    DefineDataRT(mod + "/" + subname, s, "dBB");

	// Frozen current mode
	s->DefineSetting(
		mod + "/" + subname + "/frozen_current_mode",
		"Type of transport operator to use in the frozen current mode.",
		(int_t)OptionConstants::eqterm_frozen_current_mode::EQTERM_FROZEN_CURRENT_MODE_DISABLED
	);
	s->DefineSetting(
		mod + "/" + subname + "/D_I_min",
		"Minimum value allowed for frozen current diffusion coefficient.",
		(real_t)0
	);
	s->DefineSetting(
		mod + "/" + subname + "/D_I_max",
		"Maximum value allowed for frozen current diffusion coefficient.",
		(real_t)1000
	);
	DefineDataT(mod + "/" + subname, s, "I_p_presc");
}

/**
 * Construct an equation term of the specified type and
 * return it.
 *
 * path:    Path in settings to load data from.
 * s:       Object to load settings from.
 * kinetic: If 'true', the term is assumed to be applied to a kinetic
 *          grid and the transport coefficient is expected to be 4D
 *          (time + radius + p1 + p2).
 */
template<typename T>
T *SimulationGenerator::ConstructTransportTerm_internal(
    const std::string& mod, FVM::Grid *grid,
    enum OptionConstants::momentumgrid_type momtype,
    Settings *s, bool kinetic,
    const std::string& subname
) {
    const real_t **x;
    const real_t *t, *r, *p1, *p2;
    len_t nt, nr, np1, np2;
    enum FVM::Interpolator3D::momentumgrid_type mtype = FVM::Interpolator3D::GRID_PXI;
    enum FVM::Interpolator3D::interp_method interp3d;

    if (kinetic) {
        struct dream_4d_data *d = LoadDataTR2P(mod, s, subname);

        x   = const_cast<const real_t**>(d->x);
        t   = d->t;
        r   = d->r;
        p1  = d->p1;
        p2  = d->p2;
        nt  = d->nt;
        nr  = d->nr;
        np1 = d->np1;
        np2 = d->np2;

        mtype = d->gridtype;
        interp3d = d->ps_interp;

        delete d;
    } else {
        struct dream_2d_data *d = LoadDataRT(mod, grid->GetRadialGrid(), s, subname, true);

        nt = d->nt;
        nr = d->nr;

        // Convert to 2D array
        x = new const real_t*[nt];
        for (len_t i = 0; i < nt; i++)
            x[i] = d->x + i*nr;

        r   = d->r;
        t   = d->t;
        p1  = nullptr;
        p2  = nullptr;
        np1 = 1;
        np2 = 1;

        interp3d = FVM::Interpolator3D::INTERP_LINEAR;
        
        delete d;
    }

    enum FVM::Interpolator3D::momentumgrid_type gridtype = FVM::Interpolator3D::GRID_PXI;

    // Determine grid type
    if (kinetic) {
        switch (momtype) {
            case OptionConstants::MOMENTUMGRID_TYPE_PXI: gridtype = FVM::Interpolator3D::GRID_PXI; break;
            case OptionConstants::MOMENTUMGRID_TYPE_PPARPPERP: gridtype = FVM::Interpolator3D::GRID_PPARPPERP; break;
            default: break;
        }
    }
    return new T(
        grid, nt, nr, np1, np2, x, t, r, p1, p2,
        mtype, gridtype, interp3d
    );
}


/**
 * For SvenssonTransport ....
 */
template<typename T>
T *SimulationGenerator::ConstructSvenssonTransportTerm_internal(
    const std::string& mod, FVM::Grid *grid,
    EquationSystem *eqsys, Settings *s,
    const std::string& subname
) {
    real_t pStar=s->GetReal(mod + "/pstar");

    enum OptionConstants::svensson_interp1d_param interp1dParam_input  =
        static_cast<OptionConstants::svensson_interp1d_param>(s->GetInteger(mod+"/interp1d_param"));
    enum T::svensson_interp1d_param interp1dParam;
    switch (interp1dParam_input){
        case OptionConstants::SVENSSON_INTERP1D_TIME: interp1dParam = T::svensson_interp1d_param::TIME;
            break;
        case OptionConstants::SVENSSON_INTERP1D_IP: interp1dParam = T::svensson_interp1d_param::IP;
            break;
        default: throw SettingsException("SvenssonTransport: Unrecognized value for interp1d_param: %d", interp1dParam_input);
    }

    FVM::UnknownQuantityHandler *unknowns = eqsys->GetUnknownHandler();
    RunawayFluid *REFluid = eqsys->GetREFluid();

    struct dream_4d_data *data4D = LoadDataTR2P(mod, s, subname);
    T *t = new T(grid, pStar, interp1dParam, unknowns, REFluid, data4D);

    delete data4D;

    return t;
}

/**
 * Construct the equation for the frozen current coefficient 'D_I'.
 */
void SimulationGenerator::ConstructEquation_D_I(
	EquationSystem *eqsys, Settings *s, const string& path
) {
	if (eqsys->GetUnknownHandler()->HasUnknown(OptionConstants::UQTY_D_I))
		return;

	// Define unknown quantity
	eqsys->SetUnknown(
		OptionConstants::UQTY_D_I,
		OptionConstants::UQTY_D_I_DESC,
		eqsys->GetScalarGrid()
	);

	FVM::Grid *scalarGrid = eqsys->GetScalarGrid();
	FVM::Grid *fluidGrid = eqsys->GetFluidGrid();

	FVM::Interpolator1D *I_p_presc = LoadDataT(path, s, "I_p_presc");
	real_t D_I_min = s->GetReal(path + "/D_I_min");
	real_t D_I_max = s->GetReal(path + "/D_I_max");

	FVM::Operator *eqn = new FVM::Operator(scalarGrid);
	FrozenCurrentCoefficient *fcc =
		new FrozenCurrentCoefficient(
			scalarGrid, fluidGrid, I_p_presc, eqsys->GetUnknownHandler(),
			D_I_min, D_I_max
		);
	eqn->AddTerm(fcc);

	eqsys->SetOperator(
		OptionConstants::UQTY_D_I,
		OptionConstants::UQTY_D_I,
		eqn,
		"I_p = I_p_presc",
		true	// Solved externally
	);

	real_t v = 0;
	eqsys->SetInitialValue(eqsys->GetUnknownID(OptionConstants::UQTY_D_I), &v);
}

/**
 * Construct the transport term(s) to add to the given operator.
 *
 * oprtr:        Operator to add the transport term to.
 * mod:          Name of module to load settings from.
 * grid:         Grid on which the operator will be defined.
 * momtype:      Type of momentum grid.
 * unknowns:     Unknown quantity handler.
 * s:            Object to load settings from.
 * kinetic:      If 'true', the term is assumed to be applied to a kinetic
 *               grid and the transport coefficient is expected to be 4D
 *               (time + radius + p1 + p2).
 * heat:         Indicates that the quantity to which this operator is
 *               applied represents heat (i.e. a temperature) and that
 *               operators for heat transport should be used where available.
 * advective_bc: If not 'nullptr', sets this pointer to the newly created
 *               advective B.C. (if any, otherwise nullptr). This can be used
 *               to gain access to the B.C. in, for example, the
 *               OtherQuantityHandler.
 * diffusive_bc: If not 'nullptr', sets this pointer to the newly created
 *               diffusive B.C. (if any, otherwise nullptr). This can be used
 *               to gain access to the B.C. in, for example, the
 *               OtherQuantityHandler.
 * subname:      Name of section in the settings module which the transport
 *               settings are stored.
 * 
 * returns: true if non-zero transport, otherwise false
 */
bool SimulationGenerator::ConstructTransportTerm(
    FVM::Operator *oprtr, const string& mod, FVM::Grid *grid,
    enum OptionConstants::momentumgrid_type momtype,
    EquationSystem *eqsys,
    Settings *s, bool kinetic, bool heat,
    TransportAdvectiveBC **advective_bc,
    TransportDiffusiveBC **diffusive_bc,
    struct OtherQuantityHandler::eqn_terms *oqty_terms,       // List of other quantity terms to save (only for SvenssonTransport)
    const string& subname
) {
    string path = mod + "/" + subname;

    bool hasNonTrivialTransport = false;

    auto hasCoeff = [&s,&path](const std::string& name, len_t ndim) {
        len_t ndims[4];
        const real_t *c;

        c = s->GetRealArray(path + "/" + name + "/x", ndim, ndims, false);

        return (c!=nullptr);
    };

    enum OptionConstants::eqterm_transport_bc bc =
        (enum OptionConstants::eqterm_transport_bc)s->GetInteger(path + "/boundarycondition");

    // Has advection?
    if (hasCoeff("ar", (kinetic?4:2))) {
        hasNonTrivialTransport = true;
        auto tt = ConstructTransportTerm_internal<TransportPrescribedAdvective>(
            path, grid, momtype, s, kinetic, "ar"
        );

        oprtr->AddTerm(tt);

        // Add boundary condition...
        TransportAdvectiveBC *abc=nullptr;
            ConstructTransportBoundaryCondition<TransportAdvectiveBC>(
                bc, tt, oprtr, path, grid
            );

        // Store B.C. for OtherQuantityHandler
        if (advective_bc != nullptr)
            *advective_bc = abc;
    }
    
    // Has diffusion?
    if (hasCoeff("drr", (kinetic?4:2))){
        hasNonTrivialTransport = true;
        FVM::DiffusionTerm *dt;
        if (not heat) {
            auto tt = ConstructTransportTerm_internal<TransportPrescribedDiffusive>(
                path, grid, momtype, s, kinetic, "drr"
            );

            oprtr->AddTerm(tt);
            dt = tt;
        } else {
            FVM::Interpolator1D *intp1 = LoadDataRT_intp(
                path, grid->GetRadialGrid(), s, "drr",
                true      // true: Drr is defined on r flux grid
            );
            
            HeatTransportDiffusion *tt = new HeatTransportDiffusion(
                grid, momtype, intp1, eqsys->GetUnknownHandler()
            );

            oprtr->AddTerm(tt);
            dt = tt;
        }

        // Add boundary condition...
        TransportDiffusiveBC *dbc =
            ConstructTransportBoundaryCondition<TransportDiffusiveBC>(
                bc, dt, oprtr, path, grid
            );

        // Store B.C. for OtherQuantityHandler
        if (diffusive_bc != nullptr)
            *diffusive_bc = dbc;
    }

    // Rechester-Rosenbluth diffusion?
    if (hasCoeff("dBB", 2)) {
        if (hasNonTrivialTransport)
            DREAM::IO::PrintWarning(
                DREAM::IO::WARNING_INCOMPATIBLE_TRANSPORT,
                "Rechester-Rosenbluth transport applied alongside other transport model."
            );

        if (!kinetic && !heat)
            throw SettingsException(
                "%s: Rechester-Rosenbluth diffusion can only be applied to kinetic quantities and heat.",
                path.c_str()
            );

        hasNonTrivialTransport = true;

        FVM::Interpolator1D *dBB = LoadDataRT_intp(
            mod + "/" + subname, grid->GetRadialGrid(),
            s, "dBB", true      // true: dBB is defined on r flux grid
        );

        FVM::DiffusionTerm *dt;
        if (not heat) { // Particle transport
            RechesterRosenbluthTransport *rrt = new RechesterRosenbluthTransport(
                grid, momtype, dBB
            );
            oprtr->AddTerm(rrt);

            dt = rrt;
        } else {
            HeatTransportRechesterRosenbluth *htrr = new HeatTransportRechesterRosenbluth(
                grid, momtype, dBB, eqsys->GetUnknownHandler()
            );
            oprtr->AddTerm(htrr);

            dt = htrr;
        }

        // Add boundary condition...
        TransportDiffusiveBC *dbc =
            ConstructTransportBoundaryCondition<TransportDiffusiveBC>(
                bc, dt, oprtr, path, grid
            );

        // Store B.C. for OtherQuantityHandler
        if (diffusive_bc != nullptr)
            *diffusive_bc = dbc;
    }

    bool hasSvenssonA = false;
    if (hasCoeff("s_ar", 4)) {
        hasSvenssonA = true;
        hasNonTrivialTransport = true;
        auto tt = ConstructSvenssonTransportTerm_internal<SvenssonTransportAdvectionTermA>(
            path, grid, eqsys, s, "s_ar"
            );

        if (oqty_terms != nullptr)
            oqty_terms->svensson_A = tt;

        oprtr->AddTerm(tt);

        // Add boundary condition...
		TransportAdvectiveBC *abc =
			ConstructTransportBoundaryCondition<TransportAdvectiveBC>(
				bc, tt, oprtr, path, grid
            );

        // Store B.C. for OtherQuantityHandler
        if (advective_bc != nullptr)
            *advective_bc = abc;
    }
    
    if (hasCoeff("s_drr", 4)) {
        hasNonTrivialTransport = true;
        auto tt_drr = ConstructSvenssonTransportTerm_internal<SvenssonTransportDiffusionTerm>(
            path, grid, eqsys, s, "s_drr"
            );
        auto tt_ar = ConstructSvenssonTransportTerm_internal<SvenssonTransportAdvectionTermD>(
            path, grid, eqsys, s, "s_drr"
            );

        if (oqty_terms != nullptr) {
            oqty_terms->svensson_D = tt_drr;
            oqty_terms->svensson_advD = tt_ar;
        }

        oprtr->AddTerm(tt_drr);
        oprtr->AddTerm(tt_ar);

        // Add boundary condition...
        TransportDiffusiveBC *dbc =
			ConstructTransportBoundaryCondition<TransportDiffusiveBC>(
				bc, tt_drr, oprtr, path, grid
			);

        // Store B.C. for OtherQuantityHandler
        if (diffusive_bc != nullptr)
            *diffusive_bc = dbc;

        if ( not hasSvenssonA ){
            // Add boundary condition...
			TransportAdvectiveBC *abc =
				ConstructTransportBoundaryCondition<TransportAdvectiveBC>(
					bc, tt_ar, oprtr, path, grid
                );

			// Store B.C. for OtherQuantityHandler
			if (advective_bc != nullptr)
				*advective_bc = abc;
        }
    }

	// Frozen current mode?
	if (
		s->GetInteger(path + "/frozen_current_mode", false) !=
		OptionConstants::eqterm_frozen_current_mode::EQTERM_FROZEN_CURRENT_MODE_DISABLED
	) {
		enum OptionConstants::eqterm_frozen_current_mode mode =
			(enum OptionConstants::eqterm_frozen_current_mode)
				s->GetInteger(path + "/frozen_current_mode");

		enum FrozenCurrentTransport::TransportMode ts;
		switch (mode) {
			case OptionConstants::eqterm_frozen_current_mode::EQTERM_FROZEN_CURRENT_MODE_CONSTANT:
				ts = FrozenCurrentTransport::TRANSPORT_MODE_CONSTANT;
				break;

			case OptionConstants::eqterm_frozen_current_mode::EQTERM_FROZEN_CURRENT_MODE_BETAPAR:
				ts = FrozenCurrentTransport::TRANSPORT_MODE_BETAPAR;
				break;

			default:
				throw SettingsException(
					"Frozen current transport: Unrecognized transport mode: %d.",
					mode
				);
		}

		// Add D_I to system of equations
		ConstructEquation_D_I(eqsys, s, path);

		// Add transport operator
		if (!heat) {
			FrozenCurrentTransport *fct = new FrozenCurrentTransport(
				grid, eqsys->GetUnknownHandler(), ts
			);
			oprtr->AddTerm(fct);

			// Add boundary condition
			TransportDiffusiveBC *dbc =
				ConstructTransportBoundaryCondition<TransportDiffusiveBC>(
					bc, fct, oprtr, path, grid
				);

			// Store B.C. for OtherQuantityHandler
			if (diffusive_bc != nullptr)
				*diffusive_bc = dbc;
		} else {
			// TODO Heat transport
			throw SettingsException(
				"Heat transport in frozen current mode not yet implemented."
			);
		}
	}
            
    return hasNonTrivialTransport;
}


template<class T1, class T2>
T1 *SimulationGenerator::ConstructTransportBoundaryCondition(
    enum OptionConstants::eqterm_transport_bc bc,
    T2 *transpTerm, FVM::Operator *oprtr, const string &path,
    FVM::Grid *grid
) {
    T1 *t = nullptr;
    switch (bc) {
        case OptionConstants::EQTERM_TRANSPORT_BC_CONSERVATIVE:
            // Nothing needs to be added...
            break;

        case OptionConstants::EQTERM_TRANSPORT_BC_F_0:
            t = new T1(grid, transpTerm, T1::TRANSPORT_BC_F0);
            oprtr->AddBoundaryCondition(t);
            break;

        case OptionConstants::EQTERM_TRANSPORT_BC_DF_CONST:
            t = new T1(grid, transpTerm, T1::TRANSPORT_BC_DF_CONST);
            oprtr->AddBoundaryCondition(t);
            break;

        default:
            throw SettingsException(
                "%s: Unrecognized boundary condition specified: %d.",
                path.c_str(), bc
            );
    }

    return t;
}


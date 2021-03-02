/**
 * Common routines for adding transport terms to equations.
 */

#include "DREAM/Equations/Fluid/HeatTransportDiffusion.hpp"
#include "DREAM/Equations/Fluid/HeatTransportRechesterRosenbluth.hpp"
#include "DREAM/Equations/Kinetic/RechesterRosenbluthTransport.hpp"
#include "DREAM/Equations/TransportPrescribed.hpp"
#include "DREAM/Equations/TransportBC.hpp"
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

    // Boundary condition
    s->DefineSetting(mod + "/" + subname + "/boundarycondition", "Boundary condition to use for radial transport.", (int_t)OptionConstants::EQTERM_TRANSPORT_BC_F_0);

    // Rechester-Rosenbluth diffusion
    DefineDataRT(mod + "/" + subname, s, "dBB");
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
    FVM::UnknownQuantityHandler *unknowns,
    Settings *s, bool kinetic, bool heat,
    TransportAdvectiveBC **advective_bc,
    TransportDiffusiveBC **diffusive_bc,
    const string& subname
) {
    string path = mod + "/" + subname;

    bool hasNonTrivialTransport = false;

    auto hasCoeff = [&s,&path](const std::string& name, bool kinetic) {
        len_t ndims[4];
        const real_t *c;

        c = s->GetRealArray(path + "/" + name + "/x", (kinetic?4:2), ndims, false);

        return (c!=nullptr);
    };

    enum OptionConstants::eqterm_transport_bc bc =
        (enum OptionConstants::eqterm_transport_bc)s->GetInteger(path + "/boundarycondition");

    // Has advection?
    if (hasCoeff("ar", kinetic)) {
        hasNonTrivialTransport = true;
        auto tt = ConstructTransportTerm_internal<TransportPrescribedAdvective>(
            path, grid, momtype, s, kinetic, "ar"
        );

        oprtr->AddTerm(tt);

        // Add boundary condition...
        TransportAdvectiveBC *abc=nullptr;
        switch (bc) {
            case OptionConstants::EQTERM_TRANSPORT_BC_CONSERVATIVE:
                // Nothing needs to be added...
                break;
            case OptionConstants::EQTERM_TRANSPORT_BC_F_0: {
                abc = new TransportAdvectiveBC(grid, tt);
                oprtr->AddBoundaryCondition(abc);
                break;
            }

            default:
                throw SettingsException(
                    "%s: Unrecognized boundary condition specified: %d.",
                    path.c_str(), bc
                );
        }

        // Store B.C. for OtherQuantityHandler
        if (advective_bc != nullptr)
            *advective_bc = abc;
    }
    
    // Has diffusion?
    if (hasCoeff("drr", kinetic)){
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
                grid, momtype, intp1, unknowns
            );

            oprtr->AddTerm(tt);
            dt = tt;
        }

        // Add boundary condition...
        TransportDiffusiveBC *dbc=nullptr;
        switch (bc) {
            case OptionConstants::EQTERM_TRANSPORT_BC_CONSERVATIVE:
                // Nothing needs to be added...
                break;
            case OptionConstants::EQTERM_TRANSPORT_BC_F_0:
                dbc = new TransportDiffusiveBC(grid, dt);
                oprtr->AddBoundaryCondition(dbc);
                break;

            default:
                throw SettingsException(
                    "%s: Unrecognized boundary condition specified: %d.",
                    path.c_str(), bc
                );
        }

        // Store B.C. for OtherQuantityHandler
        if (diffusive_bc != nullptr)
            *diffusive_bc = dbc;
    }

    // Rechester-Rosenbluth diffusion?
    if (hasCoeff("dBB", false)) {
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
                grid, momtype, dBB, unknowns
            );
            oprtr->AddTerm(htrr);

            dt = htrr;
        }

        // Add boundary condition...
        TransportDiffusiveBC *dbc=nullptr;
        switch (bc) {
            case OptionConstants::EQTERM_TRANSPORT_BC_CONSERVATIVE:
                // Nothing needs to be added...
                break;
            case OptionConstants::EQTERM_TRANSPORT_BC_F_0:
                dbc = new TransportDiffusiveBC(grid, dt);
                oprtr->AddBoundaryCondition(dbc);
                break;

            default:
                throw SettingsException(
                    "%s: Unrecognized boundary condition specified: %d.",
                    path.c_str(), bc
                );
        }

        // Store B.C. for OtherQuantityHandler
        if (diffusive_bc != nullptr)
            *diffusive_bc = dbc;
    }

    return hasNonTrivialTransport;
}


#ifndef _DREAM_PROCESS_SETTINGS_HPP
#define _DREAM_PROCESS_SETTINGS_HPP

#include "DREAM/ADAS.hpp"
#include "DREAM/ConvergenceChecker.hpp"
#include "DREAM/DiagonalPreconditioner.hpp"
#include "DREAM/EquationSystem.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/Equations/RunawaySourceTermHandler.hpp"
#include "DREAM/Equations/TransportBC.hpp"
#include "DREAM/MultiInterpolator1D.hpp"
#include "DREAM/NIST.hpp"
#include "DREAM/OtherQuantityHandler.hpp"
#include "DREAM/Settings/LoadData.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Settings/Settings.hpp"
#include "DREAM/Simulation.hpp"
#include "DREAM/Solver/Solver.hpp"
#include "DREAM/Solver/SolverLinearlyImplicit.hpp"
#include "DREAM/Solver/SolverNonLinear.hpp"
#include "DREAM/TimeStepper/TimeStepper.hpp"
#include "DREAM/TimeStepper/TimeStepperConstant.hpp"
#include "DREAM/TimeStepper/TimeStepperAdaptive.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/PXiGrid/PXiMomentumGrid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Interpolator1D.hpp"
#include "FVM/Interpolator3D.hpp"

namespace DREAM {
    class SimulationGenerator {
    public:
        // PUBLIC INTERFACE
        static Settings *CreateSettings() {
            Settings *s = new Settings();
            DefineOptions(s);

            return s;
        }
        static void DefineOptions(Settings*);
        static Simulation *ProcessSettings(Settings*);

        // FOR INTERNAL USE
        static EquationSystem *ConstructEquationSystem(Settings*, FVM::Grid*, FVM::Grid*,  enum OptionConstants::momentumgrid_type, FVM::Grid*, enum OptionConstants::momentumgrid_type, FVM::Grid*, ADAS*, NIST*);
        static FVM::Grid *ConstructScalarGrid();
        static FVM::Grid *ConstructHotTailGrid(Settings*, FVM::RadialGrid*, enum OptionConstants::momentumgrid_type*);
        static FVM::Grid *ConstructRunawayGrid(Settings*, FVM::RadialGrid*, FVM::Grid*, enum OptionConstants::momentumgrid_type*);
        
        static void ConstructRunawayFluid(
            FVM::Grid *g, FVM::UnknownQuantityHandler *unknowns, IonHandler *ih,
            OptionConstants::momentumgrid_type, EquationSystem*, Settings*
        );
        static RunawaySourceTermHandler *ConstructRunawaySourceTermHandler(
            FVM::Grid*, FVM::Grid*, FVM::Grid*, FVM::Grid*, FVM::UnknownQuantityHandler*,
            RunawayFluid*, IonHandler*, AnalyticDistributionHottail*, 
            struct OtherQuantityHandler::eqn_terms*, Settings *s, bool signPositive = true
        );

        static FVM::Grid *ConstructRadialGrid(Settings*);
        static FVM::RadialGrid *ConstructRadialGrid_Cylindrical(const int_t, Settings*);
        static FVM::RadialGrid *ConstructRadialGrid_ToroidalAnalytical(const int_t, Settings*);
        static FVM::RadialGrid *ConstructRadialGrid_Numerical(const int_t, Settings*);

        static FVM::PXiGrid::PXiMomentumGrid *Construct_PXiGrid(
            Settings*, const std::string&, const real_t, FVM::RadialGrid*
        );

        static enum FVM::Interpolator3D::momentumgrid_type GetInterp3DMomentumGridType(
            enum OptionConstants::momentumgrid_type
        );

        static void DefineOptions_RadialGrid(Settings*);

        static void DefineOptions_ADAS(Settings*);
        static void DefineOptions_CollisionQuantityHandler(Settings*);
        
        static void DefineOptions_RunawayFluid(Settings*);
        static void DefineOptions_EquationSystem(Settings*);
        static void DefineOptions_ElectricField(Settings*);
        static void DefineOptions_f_hot(Settings*);
        static void DefineOptions_f_general(Settings*, const std::string&);
        static void DefineOptions_f_re(Settings*);
        static void DefineOptions_f_ripple(const std::string&, Settings*);
        static void DefineOptions_HotTailGrid(Settings*);
        static void DefineOptions_Initializer(Settings*);
        static void DefineOptions_Ions(Settings*);
        static void DefineOptions_j_ohm(Settings*);
        static void DefineOptions_j_tot(Settings*);
        static void DefineOptions_KineticGrid(const std::string&, Settings*);
        static void DefineOptions_OtherQuantities(Settings*);
        static void DefineOptions_Output(Settings*);
        static void DefineOptions_n_re(Settings*);
        static void DefineOptions_RunawayGrid(Settings*);
        static void DefineOptions_Solver(Settings*);
        static void DefineOptions_T_cold(Settings*);
        static void DefineOptions_TimeStepper(Settings*);
        static void DefineOptions_Transport(const std::string&, Settings*, bool, const std::string& subname="transport");

        static void DefineToleranceSettings(const std::string&, Settings*, const std::string& name="tolerance");
        static ConvergenceChecker *LoadToleranceSettings(const std::string&, Settings*, FVM::UnknownQuantityHandler*, const std::vector<len_t>&, const std::string& name="tolerance");
        static void DefinePreconditionerSettings(Settings*);
        static DiagonalPreconditioner *LoadPreconditionerSettings(Settings*, FVM::UnknownQuantityHandler*, const std::vector<len_t>&);

        static ADAS *LoadADAS(Settings*);
        static NIST *LoadNIST(Settings*);
        static void LoadOutput(Settings*, Simulation*);
        static CollisionQuantityHandler *ConstructCollisionQuantityHandler(enum OptionConstants::momentumgrid_type, FVM::Grid *,FVM::UnknownQuantityHandler *, IonHandler *,  Settings*);
        static void ConstructEquations(EquationSystem*, Settings*, ADAS*, NIST*, struct OtherQuantityHandler::eqn_terms*);
        static real_t ConstructInitializer(EquationSystem*, Settings*);
        static void ConstructOtherQuantityHandler(EquationSystem*, Settings*, struct OtherQuantityHandler::eqn_terms*);
        static void ConstructSolver(EquationSystem*, Settings*);
        static void ConstructTimeStepper(EquationSystem*, Settings*);
        static void ConstructUnknowns(EquationSystem*, Settings*, FVM::Grid*, FVM::Grid*, FVM::Grid*, FVM::Grid*);
        
        // Routines for constructing specific equations
        static void ConstructEquation_E_field(EquationSystem*, Settings*);
        static void ConstructEquation_E_field_prescribed(EquationSystem*, Settings*);
        static void ConstructEquation_E_field_selfconsistent(EquationSystem*, Settings*);

        static void ConstructEquation_f_hot(EquationSystem*, Settings*, struct OtherQuantityHandler::eqn_terms*);
        static void ConstructEquation_f_maxwellian(const len_t, EquationSystem*, FVM::Grid*, const real_t*, const real_t*,bool);
        static void ConstructEquation_f_re(EquationSystem*, Settings*, struct OtherQuantityHandler::eqn_terms*);
        static DREAM::FVM::Operator *ConstructEquation_f_general(
            Settings*, const std::string&, DREAM::EquationSystem*, len_t, DREAM::FVM::Grid*,
            enum OptionConstants::momentumgrid_type, DREAM::CollisionQuantityHandler*,
            bool, bool, DREAM::TransportAdvectiveBC **abc=nullptr, DREAM::TransportDiffusiveBC **dbc=nullptr, 
            DREAM::RipplePitchScattering **rps=nullptr, bool rescaleMaxwellian=false
        );
        static DREAM::RipplePitchScattering *ConstructEquation_f_ripple(Settings*, const std::string&, FVM::Grid*, enum OptionConstants::momentumgrid_type);
        static void ConstructEquation_S_particle_explicit(EquationSystem*, Settings*, struct OtherQuantityHandler::eqn_terms*);
        static void ConstructEquation_S_particle_implicit(EquationSystem*, Settings*);

        static void ConstructEquation_Ions(EquationSystem*, Settings*, ADAS*);
        static void ConstructEquation_Ion_Ni(EquationSystem*, Settings*);
        static void ConstructEquation_T_i(EquationSystem*, Settings*);
        static void ConstructEquation_T_i_trivial(EquationSystem*, Settings*);
        static void ConstructEquation_T_i_selfconsistent(EquationSystem*, Settings*);

        static void ConstructEquation_n_cold(EquationSystem*, Settings*);
        static void ConstructEquation_n_cold_prescribed(EquationSystem*, Settings*);
        static void ConstructEquation_n_cold_selfconsistent(EquationSystem*, Settings*);

        static void ConstructEquation_n_hot(EquationSystem*, Settings*);
        static void ConstructEquation_j_hot(EquationSystem*, Settings*);
        static void ConstructEquation_j_hot_hottailMode(EquationSystem*, Settings*);
        static void ConstructEquation_j_re(EquationSystem*, Settings*);
        static void ConstructEquation_j_ohm(EquationSystem*, Settings*);
        static void ConstructEquation_j_tot(EquationSystem*, Settings*);

        static void ConstructEquation_psi_p(EquationSystem*, Settings*);
        static void ConstructEquation_psi_edge(EquationSystem*, Settings*);
        static void ConstructEquation_psi_wall_zero(EquationSystem*, Settings*);
        static void ConstructEquation_psi_wall_selfconsistent(EquationSystem*, Settings*);

        static void ConstructEquation_n_re(EquationSystem*, Settings*, struct OtherQuantityHandler::eqn_terms*);

        static void ConstructEquation_n_tot(EquationSystem*, Settings*);

        static void ConstructEquation_tau_coll(EquationSystem*);

        static void ConstructEquation_T_cold(EquationSystem*, Settings*, ADAS*, NIST*, struct OtherQuantityHandler::eqn_terms*);
        static void ConstructEquation_T_cold_prescribed(EquationSystem*, Settings*);
        static void ConstructEquation_T_cold_selfconsistent(EquationSystem*, Settings*, ADAS*, NIST*, struct OtherQuantityHandler::eqn_terms*);
        static void ConstructEquation_W_cold(EquationSystem*, Settings*);


        template<typename T>
        static T *ConstructTransportTerm_internal(
            const std::string&, FVM::Grid*, enum OptionConstants::momentumgrid_type,
            Settings*, bool, const std::string& subname="transport"
        );
        template<class T1, class T2>
        static T1 *ConstructTransportBoundaryCondition(
            enum OptionConstants::eqterm_transport_bc, T2*,
            FVM::Operator*, const std::string&, FVM::Grid*
        );
        template<typename T>
        static T *ConstructSvenssonTransportTerm_internal(const std::string&, FVM::Grid*, EquationSystem*, Settings*, const std::string& subname="transport");
        
        static bool ConstructTransportTerm(
            FVM::Operator*, const std::string&, FVM::Grid*,
            enum OptionConstants::momentumgrid_type, EquationSystem*,
            Settings*, bool, bool, DREAM::TransportAdvectiveBC** abc=nullptr,
            DREAM::TransportDiffusiveBC** dbc=nullptr,
            struct OtherQuantityHandler::eqn_terms *oqty_terms=nullptr,
            const std::string& subname="transport"
        );

        // Routines for constructing time steppers
        static TimeStepperConstant *ConstructTimeStepper_constant(Settings*, FVM::UnknownQuantityHandler*);
        static TimeStepperAdaptive *ConstructTimeStepper_adaptive(Settings*, FVM::UnknownQuantityHandler*, std::vector<len_t>*);

        // Data loading routines
        static void DefineDataR(const std::string&, Settings*, const std::string& name="data");
        static real_t *LoadDataR(const std::string&, FVM::RadialGrid*, Settings*, const std::string& name="data");
        static void DefineDataT(const std::string&, Settings*, const std::string& name="data");
        static FVM::Interpolator1D *LoadDataT(const std::string&, Settings*, const std::string& name="data");
        static void DefineDataRT(const std::string&, Settings*, const std::string& name="data");
        static struct dream_2d_data *LoadDataRT(const std::string&, FVM::RadialGrid*, Settings*, const std::string& name="data", const bool rFluxGrid=false);
        static FVM::Interpolator1D *LoadDataRT_intp(const std::string&, FVM::RadialGrid*, Settings*, const std::string& name="data", const bool rFluxGrid=false);
        static void DefineDataR2P(const std::string&, Settings*, const std::string& name="data");
        static FVM::Interpolator3D *LoadDataR2P(const std::string&, Settings*, const std::string& name="data");
        static void DefineDataTR2P(const std::string&, Settings*, const std::string& name="data");
        static struct dream_4d_data *LoadDataTR2P(const std::string&, Settings*, const std::string& name="data");
        static void DefineDataIonR(const std::string&, Settings*, const std::string& name="data");
        static real_t *LoadDataIonR(const std::string&, FVM::RadialGrid*, Settings*, const len_t, const std::string& name="data");
        static void DefineDataIonRT(const std::string&, Settings*, const std::string& name="data");
        static MultiInterpolator1D *LoadDataIonRT(const std::string&, FVM::RadialGrid*, Settings*, const len_t, const std::string& name="data");

        static real_t *InterpolateR(
            const len_t, const real_t*, const real_t*,
            FVM::RadialGrid*, const gsl_interp_type*
        );
        static real_t *InterpolateIonR(
            FVM::RadialGrid*, const len_t, const len_t,
            const real_t*, const real_t*, const gsl_interp_type*
        );

        static len_t GetNumberOfIonChargeStates(Settings*);
        static len_t GetNumberOfIonSpecies(Settings*);

        // Routines for constructing solvers
        static SolverLinearlyImplicit *ConstructSolver_linearly_implicit(Settings*, FVM::UnknownQuantityHandler*, std::vector<UnknownQuantityEquation*>*, EquationSystem*);
        static SolverNonLinear *ConstructSolver_nonlinear(Settings*, FVM::UnknownQuantityHandler*, std::vector<UnknownQuantityEquation*>*, EquationSystem*);
    };
}

#endif/*_DREAM_PROCESS_SETTINGS_HPP*/

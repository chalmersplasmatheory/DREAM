#ifndef _DREAM_EQUATIONS_FLUID_NBI_ELECTRON_HEAT_TERM
#define _DREAM_EQUATIONS_FLUID_NBI_ELECTRON_HEAT_TERM

#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/ADAS.hpp"
#include "DREAM/Equations/Fluid/IonRateEquation.hpp"
#include "FVM/Matrix.hpp"
#include <unordered_map>
#include "FVM/Interpolator1D.hpp"

namespace DREAM
{

    class NBIElectronHeatTerm : public FVM::EquationTerm
    {

    private:
        // Grid and unknown quantity identifiers
        len_t id_ncold, id_Tcold, id_ion_density, id_ion_temperature;
        FVM::Grid *grid;
        FVM::RadialGrid *radialGrid;
        ADAS *adas;
        FVM::UnknownQuantityHandler *unknowns;

        // Beam geometry parameters
        real_t P0[3];  // Beam entry point
        real_t n[3];   // Beam direction
        real_t a[3];   // Additional beam geometry parameter
        real_t ds;     // Step length along beam
        real_t s_max;  // Maximum beam length
        real_t r_beam; // Beam radius
        real_t s_start;

        // Beam physics parameters
        real_t Ti_beam;   // Beam ion temperature
        real_t m_i_beam;  // Beam ion mass
        real_t beamPower; // Total beam power
        FVM::Interpolator1D *j_B_profile;
        real_t I_B;  // Total beam current
        real_t Z0;   // Initial charge state
        real_t Zion; // Ion charge state

        // Tokamak parameters
        real_t R0;           // Major radius
        real_t plasmaVolume; // Plasma volume

        // Grid resolution parameters
        len_t nr;             // Number of radial points
        len_t ntheta;         // Number of poloidal points
        real_t dtheta;        // Poloidal angle step
        len_t n_beam_radius;  // Number of radial points for beam
        real_t d_beam_radius; // Radial step size for beam
        len_t n_beam_theta;   // Number of poloidal points for beam
        real_t d_beam_theta;  // Poloidal step size for beam
        len_t n_beam_s;       // Number of s points for beam
        real_t d_beam_s;      // s step size for beam

        // Cached calculations
        real_t *NBIHeatTerm; // Stored heating term values
        real_t *Deposition_profile;
        real_t *Deposition_profile_times_Vprime;
        real_t depositedFraction = 1.0;
        std::vector<real_t> dPdne, dPdTe, dPdni, dPdTi;

        // Structure for storing R,Z coordinates of flux surfaces
        struct FluxSurfacePoint
        {
            real_t R, Z; // Major radius and vertical position
        };
        std::vector<std::vector<FluxSurfacePoint>> cachedFluxSurfaces;

        // Precomputed beam geometry
        real_t e1[3], e2[3]; // Orthogonal basis vectors
        real_t n_norm[3];    // Normalized beam direction

        // Private helper methods
        void PrecomputeFluxSurfaces();
        void PrecomputeBeamBasisVectors();

    public:
        NBIElectronHeatTerm(
            FVM::Grid *, FVM::UnknownQuantityHandler *, ADAS *, real_t s_max, real_t r_beam,
            const real_t P0[3], const real_t n[3],
            real_t Ti_beam, real_t m_i_beam,
            real_t beamPower,
            FVM::Interpolator1D *j_B_profile, real_t Z0, real_t Zion, real_t R0);
        ~NBIElectronHeatTerm();
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler *unknowns) override;
        virtual void SetMatrixElements(FVM::Matrix *mat, real_t *rhs) override;
        virtual bool SetJacobianBlock(const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t *x) override;
        virtual void SetVectorElements(real_t *rhs, const real_t *x) override;
        virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return 1; }

        // Calculation methods
        real_t ComputeMeanFreePath(len_t ir);
        void ComputeDepositionProfile(FVM::UnknownQuantityHandler *unknowns);
        real_t Compute_dP_derivative(len_t ir, len_t derivId, FVM::UnknownQuantityHandler *unknowns);
        int_t CalculatePencilBeamFindFlux(real_t s_B, real_t r_B, real_t theta_B);
        void CartesianToCylindrical(real_t x, real_t y, real_t z, const real_t P0[3], const real_t n[3], const real_t a[3], real_t &r, real_t &theta, real_t &s);
    };

}

#endif /*_DREAM_EQUATIONS_FLUID_NBI_ELECTRON_HEAT_TERM */
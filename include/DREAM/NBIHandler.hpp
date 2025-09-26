#ifndef _DREAM_EQUATIONS_NBI_HANDLER_HPP
#define _DREAM_EQUATIONS_NBI_HANDLER_HPP

#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/ADAS.hpp"
#include "DREAM/Equations/Fluid/IonRateEquation.hpp"
#include "FVM/Matrix.hpp"
#include <unordered_map>
#include "FVM/Interpolator1D.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/Settings/Settings.hpp"

namespace DREAM {

    class NBIHandler {

    protected:
        // Grid and unknown quantity identifiers
        len_t id_ncold, id_Tcold, id_ion_density, id_ion_temperature, nZ;
        FVM::Grid *grid;
        FVM::RadialGrid *radialGrid;
        ADAS *adas;
        FVM::UnknownQuantityHandler *unknowns;
        IonHandler *ions;

        // Beam geometry parameters
        real_t P0[3];  // Beam entry point
        real_t n[3];   // Beam direction
        real_t a[3];   // Additional beam geometry parameter
        real_t ds;     // Step length along beam
        real_t s_max;  // Maximum beam length
        real_t r_beam; // Beam radius
        real_t s_start;
        real_t s_stop;

        // Beam physics parameters
        real_t Ti_beam;   // Beam ion temperature
        real_t m_i_beam;  // Beam ion mass
        real_t beamPower; // Total beam power
        FVM::Interpolator1D *j_B_profile;
        FVM::Interpolator1D *Power_Profile; // Beam power profile
        real_t I_B;                         // Total beam current
        real_t Z0;                          // Initial charge state
        real_t Zion;                        // Ion charge state
        bool TCVGaussian;                   // Flag for Gaussian profile

        // Tokamak parameters
        real_t R0;           // Major radius
        real_t plasmaVolume; // Plasma volume
        real_t beamVolume;   // Beam volume

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
        real_t *NBIHeatTerm_e; // Stored heating term values
        real_t *NBIHeatTerm_i; // Stored heating term values
        real_t *Deposition_profile;
        real_t *Deposition_profile_times_Vprime;
        real_t *dV_beam_prime_tot;
        real_t *H_r_dTe;
        real_t *H_r_dni;
        real_t *H_r_dTi;
        real_t *H_r_dne;
        real_t depositedFraction = 1.0;
        std::vector<real_t> dPdne, dPdTe, dPdni, dPdTi;

        // Structure for storing R,Z coordinates of flux surfaces
        struct FluxSurfacePoint{
            real_t R, Z; // Major radius and vertical position
        };
        std::vector<std::vector<FluxSurfacePoint>> cachedFluxSurfaces;

        // Precomputed beam geometry
        real_t e1[3], e2[3]; // Orthogonal basis vectors
        real_t n_norm[3];    // Normalized beam direction

        // Private helper methods
        void PrecomputeFluxSurfaces();
        void PrecomputeBeamBasisVectors();
        std::vector<len_t> ir_lut; // size = n_beam_theta * n_beam_radius * n_beam_s
        inline size_t lut_index(len_t it, len_t irad, len_t is) const {
            return (static_cast<size_t>(it)*n_beam_radius + irad)*n_beam_s + is;
        }
        void PrecomputeBeamMapLUT();

        void IonElectronFractions(FVM::UnknownQuantityHandler *unknowns, len_t ir, real_t &f_i, real_t &f_e);

    public:
        NBIHandler(FVM::Grid *, ADAS *, IonHandler *ions);

        
        void ConfigureFromSettings(
            Settings *s, FVM::UnknownQuantityHandler *unknowns,
            real_t s_max, real_t r_beam,
            const real_t P0[3], const real_t nvec[3],
            real_t Ti_beam, real_t m_i_beam, real_t beamPower,
            FVM::Interpolator1D *j_B_profile,
            real_t Z0, real_t Zion, real_t R0,
            bool TCVGaussian,
            FVM::Interpolator1D *Power_Profile
        );


        ~NBIHandler();
        void Build(const real_t t, const real_t, FVM::UnknownQuantityHandler *unknowns);


        // Calculation methods
        void ComputeMeanFreePath(len_t ir, real_t ncold, real_t Tcold, real_t ni, real_t Ti, real_t &lambda_s, real_t &dlambda_dI, real_t &dlambda_dne, real_t &dI_dni, real_t &dI_dTi, real_t &dI_dTe);
        void ComputeDepositionProfile(FVM::UnknownQuantityHandler *unknowns);
        int_t CalculatePencilBeamFindFlux(real_t s_B, real_t r_B, real_t theta_B);
        void CartesianToCylindrical(real_t x, real_t y, real_t z, const real_t P0[3], const real_t n[3], real_t &r, real_t &theta, real_t &s);
        real_t Calculate_jB_IB(real_t r_B, real_t theta_B);

        const real_t* GetNBIHeatTerm_e() const { return NBIHeatTerm_e; }
        const real_t* GetNBIHeatTerm_i() const { return NBIHeatTerm_i; }
        const real_t* GetH_r_dTe() const { return H_r_dTe; }
        const real_t* GetH_r_dni() const { return H_r_dni; }
        const real_t* GetH_r_dTi() const { return H_r_dTi; }
        const real_t* GetH_r_dne() const { return H_r_dne; }



    };

}

#endif /*_DREAM_EQUATIONS_NBI_HANDLER_HPP */
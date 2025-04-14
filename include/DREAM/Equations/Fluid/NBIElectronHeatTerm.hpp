#ifndef _DREAM_EQUATIONS_FLUID_NBI_ELECTRON_HEAT_TERM
#define _DREAM_EQUATIONS_FLUID_NBI_ELECTRON_HEAT_TERM

#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/ADAS.hpp"
#include "DREAM/Equations/Fluid/IonRateEquation.hpp"
#include "FVM/Matrix.hpp"
#include <unordered_map>
#include "FVM/Interpolator1D.hpp"



namespace DREAM{

class NBIElectronHeatTerm:public FVM::EquationTerm{

private:
         // Grid and unknown quantity identifiers
    len_t id_ncold, id_Tcold, id_ion_density, id_ion_temperature;
    FVM::Grid *grid;
    FVM::RadialGrid *radialGrid;
    ADAS *adas;
    FVM::UnknownQuantityHandler *unknowns;

    // Beam geometry parameters
    std::array<real_t, 3> P0;   // Beam entry point
    std::array<real_t, 3> n;    // Beam direction
    std::array<real_t, 3> a;    // Additional beam geometry parameter
    real_t ds;                  // Step length along beam
    real_t s_max;               // Maximum beam length
    real_t r_beam;             // Beam radius

    // Beam physics parameters
    real_t Ti_beam;            // Beam ion temperature
    real_t m_i_beam;           // Beam ion mass
    real_t beamPower;          // Total beam power
    FVM::Interpolator1D *j_B_profile;
    real_t I_B;                // Total beam current
    real_t Z0;                 // Initial charge state
    real_t Zion;               // Ion charge state
    
    // Tokamak parameters
    real_t R0;                 // Major radius
    real_t plasmaVolume;       // Total plasma volume

    // Grid resolution parameters
    len_t nr;                  // Number of radial points
    len_t ntheta;             // Number of poloidal points
    len_t nphi;               // Number of toroidal points
    real_t dtheta;            // Poloidal angle step
    real_t dphi;              // Toroidal angle step

    // Cached calculations
    real_t *NBIHeatTerm;      // Stored heating term values
    
    // Structure for storing R,Z coordinates of flux surfaces
    struct FluxSurfacePoint {
        real_t R, Z;          // Major radius and vertical position
    };
    std::vector<std::vector<FluxSurfacePoint>> cachedFluxSurfaces;

    // Precomputed beam geometry
    std::array<real_t, 3> e1, e2;    // Orthogonal basis vectors
    std::array<real_t, 3> n_norm;    // Normalized beam direction

    // Private helper methods
    void PrecomputeFluxSurfaces();
    void PrecomputeBeamBasisVectors();
        

public:
        //Needed for equationterm
        NBIElectronHeatTerm(
            FVM::Grid*, FVM::UnknownQuantityHandler*, ADAS*,
            real_t ds, real_t s_max, real_t r_beam,
            std::array<real_t, 3> P0, std::array<real_t, 3> n,
            real_t Ti_beam, real_t m_i_beam,
            real_t beamPower, real_t plasmaVolume,
            FVM::Interpolator1D *j_B_profile, real_t Z0, real_t Zion, real_t R0
        );
        ~NBIElectronHeatTerm();
        void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler* unknowns) override;
        void SetMatrixElements(FVM::Matrix *mat, real_t *rhs) override;
        bool SetJacobianBlock(const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t *x) override;
        void SetVectorElements(real_t *rhs, const real_t *x) override;
        virtual len_t GetNumberOfNonZerosPerRow() const override { return 50; }
		virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return 50; }

        //Calculation methods
        real_t ComputeMeanFreePath(len_t ir, real_t ncold, real_t Tcold, real_t Ti, real_t ni);
        real_t ComputeSurvivalProbability(real_t s_B, real_t r_B, real_t theta_B, FVM::UnknownQuantityHandler* unknowns);
        real_t ComputeDepositionProfile(len_t ir, real_t ncold, real_t Tcold, real_t Ti, real_t ni, FVM::UnknownQuantityHandler* unknowns);
        std::tuple<real_t, real_t, real_t, real_t> Compute_dP_derivative(len_t ir, FVM::UnknownQuantityHandler *unknowns);
        len_t CalculatePencilBeamFindFlux(real_t s_B, real_t r_B, real_t theta_B);
        void CartesianToCylindrical(real_t x, real_t y, real_t z, const std::array<real_t, 3>& P0, const std::array<real_t, 3>& n, const std::array<real_t, 3>& a,real_t &r, real_t &theta, real_t &s);


};

}


#endif/*_DREAM_EQUATIONS_FLUID_NBI_ELECTRON_HEAT_TERM */
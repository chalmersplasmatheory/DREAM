#ifndef _DREAM_EQUATIONS_FLUID_NBI_ELECTRON_HEAT_TERM
#define _DREAM_EQUATIONS_FLUID_NBI_ELECTRON_HEAT_TERM

#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/ADAS.hpp"
#include "DREAM/Equations/Fluid/IonRateEquation.hpp"
#include "FVM/Matrix.hpp"
#include <unordered_map>





namespace DREAM{

class NBIElectronHeatTerm:public FVM::EquationTerm{

private:
        len_t id_ncold, id_Tcold, id_ion_density, id_ion_temperature;
        FVM::Grid *grid;
        FVM::RadialGrid *radialGrid;
        ADAS *adas;
        real_t *NBIHeatTerm;
        real_t ds, s_max;
        real_t Ti_beam, m_i_beam;
        len_t nr;
        real_t Z0, Zion;
        real_t r_beam;
        std::array<real_t, 3> P0, n, a;
        std::map<len_t, std::vector<real_t>> fluxSurfaceToS;
        std::vector<real_t> fluxSurfaceIntersectionPoints;
        real_t beamPower;
        real_t plasmaVolume;
        real_t j_B;
        real_t I_B;

        void PrecomputeBeamIntersections(const std::array<real_t, 3>& P0, std::array<real_t, 3>& n);


        

        


public:
        NBIElectronHeatTerm(
            FVM::Grid*, FVM::UnknownQuantityHandler*, ADAS*,
            real_t ds, real_t s_max, real_t r_beam,
            std::array<real_t, 3> P0, std::array<real_t, 3> n,
            real_t Ti_beam, real_t m_i_beam,
            real_t beamPower, real_t plasmaVolume,
            real_t j_B, real_t I_B
        );
        
        ~NBIElectronHeatTerm();
        void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler* unknowns) override;
        void SetMatrixElements(FVM::Matrix *mat, real_t *rhs) override;
        bool SetJacobianBlock(const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t *x, FVM::UnknownQuantityHandler* unknowns);
        void SetVectorElements(real_t *rhs, const real_t *x) override;

        real_t ComputeMeanFreePath(len_t ir, real_t ncold, real_t Tcold, real_t Ti, real_t ni);
        real_t ComputeSurvivalProbability(real_t ir, real_t s_B, FVM::UnknownQuantityHandler* unknowns);
        real_t ComputeDepositionProfile(real_t ir, real_t ncold, real_t Tcold, real_t Ti, real_t ni, FVM::UnknownQuantityHandler* unknowns);
        std::tuple<real_t, real_t, real_t, real_t> Compute_dP_derivative(len_t ir, FVM::UnknownQuantityHandler *unknowns);
        real_t ComputeConfigurationSpaceJacobian(len_t ir, real_t theta);
        len_t FindFluxSurfaceIndexForS(real_t s);
        void CartesianToCylindrical(real_t x, real_t y, real_t z, const std::array<real_t, 3>& P0, const std::array<real_t, 3>& n, const std::array<real_t, 3>& a,real_t &r, real_t &theta, real_t &s);


};

}


#endif/*_DREAM_EQUATIONS_FLUID_NBI_ELECTRON_HEAT_TERM */
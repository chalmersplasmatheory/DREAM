#ifndef _DREAM_EQUATIONS_FLUID_NBI_ELECTRON_HEAT_TERM
#define _DREAM_EQUATIONS_FLUID_NBI_ELECTRON_HEAT_TERM

#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"



namespace DREAM{

class NBIElectronHeatTerm:public FVM::EquationTerm{

private:
    len_t id_ncold;    //  cold electron density id
    len_t id_Tcold;    
    len_t id_vcold;   
    len_t nr;          // Number of radial grid points

    real_t beamIntensity = 10;  //  beam intensity, change
    real_t beamCurrentDensity = 10;  //  beam current density, change

    real_t *NBIHeatTerm = nullptr; //create an array for the NBIHeatTerm

    // helper functions
    real_t ComputeCylindricalCoordinates(len_t ir);
    real_t GetTotalCrossSection(real_t s, real_t ncold, real_t Tcold, real_t vcold);
    real_t ComputeMeanFreePath(real_t ncold, real_t sigma);
    real_t ComputeSurvivalProbability(real_t s, real_t lambda_s);
    real_t ComputeDepositionProfile(real_t lambda_s, real_t survivalProb);

public:
    NBIElectronHeatTerm(FVM::Grid*); //constructor

    virtual ~NBIElectronHeatTerm(); //destructor 

    virtual len_t GetNumberOfNonZerosPerRow() const ;

    virtual len_t GetNumberOfNonZerosPerRow_jac() const;

    virtual void Rebuild(const real_t, const real_t, UnknownQuantityHandler*); //rebuild

    virtual bool SetJacobianBlock(const len_t uqtyId, const len_t derivId, Matrix*, const real_t*);

    virtual void SetMatrixElements(Matrix */*mat*/, real_t *rhs);

    virtual void SetVectorElements(real_t*, const real_t*);
};

}


#endif/*_DREAM_EQUATIONS_FLUID_NBI_ELECTRON_HEAT_TERM */
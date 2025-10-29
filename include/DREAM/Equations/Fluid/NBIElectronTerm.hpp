#ifndef _DREAM_EQUATIONS_FLUID_NBI_ELECTRON_TERM_HPP
#define _DREAM_EQUATIONS_FLUID_NBI_ELECTRON_TERM_HPP

#include "FVM/Equation/EquationTerm.hpp"
#include "DREAM/NBIHandler.hpp"
#include "DREAM/IonHandler.hpp"

//#include "DREAM/UnknownQuantityHandler.hpp"


namespace DREAM {

class NBIElectronTerm : public FVM::EquationTerm {

protected:
    NBIHandler *handler;    
    len_t nr;
    FVM::UnknownQuantityHandler *unknowns;
    len_t id_ncold, id_Tcold, id_ion_density, id_ion_temperature ,id_ni;
    IonHandler* ions;
    FVM::RadialGrid *radialGrid;
    real_t* Qe_1;
    real_t* Qe_2;
    real_t* Qe_3;
    real_t* d_Qe1_d_Te;
    real_t* d_Qe1_d_ne;
    real_t* d_Qe1_d_n_ij;   
    real_t* d_Qe1_d_T_ij;
    real_t* d_Qe2_d_Te;
    real_t* d_Qe2_d_ne;
    real_t* d_Qe2_d_n_ij;   
    real_t* d_Qe2_d_T_ij;
    real_t* d_Qe3_d_Te;
    real_t* d_Qe3_d_ne;
    real_t* d_Qe3_d_n_ij;   
    real_t* d_Qe3_d_T_ij;

    len_t CalculateNonZeros() const {
        len_t nnz = 2;  // T_cold and n_cold
        for (len_t iz = 0; iz < ions->GetNZ(); iz++)
            nnz += ions->GetZ(iz) + 1;  // All charge states
        
        nnz += ions->GetNZ();  // Ion temperatures
        return nnz;
    }
public:
    NBIElectronTerm(NBIHandler *h, FVM::Grid *grid, FVM::UnknownQuantityHandler *unknowns, IonHandler* ionHandler);

    ~NBIElectronTerm() override {
        delete[] Qe_1;
        delete[] Qe_2;
        delete[] Qe_3;
        delete[] d_Qe1_d_Te;
        delete[] d_Qe1_d_ne;
        delete[] d_Qe1_d_n_ij;
        delete[] d_Qe1_d_T_ij;
        delete[] d_Qe2_d_Te;
        delete[] d_Qe2_d_ne;
        delete[] d_Qe2_d_n_ij;
        delete[] d_Qe2_d_T_ij;
        delete[] d_Qe3_d_Te;
        delete[] d_Qe3_d_ne;
        delete[] d_Qe3_d_n_ij;
        delete[] d_Qe3_d_T_ij;

    }

    virtual void Rebuild(const real_t t, const real_t dt, FVM::UnknownQuantityHandler *unknowns) override;
    virtual void SetVectorElements(real_t *rhs, const real_t *x) override;
    virtual void SetMatrixElements(FVM::Matrix *mat, real_t *rhs) override;
    virtual bool SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
    virtual len_t GetNumberOfNonZerosPerRow() const override {return CalculateNonZeros();}
    virtual len_t GetNumberOfNonZerosPerRow_jac() const override {return CalculateNonZeros();}
};
}
#endif /*_DREAM_EQUATIONS_FLUID_NBI_ELECTRON_TERM_HPP */




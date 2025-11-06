#ifndef _DREAM_EQUATIONS_FLUID_NBI_ION_TERM_HPP
#define _DREAM_EQUATIONS_FLUID_NBI_ION_TERM_HPP

#include "DREAM/NBIHandler.hpp"
#include "DREAM/Equations/Fluid/IonEquationTerm.hpp"  
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "DREAM/IonHandler.hpp"


namespace DREAM {

class NBIIonTerm : public IonEquationTerm<FVM::EquationTerm> {

    protected:
    len_t id_ncold, id_Tcold, id_ion_density, id_ion_temperature, id_ni;
    len_t iz;
    NBIHandler *handler;  
    IonHandler *ions;
    FVM::RadialGrid *radialGrid;
    len_t nr;
    FVM::UnknownQuantityHandler *unknowns;
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
    real_t ComputeWeightFactor(len_t ir, len_t iIon);


public:
    
    NBIIonTerm(NBIHandler* h, FVM::Grid* grid, IonHandler* ionHandler,
               FVM::UnknownQuantityHandler* unknowns, const len_t iIon);
    
    virtual ~NBIIonTerm() override {
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
    virtual bool SetCSJacobianBlock(const len_t /*uqtyId*/, const len_t derivId,
                        FVM::Matrix *jac, const real_t* /*x*/,
                        const len_t iIon, const len_t Z0, const len_t /*rOffset*/
    ) override;
    virtual void SetCSMatrixElements(
        FVM::Matrix *mat, real_t *rhs,
        const len_t iIon, const len_t Z0, const len_t /*rOffset*/
    ) override;
    virtual void SetCSVectorElements(
        real_t *vec, const real_t* /*x*/,
        const len_t iIon, const len_t Z0, const len_t /*rOffset*/
    ) override;
    virtual len_t GetNumberOfNonZerosPerRow() const override {return CalculateNonZeros();}
    virtual len_t GetNumberOfNonZerosPerRow_jac() const override {return CalculateNonZeros();}


};


} 
#endif /*_DREAM_EQUATIONS_FLUID_NBI_ION_TERM_HPP */


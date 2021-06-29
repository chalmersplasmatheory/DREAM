#include <iostream>

#include "FVM/config.h"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Equation/PredeterminedParameter.hpp"
#include "FVM/Grid/Grid.hpp"

#ifndef _DREAM_FVM_EQUATION_CONSTANT_SPI_VELOCITY_TERM_HPP
#define _DREAM_FVM_EQUATION_CONSTANT_SPI_VELOCITY_TERM_HPP

namespace DREAM {
    class ConstantSPIVelocityTerm : public FVM::PredeterminedParameter {
    private:
        const real_t *data=nullptr;
        len_t nShard;
    public:
        ConstantSPIVelocityTerm(FVM::Grid *g, FVM::UnknownQuantityHandler* u, const real_t *vp) : FVM::PredeterminedParameter(g){
            len_t id_Yp=u->GetUnknownID(OptionConstants::UQTY_Y_P);
            nShard=u->GetUnknown(id_Yp)->NumberOfMultiples();
            this->SetData(vp);
        }
        ~ConstantSPIVelocityTerm(){delete [] this->data;}


        const real_t *GetData() { return data; }
        void SetData(const real_t *vp, bool copy=true){
            if(copy){
                real_t *data_vp=new real_t[3*nShard];
                for(len_t ip=0;ip<nShard;ip++){
                    data_vp[3*ip]=vp[3*ip];
                    data_vp[3*ip+1]=vp[3*ip+1];
                    data_vp[3*ip+2]=vp[3*ip+2];
                }
                this->data=data_vp;
            }else{
                this->data=vp;
            }
        }
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override{}
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override{}
        virtual void SetVectorElements(real_t*, const real_t*) override{}
    };
}

#endif/*_DREAM_FVM_EQUATION_CONSTANT_SPI_VELOCITY_TERM_HPP*/

#include <iostream>

#include "FVM/config.h"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Equation/PredeterminedParameter.hpp"
#include "FVM/Grid/Grid.hpp"

#ifndef _DREAM_FVM_EQUATION_CONSTANT_SPI_VELOCITY_POSITION_TERM_HPP
#define _DREAM_FVM_EQUATION_CONSTANT_SPI_VELOCITY_POSITION_TERM_HPP

namespace DREAM {
    class ConstantSPIVelocityPositionTerm : public FVM::PredeterminedParameter {
    private:
        const real_t *initial_data=nullptr;
        real_t *currentData=nullptr;

        len_t id_vp;
        len_t nShard;
        real_t *vp=nullptr;
        
        const real_t *t_delay = nullptr;

    public:
        ConstantSPIVelocityPositionTerm(FVM::Grid *g, FVM::UnknownQuantityHandler* u, const real_t *xp,  const real_t *t_delay = nullptr): FVM::PredeterminedParameter(g){
            id_vp=u->GetUnknownID(OptionConstants::UQTY_V_P);
            len_t id_Yp=u->GetUnknownID(OptionConstants::UQTY_Y_P);
            nShard=u->GetUnknown(id_Yp)->NumberOfMultiples();
            this->currentData = new real_t[3*nShard];
            this->SetData(xp);
            if(t_delay != nullptr)
                this->SetTimeDelay(t_delay);
        }
        ~ConstantSPIVelocityPositionTerm(){
            delete [] currentData;
            delete [] initial_data;
        }


        const real_t *GetData() { return currentData; }
        void SetData(const real_t* xp_init, bool copy=true){
            if(copy){
                real_t *initial_data_xp=new real_t[3*nShard];
                for(len_t ip=0;ip<nShard;ip++){
                    initial_data_xp[3*ip]=xp_init[3*ip];
                    initial_data_xp[3*ip+1]=xp_init[3*ip+1];
                    initial_data_xp[3*ip+2]=xp_init[3*ip+2];
                }
                this->initial_data=initial_data_xp;
            }else{
                this->initial_data=xp_init;
            }
        }
        void SetTimeDelay(const real_t *t_delay, bool copy=true){
            if(copy){
                real_t *data_tmp=new real_t[nShard];
                for(len_t ip=0;ip<nShard;ip++){
                    data_tmp[ip]=t_delay[ip];
                }
                this->t_delay=data_tmp;
            }else{
                this->t_delay=t_delay;
            }
        }
        virtual void Rebuild(const real_t t, const real_t dt, FVM::UnknownQuantityHandler* u) override{
            vp=u->GetUnknownData(id_vp);
            for(len_t ip=0;ip<nShard;ip++){
                if(t+dt>t_delay[ip]){
                    currentData[3*ip]=initial_data[3*ip]+vp[3*ip]*(t+dt-t_delay[ip]);
                    currentData[3*ip+1]=initial_data[3*ip+1]+vp[3*ip+1]*(t+dt-t_delay[ip]);
                    currentData[3*ip+2]=initial_data[3*ip+2]+vp[3*ip+2]*(t+dt-t_delay[ip]);
                }else{
                    currentData[3*ip]=initial_data[3*ip];
                    currentData[3*ip+1]=initial_data[3*ip+1];
                    currentData[3*ip+2]=initial_data[3*ip+2];
                }
            }
        }
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override{}
        virtual void SetVectorElements(real_t*, const real_t*) override{}
    };
}

#endif/*_DREAM_FVM_EQUATION_CONSTANT_SPI_VELOCITY_TERM_HPP*/

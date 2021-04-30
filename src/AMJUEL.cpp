
#include <cmath>
#include "DREAM/AMJUEL.hpp"

using namespace DREAM;
using namespace std;


/**
 * ???
 *
 * Z0: ??? charge state.
 * n:  Background electron density.
 * T:  Background electron temperature.
 */
real_t AMJUEL::getRecLyOpaque(len_t Z0, real_t n, real_t T){
    if(Z0==1){
        real_t lnRecAMJUEL=0;
        for(len_t i=0;i<9;i++)
            for(len_t j=0;j<9;j++)
                lnRecAMJUEL+=recLyOpaque[i*9+j]*pow(log(T),i)*pow(log(n/1e14),j);
        return exp(lnRecAMJUEL)/1e6;
    }else{
        return 0;
    }        			
}

/**
 * ???
 *
 * Z0: ??? charge state.
 * n:  Background electron density.
 * T:  Background electron temperature.
 */
real_t AMJUEL::getRecLyOpaque_deriv_n(len_t Z0, real_t n, real_t T){
    if(Z0==1){
        real_t lnRecAMJUEL=0;
        real_t derivFactor=0;
        for(len_t i=0;i<9;i++)
            for(len_t j=0;j<9;j++){
                lnRecAMJUEL+=recLyOpaque[i*9+j]*pow(log(T),i)*pow(log(n/1e14),j);
                if(j>0)
                    derivFactor+=recLyOpaque[i*9+j]*pow(log(T),i)*j*pow(log(n/1e14),j-1)/n;
            }
        return derivFactor*exp(lnRecAMJUEL)/1e6;
    }else{
        return 0;
    }        			
}

/**
 * ???
 *
 * Z0: ??? charge state.
 * n:  Background electron density.
 * T:  Background electron temperature.
 */
real_t AMJUEL::getRecLyOpaque_deriv_T(len_t Z0, real_t n, real_t T){
    if(Z0==1){
        real_t lnRecAMJUEL=0;
        real_t derivFactor=0;
        for(len_t i=0;i<9;i++)
            for(len_t j=0;j<9;j++){
                lnRecAMJUEL+=recLyOpaque[i*9+j]*pow(log(T),i)*pow(log(n/1e14),j);
                if(i>0)
                    derivFactor+=recLyOpaque[i*9+j]*i*pow(log(T),i-1)/T*pow(log(n/1e14),j);
            }
        return derivFactor*exp(lnRecAMJUEL)/1e6;
    }else{
        return 0;
    }        			
}

/**
 * ???
 *
 * Z0: ??? charge state.
 * n:  Background electron density.
 * T:  Background electron temperature.
 */
real_t AMJUEL::getRecRadLyOpaque(len_t Z0, real_t n, real_t T){
    if(Z0==1){
        real_t lnRecRadAMJUEL=0;
        for(len_t i=0;i<9;i++)
            for(len_t j=0;j<9;j++)
                lnRecRadAMJUEL+=recRadLyOpaque[i*9+j]*pow(log(T),i)*pow(log(n/1e14),j);
        return Constants::ec*exp(lnRecRadAMJUEL)/1e6;
    }else{
        return 0;
    }        			
}

/**
 * ???
 *
 * Z0: ??? charge state.
 * n:  Background electron density.
 * T:  Background electron temperature.
 */
real_t AMJUEL::getRecRadLyOpaque_deriv_n(len_t Z0, real_t n, real_t T){
    if(Z0==1){
        real_t lnRecRadAMJUEL=0;
        real_t derivFactor=0;
        for(len_t i=0;i<9;i++)
            for(len_t j=0;j<9;j++){
                lnRecRadAMJUEL+=recRadLyOpaque[i*9+j]*pow(log(T),i)*pow(log(n/1e14),j);
                if(j>0)
                    derivFactor+=recRadLyOpaque[i*9+j]*pow(log(T),i)*j*pow(log(n/1e14),j-1)/n;
            }
        return Constants::ec*derivFactor*exp(lnRecRadAMJUEL)/1e6;
    }else{
        return 0;
    }        			
}

/**
 * ???
 *
 * Z0: ??? charge state.
 * n:  Background electron density.
 * T:  Background electron temperature.
 */
real_t AMJUEL::getRecRadLyOpaque_deriv_T(len_t Z0, real_t n, real_t T){
    if(Z0==1){
        real_t lnRecRadAMJUEL=0;
        real_t derivFactor=0;
        for(len_t i=0;i<9;i++)
            for(len_t j=0;j<9;j++){
                lnRecRadAMJUEL+=recRadLyOpaque[i*9+j]*pow(log(T),i)*pow(log(n/1e14),j);
                if(i>0)
                    derivFactor+=recRadLyOpaque[i*9+j]*i*pow(log(T),i-1)/T*pow(log(n/1e14),j);
            }
        return Constants::ec*derivFactor*exp(lnRecRadAMJUEL)/1e6;
    }else{
        return 0;
    }        			
}

/**
 * ???
 *
 * Z0: ??? charge state.
 * n:  Background electron density.
 * T:  Background electron temperature.
 */
real_t AMJUEL::getIonizLyOpaque(len_t Z0, real_t n, real_t T){
    if(Z0==0){
        real_t lnIonizAMJUEL=0;
        for(len_t i=0;i<9;i++)
            for(len_t j=0;j<9;j++)
                lnIonizAMJUEL+=ionizLyOpaque[i*9+j]*pow(log(T),i)*pow(log(n/1e14),j);
        return exp(lnIonizAMJUEL)/1e6;
    }else{
        return 0;
    }			
}

/**
 * ???
 *
 * Z0: ??? charge state.
 * n:  Background electron density.
 * T:  Background electron temperature.
 */
real_t AMJUEL::getIonizLyOpaque_deriv_n(len_t Z0, real_t n, real_t T){
    if(Z0==0){
        real_t lnIonizAMJUEL=0;
        real_t derivFactor=0;
        for(len_t i=0;i<9;i++)
            for(len_t j=0;j<9;j++){
                lnIonizAMJUEL+=ionizLyOpaque[i*9+j]*pow(log(T),i)*pow(log(n/1e14),j);
                if(j>0)
                    derivFactor+=ionizLyOpaque[i*9+j]*pow(log(T),i)*j*pow(log(n/1e14),j-1)/n;
            }
        return derivFactor*exp(lnIonizAMJUEL)/1e6;
    }else{
        return 0;
    }			
}

/**
 * ???
 *
 * Z0: ??? charge state.
 * n:  Background electron density.
 * T:  Background electron temperature.
 */
real_t AMJUEL::getIonizLyOpaque_deriv_T(len_t Z0, real_t n, real_t T){
    if(Z0==0){
        real_t lnIonizAMJUEL=0;
        real_t derivFactor=0;
        for(len_t i=0;i<9;i++)
            for(len_t j=0;j<9;j++){
                lnIonizAMJUEL+=ionizLyOpaque[i*9+j]*pow(log(T),i)*pow(log(n/1e14),j);
                if(i>0)
                    derivFactor+=ionizLyOpaque[i*9+j]*i*pow(log(T),i-1)/T*pow(log(n/1e14),j);
            }
        return derivFactor*exp(lnIonizAMJUEL)/1e6;
    }else{
        return 0;
    }			
}

/**
 * ???
 *
 * Z0: ??? charge state.
 * n:  Background electron density.
 * T:  Background electron temperature.
 */
real_t AMJUEL::getIonizLossLyOpaque(len_t Z0, real_t n, real_t T){
    if(Z0==0){
        real_t lnIonizLossAMJUEL=0;
        for(len_t i=0;i<9;i++)
            for(len_t j=0;j<9;j++)
                lnIonizLossAMJUEL+=ionizLossLyOpaque[i*9+j]*pow(log(T),i)*pow(log(n/1e14),j);
        return Constants::ec*exp(lnIonizLossAMJUEL)/1e6;
    }else{
        return 0;
    }        			
}

/**
 * ???
 *
 * Z0: ??? charge state.
 * n:  Background electron density.
 * T:  Background electron temperature.
 */
real_t AMJUEL::getIonizLossLyOpaque_deriv_n(len_t Z0, real_t n, real_t T){
    if(Z0==0){
        real_t lnIonizLossAMJUEL=0;
        real_t derivFactor=0;
        for(len_t i=0;i<9;i++)
            for(len_t j=0;j<9;j++){
                lnIonizLossAMJUEL+=ionizLossLyOpaque[i*9+j]*pow(log(T),i)*pow(log(n/1e14),j);
                if(j>0)
                    derivFactor+=ionizLossLyOpaque[i*9+j]*pow(log(T),i)*j*pow(log(n/1e14),j-1)/n;
            }
        return Constants::ec*derivFactor*exp(lnIonizLossAMJUEL)/1e6;
    }else{
        return 0;
    }        			
}

/**
 * ???
 *
 * Z0: ??? charge state.
 * n:  Background electron density.
 * T:  Background electron temperature.
 */
real_t AMJUEL::getIonizLossLyOpaque_deriv_T(len_t Z0, real_t n, real_t T){
    if(Z0==0){
        real_t lnIonizLossAMJUEL=0;
        real_t derivFactor=0;
        for(len_t i=0;i<9;i++)
            for(len_t j=0;j<9;j++){
                lnIonizLossAMJUEL+=ionizLossLyOpaque[i*9+j]*pow(log(T),i)*pow(log(n/1e14),j);
                if(i>0)
                    derivFactor+=ionizLossLyOpaque[i*9+j]*i*pow(log(T),i-1)/T*pow(log(n/1e14),j);
            }
        return Constants::ec*derivFactor*exp(lnIonizLossAMJUEL)/1e6;
    }else{
        return 0;
    }        			
}


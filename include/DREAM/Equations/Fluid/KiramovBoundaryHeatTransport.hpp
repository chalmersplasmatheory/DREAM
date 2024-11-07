#ifndef _DREAM_EQUATIONS_FLUID_KIRAMOV_BOUNDARY_HEAT_TRANSPORT_HPP
#define _DREAM_EQUATIONS_FLUID_KIRAMOV_BOUNDARY_HEAT_TRANSPORT_HPP

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/config.h"
#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Equations/TransportBC.hpp"



namespace DREAM {
    class KiramovBoundaryHeatTransport 
        : public FVM::AdvectionTerm {
    private:
        len_t id_T_cold, id_n_cold;
        real_t *T_cold, *n_cold;
        real_t mi; // Ion mass [kg]
        len_t Z; // Ion charge number 

        int_t minIndex = -1;
        int minZ = std::numeric_limits<int>::max(); // Initialize minZ to the maximum possible integer

        const real_t gamma = 7.0; // Heat transmission coefficient in the Kiramov model [-]
        const real_t adb_index = 1.0; // Adiabatic index [-]
        virtual void SetPartialAdvectionTerm(len_t derivId, len_t nMultiples) override;
    public:
        KiramovBoundaryHeatTransport(FVM::Grid*,FVM::UnknownQuantityHandler*, IonHandler*);
        
        
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
    }; 
    class KiramovBoundaryHeatTransportBC: public TransportAdvectiveBC{
        public:
            KiramovBoundaryHeatTransportBC(FVM::Grid *grid, FVM::UnknownQuantityHandler *unknowns, IonHandler *ions)
            : TransportAdvectiveBC(grid, new KiramovBoundaryHeatTransport(grid, unknowns, ions) , TRANSPORT_BC_F0){
            }
            virtual bool Rebuild(const real_t time, FVM::UnknownQuantityHandler *unknowns) override{
                this -> transportOperator -> Rebuild(time, 0, unknowns); 

                return true; 
            }

    };

}

#endif/*_DREAM_EQUATIONS_FLUID_KIRAMOV_BOUNDARY_HEAT_TRANSPORT_HPP*/

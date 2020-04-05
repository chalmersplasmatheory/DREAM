#ifndef _DREAM_SOLVER_SNES_HPP
#define _DREAM_SOLVER_SNES_HPP

#include <vector>
#include "DREAM/Solver/Solver.hpp"
#include "DREAM/UnknownQuantityEquation.hpp"

namespace DREAM {
    class SolverSNES : public Solver {
    private:
        // Matrix to store Jacobian in
        FVM::BlockMatrix* jacobian = nullptr;
        // Vector to store function evaluations in
        real_t* F = nullptr;

    public:
        SolverSNES(FVM::UnknownQuantityHandler*, std::vector<UnknownQuantityEquation*>*);
        virtual ~SolverSNES();

        virtual const real_t *GetSolution() const override { return F; }
        virtual void Initialize(const len_t, std::vector<len_t>&);

        virtual void Solve(const real_t, const real_t) override;
    };
}

#endif/*_DREAM_SOLVER_SNES_HPP*/

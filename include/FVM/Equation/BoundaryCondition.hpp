#ifndef _TQS_FVM_BOUNDARY_CONDITION_HPP
#define _TQS_FVM_BOUNDARY_CONDITION_HPP

namespace TQS::FVM::BC {
    class BoundaryCondition {
    protected:
        RadialGrid *grid;

    public:
        BoundaryCondition(RadialGrid *g) : grid(g);

        virtual bool GridRebuilt() {}

        virtual void Rebuild(const real_t t) = 0;
        virtual void SetMatrixElements(Matrix*) = 0;
    };
}

#endif/*_TQS_FVM_BOUNDARY_CONDITION_HPP*/

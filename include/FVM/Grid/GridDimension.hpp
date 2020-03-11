#ifndef _TQS_FVM_GRID_DIMENSION_HPP
#define _TQS_FVM_GRID_DIMENSION_HPP


namespace TQS::FVM {
    class GridDimension {
    private:
        len_t size=0;
        real_t *x=nullptr,
               *x_f=nullptr,
               *dx=nullptr,
               *dx_f=nullptr;
        real_t *weights=nullptr;

    public:
        GridDimension(/* TODO */);
        GridDimension(GridDimension*);
        ~GridDimension();

        virtual bool Rebuild(const real_t);

        len_t Size() const { return this->size; }
    };
}

#endif/*_TQS_FVM_GRID_DIMENSION_HPP*/

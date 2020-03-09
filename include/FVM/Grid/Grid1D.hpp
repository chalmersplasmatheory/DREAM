#ifndef _TQS_FVM_GRID1D_HPP
#define _TQS_FVM_GRID1D_HPP


namespace TQS::FVM {
    class Grid1D {
    private:
        len_t size=0;
        real_t *x=nullptr;
        real_t *weights=nullptr;

    public:
        Grid1D(/* TODO */);
        Grid1D(Grid1D*);
        ~Grid1D();

        virtual bool Rebuild(const real_t);

        len_t Size() const { return this->size; }
    };
}

#endif/*_TQS_FVM_GRID1D_HPP*/

#ifndef _TQS_FVM_RADIAL_GRID_GENERATOR_HPP
#define _TQS_FVM_RADIAL_GRID_GENERATOR_HPP

/***************************************************
 * Abstract base class for radial grid generators. *
 ***************************************************/

namespace TQS::FVM {
    class RadialGridGenerator {
    private:
        len_t nr=0;

    protected:
        void SetNr(const len_t n) { this->nr = n; }

    public:
        
        void GetNr() const { return this->nr; }

        virtual bool NeedsRebuild(const real_t t) = 0;
        virtual bool Rebuild(const real_t t, RadialGrid*) = 0;
    };
}

#endif/*_TQS_FVM_RADIAL_GRID_GENERATOR_HPP*/

#ifndef _DREAM_FVM_NUMERIC_B_RADIAL_GRID_GENERATOR_HPP
#define _DREAM_FVM_NUMERIC_B_RADIAL_GRID_GENERATOR_HPP

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline2d.h>
#include <string>
#include <softlib/SFile.h>
#include "FVM/Grid/RadialGridGenerator.hpp"

namespace DREAM::FVM {
    class NumericBRadialGridGenerator : public RadialGridGenerator {
    public:
        enum file_format {
            FILE_FORMAT_LUKE
        };

    private:
        bool isBuilt = false;

        real_t rMin, rMax;
        real_t *rf_provided=nullptr;
        real_t *r=nullptr, *r_f=nullptr;

        // Input data
        len_t npsi, ntheta;
        real_t Rp, Zp;      // (R,Z) coordinates of magnetic axis
        real_t *input_r=nullptr;    // Calculated minor radius in outer midplane (same size as psi)
        real_t *psi=nullptr, *theta=nullptr;        // Poloidal flux and poloidal angle grids (1D)
        real_t *R=nullptr, *Z=nullptr;          // (R,Z) coordinates of flux surfaces (2D: ntheta-by-npsi)
        real_t *dataBR=nullptr, *dataBZ=nullptr, *dataBphi=nullptr;     // Magnetic field data (2D)

        std::string name;

        // Interpolation objects for interpolating in input data
        gsl_spline2d *spline_R, *spline_Z, *spline_BR, *spline_BZ, *spline_Bphi;
        gsl_interp_accel *acc_r, *acc_theta;

		real_t *addR0DataPoint(const real_t*, const real_t*, const len_t, const len_t);
		real_t *addThetaDataPoint(const real_t*, const len_t, const len_t);

	protected:
		virtual real_t FindMagneticFieldExtremum(len_t ir, int_t sign, enum fluxGridType) override;
	
    public:
        NumericBRadialGridGenerator(
            const len_t nr, const real_t r0, const real_t ra,
            const std::string&, enum file_format frmt=FILE_FORMAT_LUKE,
			const len_t ntheta_interp=30
        );
        NumericBRadialGridGenerator(
            const real_t *r_f, const len_t nr, const std::string&,
            enum file_format frmt=FILE_FORMAT_LUKE,
			const len_t ntheta_interp=30
        );
        ~NumericBRadialGridGenerator();

        void LoadMagneticFieldData(const std::string&, enum file_format frmt=FILE_FORMAT_LUKE);
        void LoadMagneticFieldData(SFile*, enum file_format frmt=FILE_FORMAT_LUKE);

        virtual bool NeedsRebuild(const real_t) const override { return (!isBuilt); }
        virtual bool Rebuild(const real_t, RadialGrid*) override;

		real_t EvalB(const real_t, const real_t);
        
        // Override virtual methods
		virtual real_t BAtTheta(const len_t, const real_t) override;
		virtual real_t BAtTheta_f(const len_t, const real_t) override;

        virtual real_t JacobianAtTheta(const len_t ir, const real_t theta) override { return JacobianAtTheta(ir, theta, 0, 0); }
        virtual real_t JacobianAtTheta(const len_t ir, const real_t theta, const real_t, const real_t) { return JacobianAtTheta(r[ir], theta); }
        virtual real_t JacobianAtTheta_f(const len_t ir, const real_t theta) override { return JacobianAtTheta_f(ir, theta, 0, 0); }
        virtual real_t JacobianAtTheta_f(const len_t ir, const real_t theta, const real_t, const real_t) { return JacobianAtTheta(r_f[ir], theta); };
        real_t JacobianAtTheta(const real_t, const real_t, real_t *R=nullptr, real_t *dRdt=nullptr, real_t *dZdt=nullptr);

        virtual real_t ROverR0AtTheta(const len_t ir, const real_t theta) override { return ROverR0AtTheta(ir, theta, 0, 0); }
        virtual real_t ROverR0AtTheta(const len_t ir, const real_t theta, const real_t, const real_t) { return ROverR0AtTheta(r[ir], theta); }
        virtual real_t ROverR0AtTheta_f(const len_t ir, const real_t theta) override { return ROverR0AtTheta_f(ir, theta, 0, 0); }
        virtual real_t ROverR0AtTheta_f(const len_t ir, const real_t theta, const real_t, const real_t) { return ROverR0AtTheta(r_f[ir], theta); }
        real_t ROverR0AtTheta(const real_t, const real_t);

        virtual real_t NablaR2AtTheta(const len_t ir, const real_t theta) override { return NablaR2AtTheta(ir, theta, 0, 0); }
        virtual real_t NablaR2AtTheta(const len_t ir, const real_t theta, const real_t, const real_t) { return NablaR2AtTheta(r[ir], theta); }
        virtual real_t NablaR2AtTheta_f(const len_t ir, const real_t theta) override { return NablaR2AtTheta_f(ir, theta, 0, 0); }
        virtual real_t NablaR2AtTheta_f(const len_t ir, const real_t theta, const real_t, const real_t) { return NablaR2AtTheta(r_f[ir], theta); }
        real_t NablaR2AtTheta(const real_t, const real_t);

        virtual void EvaluateGeometricQuantities(
            const len_t ir, const real_t theta, real_t &B,
            real_t &Jacobian, real_t &ROverR0, real_t &NablaR2
        ) override { EvaluateGeometricQuantities(r_f[ir], theta, B, Jacobian, ROverR0, NablaR2); }
        virtual void EvaluateGeometricQuantities_fr(
            const len_t ir, const real_t theta, real_t &B,
            real_t &Jacobian, real_t &ROverR0, real_t &NablaR2
        ) override { EvaluateGeometricQuantities(r_f[ir], theta, B, Jacobian, ROverR0, NablaR2); }

        void EvaluateGeometricQuantities(
            const real_t r, const real_t theta, real_t &B,
            real_t &Jacobian, real_t &ROverR0, real_t &NablaR2
        );
    };

    struct NumericBData {
        std::string name;
        sfilesize_t npsi, ntheta;
        double Rp, Zp;
        double *psi=nullptr, *theta=nullptr;    // Coordinate grids
        double *R=nullptr, *Z=nullptr;          // Cartesian meshgrids
        double *Br=nullptr, *Bz=nullptr, *Bphi=nullptr; // Magnetic field components (meshgrids)

        ~NumericBData() {
            if (psi!=nullptr) delete [] psi;
            if (theta!=nullptr) delete [] theta;
            if (R!=nullptr) delete [] R;
            if (Z!=nullptr) delete [] Z;
            if (Br!=nullptr) delete [] Br;
            if (Bz!=nullptr) delete [] Bz;
            if (Bphi!=nullptr) delete [] Bphi;
        }
    };

    struct NumericBData *LoadNumericBFromLUKE(SFile*, const std::string& path="");
}

#endif/*_DREAM_FVM_NUMERIC_B_RADIAL_GRID_GENERATOR_HPP*/

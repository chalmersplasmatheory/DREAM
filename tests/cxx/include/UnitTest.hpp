#ifndef _UNITTEST_H
#define _UNITTEST_H

#include <string>
#include "tests/cxx/config.h"
#include "FVM/Grid/Grid.hpp"

namespace DREAMTESTS {
	class UnitTest {
		protected:
			std::string name;

		public:
            struct gridcontainer {
                std::string name;
                DREAM::FVM::Grid *grid;

                ~gridcontainer() {
                    delete this->grid;
                }
            };

			UnitTest(const std::string&);
			std::string& GetName();
			bool HasName(const std::string&);

			virtual bool Run(bool) = 0;

            struct gridcontainer *GetNextGrid(const len_t);
			//virtual DREAM::FVM::RadialGrid *InitializeGeneralGridPXi(len_t nr=10, len_t np=50, len_t nxi=30);
            virtual DREAM::FVM::Grid *InitializeGridRCylPXi(
				len_t nr=10, len_t np=50, len_t nxi=30, 
				real_t B0 = 2, real_t pMin=0, real_t pMax=10
			);
            virtual DREAM::FVM::Grid *InitializeFluidGrid(len_t nr=10, real_t B0 = 2);
			virtual DREAM::FVM::Grid *InitializeGridGeneralRPXi(
				const len_t nr, const len_t np, const len_t nxi,
                const len_t ntheta_interp = 20,
                const len_t nrProfiles=20, const real_t pMin=0,
                const real_t pMax=10
			);
            virtual DREAM::FVM::Grid *InitializeGridGeneralFluid(
                const len_t nr, const len_t ntheta_interp=20, const len_t nrProfiles=20
            );
			void PrintError(const std::string&, ...);
			void PrintOK(const std::string&, ...);
			void PrintStatus(const std::string&, ...);
			void PrintWarning(const std::string&, ...);
			real_t Rand();
			void SeedRand();

			void SaveVariable(const std::string&, real_t*, const len_t);
			void SaveVariable(const std::string&, real_t**, const len_t, const len_t);
	};
}

#endif/*_UNITTEST_H*/

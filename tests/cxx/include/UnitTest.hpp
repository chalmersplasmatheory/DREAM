#ifndef _UNITTEST_H
#define _UNITTEST_H

#include <string>
#include "tests/cxx/config.h"
#include "FVM/Grid/RadialGrid.hpp"

namespace TQSTESTS {
	class UnitTest {
		protected:
			std::string name;

		public:
            struct gridcontainer {
                std::string name;
                TQS::FVM::RadialGrid *grid;

                ~gridcontainer() {
                    delete grid;
                }
            };

			UnitTest(const std::string&);
			std::string& GetName();
			bool HasName(const std::string&);

			virtual bool Run(bool) = 0;

            struct gridcontainer *GetNextGrid(const len_t);
			virtual TQS::FVM::RadialGrid *InitializeGeneralGridPXi(len_t nr=10, len_t np=50, len_t nxi=30);

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

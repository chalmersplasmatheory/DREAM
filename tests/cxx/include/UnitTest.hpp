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
			UnitTest(const std::string&);
			std::string& GetName();
			bool HasName(const std::string&);

			virtual bool Run(bool) = 0;

			//virtual NORSE::Grid *InitializeGrid(len_t nP=50, len_t nXi=50);
			//virtual NORSE::PlasmaParameters *InitializePlasma();
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

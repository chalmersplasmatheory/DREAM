#ifndef _DREAM_STATIC_QUANTITY_HANDLER_HPP
#define _DREAM_STATIC_QUANTITY_HANDLER_HPP

#include <string>
#include <vector>
#include "FVM/Grid/Grid.hpp"
#include "DREAM/MollerKernelHandler.hpp"
#include "DREAM/StaticQuantity.hpp"


namespace DREAM {
	class StaticQuantityHandler {
	public:
		struct static_qtys {
			MollerKernelHandler *mollerKernelHandler=nullptr;
		};

	protected:
		std::vector<std::string> includeNames;
		struct static_qtys* staticData = nullptr;

		std::vector<StaticQuantity*> staticQuantities;

		FVM::Grid *fluidGrid, *hotTailGrid, *runawayGrid;

		// Helper routines
		void WriteDeltaAngleKernel(SFile*, const std::string&, const MollerDeltaAngleKernel*);

	public:
		StaticQuantityHandler(
			std::vector<std::string>&,
			struct static_qtys*,
			FVM::Grid*, FVM::Grid*, FVM::Grid*
		);
		~StaticQuantityHandler();

		void ActivateQuantities();
		StaticQuantityHandler *Clone();
		void InitializeQuantities();
		void WriteQuantities(SFile*, const std::string&);
	};
}

#endif/* _DREAM_STATIC_QUANTITY_HANDLER_HPP*/

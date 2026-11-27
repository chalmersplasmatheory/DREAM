/**
 * Implementation of a static (time-independent) quantity handler.
 */

#include <vector>
#include "DREAM/StaticQuantityHandler.hpp"
#include "DREAM/MollerKernelHandler.hpp"
#include "DREAM/Equations/Kinetic/MollerDeltaAngleKernel.hpp"
#include "DREAM/IO.hpp"


using namespace DREAM;
using namespace std;


/**
 * Constructor.
 */
StaticQuantityHandler::StaticQuantityHandler(
	vector<string> &includeNames,
	struct static_qtys *staticData,
	FVM::Grid *fluidGrid, FVM::Grid *hotTailGrid,
	FVM::Grid *runawayGrid
) : includeNames(includeNames), staticData(staticData),
	fluidGrid(fluidGrid), hotTailGrid(hotTailGrid), runawayGrid(runawayGrid) {

	this->InitializeQuantities();
}


/**
 * Destructor.
 */
StaticQuantityHandler::~StaticQuantityHandler() {
	for (auto qty : staticQuantities)
		delete qty;

	delete staticData;
}


/**
 * Activate quantities to include in the output.
 */
void StaticQuantityHandler::ActivateQuantities() {
	for (auto name : this->includeNames) {
		bool added = false;
		for (auto qty : this->staticQuantities) {
			if (qty->GetName() == name) {
				qty->Activate();
				added = true;
				break;
			}
		}

		if (!added)
			IO::PrintWarning(
				IO::WARNING_STATIC_QUANTITY_NOT_FOUND,
				"The static quantity '%s' was requested, but not defined in "
				"code. Please check the spelling of the quantity name.",
				name.c_str()
			);
	}
}


/**
 * Clone this StaticQuantityHandler.
 */
StaticQuantityHandler *StaticQuantityHandler::Clone() {
	struct static_qtys *staticData = new struct static_qtys;
	*staticData = *(this->staticData);

	return new StaticQuantityHandler(
		this->includeNames, staticData,
		this->fluidGrid, this->hotTailGrid, this->runawayGrid
	);
}


/**
 * Initialize the list of static quantities.
 */
void StaticQuantityHandler::InitializeQuantities() {
	#define DEF_QTY(NAME, DESC, FUNC) \
		this->staticQuantities.push_back(new StaticQuantity((NAME), (DESC), [this](SFile *sf, const string &path) {FUNC}))
	
	MollerKernelHandler *MKH = this->staticData->mollerKernelHandler;
	if (MKH != nullptr && MKH->HasHotHot()) {
		DEF_QTY("boltzmann/deltaAngleKernel_hot_hot", "Delta-angle kernel for hot-hot collisions",
			MollerKernelHandler *MKH = this->staticData->mollerKernelHandler;
			this->WriteDeltaAngleKernel(sf, path, MKH->AngleHotHot());
		);
	}
	
	if (MKH != nullptr && MKH->HasHotRe()) {
		DEF_QTY("boltzmann/deltaAngleKernel_hot_re", "Delta-angle kernel for hot-re collisions",
			MollerKernelHandler *MKH = this->staticData->mollerKernelHandler;
			this->WriteDeltaAngleKernel(sf, path, MKH->AngleHotRe());
		);
	}

	if (MKH != nullptr && MKH->HasReRe()) {
		DEF_QTY("boltzmann/deltaAngleKernel_re_re", "Delta-angle kernel for re-re collisions",
			MollerKernelHandler *MKH = this->staticData->mollerKernelHandler;
			this->WriteDeltaAngleKernel(sf, path, MKH->AngleReRe());
		);
	}

	this->ActivateQuantities();
}


/**
 * Write the selected static quantities to the output file.
 */
void StaticQuantityHandler::WriteQuantities(
	SFile *sf, const string &group
) {
	vector<string> existingGroups;

	if (this->includeNames.size() == 0)
		return;

	sf->CreateStruct(group);

	for (auto qty : this->staticQuantities) {
		if (!qty->IsActive())
			continue;

		// Create local path for quantity?
		auto slash = qty->GetName().find_last_of('/');
		if (slash != string::npos) {
			string groupname = qty->GetName().substr(0, slash);

			if (find(existingGroups.begin(), existingGroups.end(), groupname) == existingGroups.end()) {
				sf->CreateStruct(group+groupname);
				existingGroups.push_back(groupname);
			}
		}

		// Save quantity
		qty->WriteSFile(sf, group);
	}
}


/******************************************
 * HELPER ROUTINES                        *
 ******************************************/
/**
 * Write the delta-angle kernel to the output file.
 */
void StaticQuantityHandler::WriteDeltaAngleKernel(
	SFile *sf, const string &path, const MollerDeltaAngleKernel *kernel
) {
	len_t N = kernel->GetDeltaTableSize();
	len_t Nr = this->fluidGrid->GetNr();
	len_t NxiStarsTabulate = kernel->GetNxiStarsTabulate();
	len_t sliceSize = N/(Nr*NxiStarsTabulate);
	real_t *storage = new real_t[N];
	real_t ***deltaTable = kernel->GetDeltaTable();

	for (len_t i = 0, ir = 0; ir < Nr; ir++)
		for (len_t m = 0; m < NxiStarsTabulate; ++m)
			for (len_t j = 0; j < sliceSize; ++j, ++i)
				storage[i] = deltaTable[ir][m][j];

	sfilesize_t dims[3] = {Nr, NxiStarsTabulate, sliceSize};
	sf->WriteMultiArray(path, storage, 3, dims);

	delete [] storage;
}



/**
 * OutputGenerator implementation for saving output using an
 * SFile object.
 */

#include "DREAM/OutputGeneratorSFile.hpp"


using namespace DREAM;

/**
 * Constructor.
 */
OutputGeneratorSFile::OutputGeneratorSFile(
	FVM::Grid *grid, FVM::UnknownQuantityHandler *unknowns,
	IonHandler *ions, OtherQuantityHandler *oqty,
	EquationSystem *eqsys, const std::string& filename
) : OutputGeneratorSFile(grid, unknowns, ions, oqty, eqsys, SFile::Create(filename, SFILE_MODE_WRITE)) {}

OutputGeneratorSFile::OutputGeneratorSFile(
	FVM::Grid *grid, FVM::UnknownQuantityHandler *unknowns,
	IonHandler *ions, OtherQuantityHandler *oqty,
	EquationSystem *eqsys, SFile *sf
) : OutputGenerator(grid, unknowns, ions, oqty, eqsys), sf(sf) {}

/**
 * Destructor.
 */
OutputGeneratorSFile::~OutputGeneratorSFile() {
	delete this->sf;
}


/**
 * Save grid data.
 */
OutputGeneratorSFile::SaveGrids(const std::string& name) {
    // Save grids
    this->sf->CreateStruct(name);
    this->eqsys->SaveGrids(this->sf, name);
}
    
/**
 * Save ion meta data.
 */
OutputGeneratorSFile::SaveIonMetaData(const std::string& name) {
    // Save ion metadata
    this->sf->CreateStruct(name);
    this->eqsys->SaveIonMetaData(this->sf, name);
}

/**
 * Save other quantities.
 */
OutputGeneratorSFile::SaveOtherQuantities(const std::string& name) {
	this->oqty->SaveSFile(sf, name);
}

/**
 * Save timing information.
 */
OutputGeneratorSFile::SaveTimings(const std::string& name) {
    // Save timing information
    this->eqsys->SaveTimings(this->sf, name);
}

/**
 * Save unknown quantities.
 */
OutputGeneratorSFile::SaveUnknowns(const std::string& name) {
    // Save unknowns
    this->sf->CreateStruct(name);
    this->unknowns->SaveSFile(this->sf, name, false);
}


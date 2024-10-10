/**
 * OutputGenerator implementation for saving output using an
 * SFile object.
 */

#include <string>
#include <DREAM/IO.hpp>
#include "DREAM/OutputGeneratorSFile.hpp"
#include "DREAM/Settings/Settings.hpp"
#include "DREAM/Settings/SFile.hpp"


using namespace DREAM;
using namespace std;

/**
 * Constructor.
 */
OutputGeneratorSFile::OutputGeneratorSFile(
	EquationSystem *eqsys, const std::string& filename,
    bool savesettings
) : OutputGenerator(eqsys, savesettings), filename(filename) {}

OutputGeneratorSFile::OutputGeneratorSFile(
	EquationSystem *eqsys, SFile *sf,
    bool savesettings
) : OutputGenerator(eqsys, savesettings), sf(sf) {}

/**
 * Destructor.
 */
OutputGeneratorSFile::~OutputGeneratorSFile() {}


/**
 * Save all quantities.
 */
void OutputGeneratorSFile::Save(bool current) {
    bool close = false;
    if (this->sf == nullptr) {
        this->sf = SFile::Create(this->filename, SFILE_MODE_WRITE);
        close = true;
    }

    this->OutputGenerator::Save(current);

    if (close) {
        this->sf->Close();
        delete this->sf;

		IO::PrintInfo("Saved output file to '%s'.", this->filename.c_str());
    }
}


/**
 * Save grid data.
 *
 * name:    Name of section under which the grid data should be saved.
 * current: If true, saves only data for the current iteration/time step.
 */
void OutputGeneratorSFile::SaveGrids(const std::string& name, bool current) {
    this->sf->CreateStruct(name);

    string group;
    if (name.back() == '/')
        group = name;
    else
        group = name + "/";

    // Time grid
    const real_t *t = this->eqsys->GetTimes().data();
    if (current)
        this->sf->WriteList(group + "t", t+(this->eqsys->GetTimes().size()-1), 1);
    else
        this->sf->WriteList(group + "t", t, this->eqsys->GetTimes().size());

    FVM::RadialGrid *rgrid = this->fluidGrid->GetRadialGrid();

    // Radial grid
    const len_t nr = this->fluidGrid->GetNr();
    const real_t *r   = rgrid->GetR();
    const real_t *r_f = rgrid->GetR_f();
    const real_t *dr  = rgrid->GetDr();
    this->sf->WriteList(group + "r", r, nr);
    this->sf->WriteList(group + "r_f", r_f, nr+1);
    this->sf->WriteList(group + "dr", dr, nr);

    // Volume elements
    const real_t *VpVol = this->fluidGrid->GetVpVol();
    this->sf->WriteList(group + "VpVol", VpVol, nr);

    // Plasma size
    const real_t R0 = rgrid->GetR0();
    this->sf->WriteScalar(group + "R0", R0);
    const real_t a = rgrid->GetMinorRadius();
    this->sf->WriteScalar(group + "a", a);
    

    this->sf->CreateStruct(group + "geometry");
    string geom = group + "geometry/";
    // Geometric quantities
    const real_t *effectivePassingFraction = rgrid->GetEffPassFrac();
    this->sf->WriteList(geom + "effectivePassingFraction", effectivePassingFraction, nr);
    const real_t *xi0TrappedBoundary = rgrid->GetXi0TrappedBoundary();
    this->sf->WriteList(geom + "xi0TrappedBoundary", xi0TrappedBoundary, nr);
    const real_t *toroidalFlux = rgrid->GetToroidalFlux();
    this->sf->WriteList(geom + "toroidalFlux", toroidalFlux, nr);
    const real_t *BTorGOverR0 = rgrid->GetBTorG();
    this->sf->WriteList(geom + "GR0", BTorGOverR0, nr);
    const real_t *Bmin = rgrid->GetBmin();
    this->sf->WriteList(geom + "Bmin", Bmin, nr);
    const real_t *Bmax = rgrid->GetBmax();
    this->sf->WriteList(geom + "Bmax", Bmax, nr);
    const real_t *FSA_B2 = rgrid->GetFSA_B2();
    this->sf->WriteList(geom + "FSA_BOverBmin2", FSA_B2, nr);
    const real_t *FSA_B = rgrid->GetFSA_B();
    this->sf->WriteList(geom + "FSA_BOverBmin", FSA_B, nr);
    const real_t *FSA_1OverR2 = rgrid->GetFSA_1OverR2();
    this->sf->WriteList(geom + "FSA_R02OverR2", FSA_1OverR2, nr);
    const real_t *FSA_NablaR2OverR2 = rgrid->GetFSA_NablaR2OverR2();
    this->sf->WriteList(geom + "FSA_NablaR2_R02OverR2", FSA_NablaR2OverR2, nr);



    // Hot-tail grid
    if (this->hottailGrid != nullptr) {
        this->sf->CreateStruct(group + "hottail");
        SaveMomentumGrid(this->sf, group + "hottail/", this->hottailGrid, this->eqsys->GetHotTailGridType());
    }

    // Runaway grid
    if (this->runawayGrid != nullptr) {
        this->sf->CreateStruct(group + "runaway");
        SaveMomentumGrid(this->sf, group + "runaway/", this->runawayGrid, this->eqsys->GetRunawayGridType());
    }
}

/**
 * XXX Here we assume that all momentum grids are the same
 * Save a momentum grid to the given SFile object.
 *
 * sf:       SFile object to write grids to.
 * gridname: Full path to grid data in output file.
 * g:        Grid to save.
 */
void OutputGeneratorSFile::SaveMomentumGrid(
    SFile *sf, const string& gridname, FVM::Grid *g,
    enum OptionConstants::momentumgrid_type tp
) {
    const real_t *p1   = g->GetMomentumGrid(0)->GetP1();
    const real_t *p2   = g->GetMomentumGrid(0)->GetP2();
    const real_t *p1_f = g->GetMomentumGrid(0)->GetP1_f();
    const real_t *p2_f = g->GetMomentumGrid(0)->GetP2_f();
    const real_t *dp1  = g->GetMomentumGrid(0)->GetDp1();
    const real_t *dp2  = g->GetMomentumGrid(0)->GetDp2();
    const len_t np1    = g->GetMomentumGrid(0)->GetNp1();
    const len_t np2    = g->GetMomentumGrid(0)->GetNp2();
    const len_t nr     = g->GetRadialGrid()->GetNr();

    // Write grid type
    sf->WriteInt32List(gridname + "type", (int32_t*)&tp, 1);

    // Grid coordinates
    sf->WriteList(gridname + "p1", p1, np1);
    sf->WriteList(gridname + "p2", p2, np2);

    // Flux grid coordinates
    sf->WriteList(gridname + "p1_f", p1_f, np1+1);
    sf->WriteList(gridname + "p2_f", p2_f, np2+1);

    // Grid cell sizes
    sf->WriteList(gridname + "dp1", dp1, np1);
    sf->WriteList(gridname + "dp2", dp2, np2);

    // Grid volumes
    sfilesize_t dims[3] = {nr, np2, np1};
    sfilesize_t dims_f2[3] = {nr, np2+1, np1};
    const real_t *const* Vp = g->GetVp();
    const real_t *const* Vp_f2 = g->GetVp_f2();
    WriteCopyMultiArray(sf, gridname + "Vprime", Vp, 3, dims);
    WriteCopyMultiArray(sf, gridname + "Vprime_f2", Vp_f2, 3, dims_f2);
}

/**
 * Save ion meta data.
 */
void OutputGeneratorSFile::SaveIonMetaData(const std::string& name) {
    this->sf->CreateStruct(name);

    string group;
    if (name.back() == '/')
        group = name;
    else
        group = name + "/";

    // Get list of charges
    const len_t nZ = this->ions->GetNZ();
    const len_t *Z = this->ions->GetZs();
    sf->WriteList(group + "Z", Z, nZ);

    // Construct name vector
    const vector<string> namelist = this->ions->GetNameList();

    string names = "";
    for (auto it = namelist.begin(); it != namelist.end(); it++)
        names += *it + ";";

    sf->WriteString(group + "names", names);
}

/**
 * Save other quantities.
 */
void OutputGeneratorSFile::SaveOtherQuantities(const std::string& name) {
    this->sf->CreateStruct(name);
	this->oqty->SaveSFile(sf, name);
}

/**
 * Save settings used for this simulation
 */
void OutputGeneratorSFile::SaveSettings(const std::string& name) {
    this->sf->CreateStruct(name);
    SettingsSFile::SaveSettings(this->eqsys->GetSettings(), this->sf, name);
}

/**
 * Save data from the solver.
 */
void OutputGeneratorSFile::SaveSolverData(const std::string& name) {
    this->eqsys->SaveSolverData(this->sf, name);
}

/**
 * Save timing information.
 */
void OutputGeneratorSFile::SaveTimings(const std::string& name) {
    this->eqsys->SaveTimings(this->sf, name);
}

/**
 * Save unknown quantities.
 *
 * name:    Name of section under which the unknown data should be saved.
 * current: If true, saves only data for the current iteration/time step.
 */
void OutputGeneratorSFile::SaveUnknowns(const std::string& name, bool current) {
    this->sf->CreateStruct(name);

    if (current)
        this->unknowns->SaveSFileCurrent(this->sf, name, false);
    else
        this->unknowns->SaveSFile(this->sf, name, false);
}

/**
 * Save the given array to the specified SFile with the
 * given name. The array is copied before writing so that
 * 'v' need not be stored contiguously in memory.
 *
 * sf:   SFile object to write array to.
 * name: Name of variable in file.
 * v:    2D array data to write.
 * m:    Length of first dimension of 'v'.
 * n:    Length of second dimension of 'v'.
 */
void OutputGeneratorSFile::WriteCopyArray(
    SFile *sf, const string& name, const real_t *const* v,
    const len_t m, const len_t n
) {
    // First, copy array to make it contiguous
    real_t **t = new real_t*[m];
    t[0] = new real_t[m*n];

    for (len_t i = 1; i < m; i++)
        t[i] = t[i-1] + n;

    for (len_t i = 0; i < m; i++)
        for (len_t j = 0; j < n; j++)
            t[i][j] = v[i][j];

    sf->WriteArray(name, t, m, n);
    
    delete [] t[0];
    delete [] t;
}

/**
 * Save the given multi-dimensional array to the specified
 * SFile with the given name. The array is copied before
 * writing so that the first dimension of 'v' need not be
 * stored contiguously in memory.
 *
 * sf:   SFile object to write array to.
 * name: Name of variable in file.
 * v:    Multi-dimensional array data to write.
 * m:    Length of first dimension of 'v'.
 * n:    Length of second dimension of 'v'.
 */
void OutputGeneratorSFile::WriteCopyMultiArray(
    SFile *sf, const string& name, const real_t *const* v,
    const sfilesize_t ndim, const sfilesize_t dims[]
) {
    if (ndim < 2)
        throw OutputGeneratorException("Saving '%s': Invalid number of dimensions of multi-dimensional array. Unless this is a bug, 'SFile::WriteList()' should be used instead.", name.c_str());

    // Calculate total number of elements in array
    len_t totSize = 1;
    for (len_t i = 0; i < ndim; i++)
        totSize *= dims[i];

    // Total size of all but the first dimension
    // (e.g. of the momentum grid at each radius)
    len_t otherSize = totSize / dims[0];

    // Copy array to make it contiguous
    real_t *t = new real_t[totSize];

    for (len_t i = 0; i < dims[0]; i++)
        for (len_t j = 0; j < otherSize; j++)
            t[i*otherSize + j] = v[i][j];

    sf->WriteMultiArray(name, t, ndim, dims);
    
    delete [] t;
}


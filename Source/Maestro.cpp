
#include <Maestro.H>
#include <maestro_defaults.H>

using namespace amrex;

int Maestro::NUM_STATE = -1;
int Maestro::Rho       = -1;
int Maestro::RhoH      = -1;
int Maestro::FirstSpec = -1;
int Maestro::NumSpec   = -1;
int Maestro::Temp      = -1;
int Maestro::Pi        = -1;

// constructor - reads in parameters from inputs file
//             - sizes multilevel vectors and data structures
Maestro::Maestro ()
{

    // Geometry on all levels has been defined already.

    // read in fortran90 parameters
    ca_set_maestro_method_params();

    // FIXME
    // ca_set_method_params(QRHO,...

    // read in C++ parameters
    ReadParameters();

    // define variable mappings (Rho, RhoH, ..., NUM_STATE, etc.)
    VariableSetup();

    // set up BCRec definitions for BC types
    BCSetup();

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    int nlevs_max = maxLevel() + 1;

    istep = 0;

    t_new = 0.0;
    t_old = -1.e100;

    // set this to a large number so change_max doesn't affect the first time step
    dt = 1.e100;

    snew.resize(nlevs_max);
    sold.resize(nlevs_max);

    // stores fluxes at coarse-fine interface for synchronization
    // this will be sized "nlevs_max+1"
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
    // therefore flux_reg[0] is never actually used in the reflux operation
    flux_reg_s.resize(nlevs_max);
    flux_reg_u.resize(nlevs_max);
}

Maestro::~Maestro ()
{

}

// set covered coarse cells to be the average of overlying fine cells
void
Maestro::AverageDown (Vector<std::unique_ptr<MultiFab> >& mf)
{
    for (int lev = finest_level-1; lev >= 0; --lev)
    {
        average_down(*mf[lev+1], *mf[lev],
                     geom[lev+1], geom[lev],
                     0, mf[lev]->nComp(), refRatio(lev));
    }
}

// more flexible version of AverageDown() that lets you average down across multiple levels
void
Maestro::AverageDownTo (int crse_lev, Vector<std::unique_ptr<MultiFab> >& mf)
{
    average_down(*mf[crse_lev+1], *mf[crse_lev],
                 geom[crse_lev+1], geom[crse_lev],
                 0, mf[crse_lev]->nComp(), refRatio(crse_lev));
}

// compute the number of cells at a level
long
Maestro::CountCells (int lev)
{
    const int N = grids[lev].size();

    long cnt = 0;

#ifdef _OPENMP
#pragma omp parallel for reduction(+:cnt)
#endif
    for (int i = 0; i < N; ++i) {
        cnt += grids[lev][i].numPts();
    }

    return cnt;
}


// Delete level data
// overrides the pure virtual function in AmrCore
void
Maestro::ClearLevel (int lev)
{
    snew[lev].reset(nullptr);
    sold[lev].reset(nullptr);
    flux_reg_s[lev].reset(nullptr);
    flux_reg_u[lev].reset(nullptr);
}

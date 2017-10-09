

#include <Maestro.H>
#include <maestro_defaults.H>

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

    // ca_set_method_params(QRHO,...

    // read in C++ parameters
    ReadParameters();

    // define variable mappings (Rho, RhoH, ..., NUM_STATE, etc.)
    VariableSetup();

    bcs.resize(NUM_STATE);

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
    flux_reg.resize(nlevs_max);
}

Maestro::~Maestro ()
{

}

// read in some parameters from inputs file
void
Maestro::ReadParameters ()
{

    ParmParse pp("maestro");

#include <maestro_queries.H>

    // Get boundary conditions
    Array<int> lo_bc(AMREX_SPACEDIM), hi_bc(AMREX_SPACEDIM);
    pp.getarr("lo_bc",lo_bc,0,AMREX_SPACEDIM);
    pp.getarr("hi_bc",hi_bc,0,AMREX_SPACEDIM);
    
    // Check phys_bc against possible periodic geometry
    // if periodic, must have internal BC marked.
    //
    if (Geometry::isAnyPeriodic())
    {
        //
        // Do idiot check.  Periodic means interior in those directions.
        //
        for (int dir = 0; dir<AMREX_SPACEDIM; dir++)
        {
            if (Geometry::isPeriodic(dir))
            {
                if (lo_bc[dir] != Interior)
                {
                    std::cerr << "Maestro::ReadParameters:periodic in direction "
                              << dir << " but low BC is not Interior\n";
                    Error();
                }
                if (hi_bc[dir] != Interior)
                {
                    std::cerr << "Maestro::ReadParameters:periodic in direction "
                              << dir << " but high BC is not Interior\n";
                    Error();
                }
            }
        }
    }
    else
    {
        //
        // Do idiot check.  If not periodic, should be no interior.
        //
        for (int dir=0; dir<AMREX_SPACEDIM; dir++)
        {
            if (lo_bc[dir] == Interior)
            {
                std::cerr << "Maestro::ReadParameters:interior bc in direction "
                          << dir << " but not periodic\n";
                Error();
            }
            if (hi_bc[dir] == Interior)
            {
                std::cerr << "Maestro::ReadParameters:interior bc in direction "
                          << dir << " but not periodic\n";
                Error();
            }
        }
    }

    for (int n = 0; n < NUM_STATE; ++n)
    {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            
            // lo-side BCs
            if (lo_bc[idim] == Interior) {
                bcs[n].setLo(idim, BCType::int_dir);  // periodic uses "internal Dirichlet"
            }
            else if (lo_bc[idim] == Inflow) {
                bcs[n].setLo(idim, BCType::ext_dir);  // external Dirichlet
            }
            else if (lo_bc[idim] == Outflow) {
                bcs[n].setLo(idim, BCType::foextrap); // first-order extrapolation
            }
            else if (lo_bc[idim] == Symmetry) {
                bcs[n].setLo(idim, BCType::reflect_even);
            }
            else if (lo_bc[idim] == SlipWall) {
                bcs[n].setLo(idim, BCType::foextrap); // first-order extrapolation
            }
            else if (lo_bc[idim] == NoSlipWall) {
                bcs[n].setLo(idim, BCType::foextrap); // first-order extrapolation
            }
            else {
                Abort("Invalid lo_bc");
            }

            // hi-side BCSs
            if (hi_bc[idim] == Interior) {
                bcs[n].setHi(idim, BCType::int_dir);  // periodic uses "internal Dirichlet"
            }
            else if (hi_bc[idim] == Inflow) {
                bcs[n].setHi(idim, BCType::ext_dir);  // external Dirichlet
            }
            else if (hi_bc[idim] == Outflow) {
                bcs[n].setHi(idim, BCType::foextrap); // first-order extrapolation
            }
            else if (hi_bc[idim] == Symmetry) {
                bcs[n].setHi(idim, BCType::reflect_even);
            }
            else if (hi_bc[idim] == SlipWall) {
                bcs[n].setHi(idim, BCType::foextrap); // first-order extrapolation
            }
            else if (hi_bc[idim] == NoSlipWall) {
                bcs[n].setHi(idim, BCType::foextrap); // first-order extrapolation
            }
            else {
                Abort("Invalid hi_bc");
            }

        }
    }
}

// set covered coarse cells to be the average of overlying fine cells
void
Maestro::AverageDown ()
{
    for (int lev = finest_level-1; lev >= 0; --lev)
    {
        average_down(*snew[lev+1], *snew[lev],
                     geom[lev+1], geom[lev],
                     0, snew[lev]->nComp(), refRatio(lev));
    }
}

// more flexible version of AverageDown() that lets you average down across multiple levels
void
Maestro::AverageDownTo (int crse_lev)
{
    average_down(*snew[crse_lev+1], *snew[crse_lev],
                 geom[crse_lev+1], geom[crse_lev],
                 0, snew[crse_lev]->nComp(), refRatio(crse_lev));
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
    flux_reg[lev].reset(nullptr);
}



#include <Maestro.H>
#include <maestro_defaults.H>

// constructor - reads in parameters from inputs file
//             - sizes multilevel arrays and data structures
Maestro::Maestro ()
{
    ReadParameters();

    ca_set_maestro_method_params();

    // Geometry on all levels has been defined already.

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    int nlevs_max = maxLevel() + 1;

    istep = 0;

    t_new = 0.0;
    t_old = -1.e100;

    // set this to a large number so change_max doesn't affect the first time step
    dt = 1.e100;

    phi_new.resize(nlevs_max);
    phi_old.resize(nlevs_max);

    const int ncomp = 2;
    bcs.resize(ncomp);
    
    for (int n = 0; n < 1; ++n)
    {

        int bc_lo[AMREX_SPACEDIM];
        int bc_hi[AMREX_SPACEDIM];

        bc_lo[0] = INT_DIR;
        bc_hi[0] = INT_DIR;

        bc_lo[1] = INT_DIR;
        bc_hi[1] = INT_DIR;

        if (AMREX_SPACEDIM == 3)
        {
            bc_lo[2] = INT_DIR;
            bc_hi[2] = INT_DIR;
        }

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            
            // lo-side BCs
            if (bc_lo[idim] == INT_DIR) {
                bcs[n].setLo(idim, BCType::int_dir);  // periodic uses "internal Dirichlet"
            }
            else if (bc_lo[idim] == FOEXTRAP) {
                bcs[n].setLo(idim, BCType::foextrap); // first-order extrapolation
            }
            else if (bc_lo[idim] == EXT_DIR) {
                bcs[n].setLo(idim, BCType::ext_dir);  // external Dirichlet
            }
            else {
                amrex::Abort("Invalid bc_lo");
            }

            // hi-side BCSs
            if (bc_hi[idim] == INT_DIR) {
                bcs[n].setHi(idim, BCType::int_dir);  // periodic uses "internal Dirichlet"
            }
            else if (bc_hi[idim] == FOEXTRAP) {
                bcs[n].setHi(idim, BCType::foextrap); // first-order extrapolation
            }
            else if (bc_hi[idim] == EXT_DIR) {
                bcs[n].setHi(idim, BCType::ext_dir);  // external Dirichlet
            }
            else {
                amrex::Abort("Invalid bc_hi");
            }

        }
    }

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

    {
        ParmParse pp("maestro");

#include <maestro_queries.H>

    }

    {
        ParmParse pp;  // Traditionally, max_step and stop_time do not have prefix.
        pp.query("max_step", max_step);
        pp.query("stop_time", stop_time);
    }

    {
        ParmParse pp("amr"); // Traditionally, these have prefix, amr.

        pp.query("regrid_int", regrid_int);
        pp.query("plot_file", plot_file);
        pp.query("plot_int", plot_int);
    }
}

// set covered coarse cells to be the average of overlying fine cells
void
Maestro::AverageDown ()
{
    for (int lev = finest_level-1; lev >= 0; --lev)
    {
        amrex::average_down(*phi_new[lev+1], *phi_new[lev],
                            geom[lev+1], geom[lev],
                            0, phi_new[lev]->nComp(), refRatio(lev));
    }
}

// more flexible version of AverageDown() that lets you average down across multiple levels
void
Maestro::AverageDownTo (int crse_lev)
{
    amrex::average_down(*phi_new[crse_lev+1], *phi_new[crse_lev],
                        geom[crse_lev+1], geom[crse_lev],
                        0, phi_new[crse_lev]->nComp(), refRatio(crse_lev));
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
    phi_new[lev].reset(nullptr);
    phi_old[lev].reset(nullptr);
    flux_reg[lev].reset(nullptr);
}

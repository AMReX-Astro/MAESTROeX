
#include <Maestro.H>
#include <maestro_defaults.H>

using namespace amrex;

int Maestro::Rho       = -1;
int Maestro::RhoH      = -1;
int Maestro::FirstSpec = -1;
int Maestro::NumSpec   = -1;
int Maestro::Temp      = -1;
int Maestro::Pi        = -1;
int Maestro::NSCAL     = -1;

#if (AMREX_SPACEDIM == 2)
IntVect Maestro::nodal_flag(1,1);
IntVect Maestro::nodal_flag_x(1,0);
IntVect Maestro::nodal_flag_y(0,1);
#elif (AMREX_SPACEDIM == 3)
IntVect Maestro::nodal_flag(1,1,1);
IntVect Maestro::nodal_flag_x(1,0,0);
IntVect Maestro::nodal_flag_y(0,1,0);
IntVect Maestro::nodal_flag_z(0,0,1);
#endif

// constructor - reads in parameters from inputs file
//             - sizes multilevel vectors and data structures
Maestro::Maestro ()
{

    ///////
    // Geometry on all levels has been defined already.
    ///////

    // read in C++ parameters in maestro_queries.H using ParmParse pp("maestro");
    ReadParameters();

    // read in F90 parameters in meth_params.F90 that are defined
    // in _cpp_parameters
    read_method_params();

    // define (Rho, RhoH, etc.)
    // calls network_init
    VariableSetup();

    // define additional module variables in meth_params.F90 that are defined
    // at the top of meth_params.template
    set_method_params(Rho,RhoH,FirstSpec,Temp,Pi,NSCAL,
                      geom[0].ProbLo(),geom[0].ProbHi());

    // 
    Box fineBox = geom[max_level].Domain();
    init_base_state_geometry(max_level+1, 
                             geom[max_level].CellSize(),
                             fineBox.hiVect());

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

    sold    .resize(nlevs_max);
    snew    .resize(nlevs_max);
    uold    .resize(nlevs_max);
    unew    .resize(nlevs_max);
    S_cc_old.resize(nlevs_max);
    S_cc_new.resize(nlevs_max);
    gpi     .resize(nlevs_max);
    dSdt    .resize(nlevs_max);

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
    destroy_base_state_geometry();
}

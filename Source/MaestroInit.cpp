
#include <Maestro.H>

using namespace amrex;

// initializes multilevel data
void
Maestro::InitData ()
{
    const Real time = 0.0;

    // here we need to allocate and fill s0_init and p0_init
    const Vector<Geometry>& geom = Geom();
    const Box& domain = geom[max_level].Domain();
    const int& nr_fine = domain.bigEnd()[AMREX_SPACEDIM-1] + 1;

    s0_init  .resize( (max_level+1)*nr_fine*NSCAL );
    p0_init  .resize( (max_level+1)*nr_fine );
    rho0_old .resize( (max_level+1)*nr_fine );
    rho0_new .resize( (max_level+1)*nr_fine );
    rhoh0_old.resize( (max_level+1)*nr_fine );
    rhoh0_new.resize( (max_level+1)*nr_fine );
    p0_old   .resize( (max_level+1)*nr_fine );
    p0_new   .resize( (max_level+1)*nr_fine );

    // read in model file and fill in s0_init and p0_init for all levels
    init_base_state(s0_init.dataPtr(),p0_init.dataPtr(),max_level+1,nr_fine);

    // calls AmrCore::InitFromScratch(), which calls a MakeNewGrids() function 
    // that repeatedly calls Maestro::MakeNewLevelFromScratch() to build and initialize
    InitFromScratch(time);

    // synchronize levels
    AverageDown(snew);
    AverageDown(unew);

    // now fill in rho0, rhoh0, and p0
    // FIXME



    // free memory in s0_init and p0_init by swapping it
    // with an empty vector that will go out of scope
    Vector<Real> s0_swap, p0_swap;
    std::swap(s0_swap,s0_init);
    std::swap(p0_swap,p0_init);

    if (plot_int > 0) {
        WritePlotFile(0);
    }
}

// read in some parameters from inputs file
void
Maestro::ReadParameters ()
{

    ParmParse pp("maestro");

#include <maestro_queries.H>

    // read in boundary conditions
    lo_bc.resize(AMREX_SPACEDIM);
    hi_bc.resize(AMREX_SPACEDIM);
    pp.getarr("lo_bc",lo_bc,0,AMREX_SPACEDIM);
    pp.getarr("hi_bc",hi_bc,0,AMREX_SPACEDIM);

    // read in tagging criteria
    int n = pp.countval("phierr");
    if (n > 0) {
        pp.getarr("phierr", phierr, 0, n);
    }


}

// define variable mappings (Rho, RhoH, ..., NSCAL, etc.)
void Maestro::VariableSetup ()
{

    int cnt = 0;
    Rho = cnt++;
    RhoH = cnt++;

    FirstSpec = cnt;
    get_num_spec(&NumSpec);
    cnt += NumSpec;

    Temp = cnt++;
    Pi = cnt++;

    NSCAL = cnt;  // NumSpec + 4 (Rho, RhoH, Temp, Pi)

    maestro_network_init();

}

// set up BCRec definitions for BC types
void
Maestro::BCSetup()
{
    bcs_s.resize(NSCAL);
    bcs_u.resize(AMREX_SPACEDIM);
    bcs_S.resize(1);
    bcs_g.resize(AMREX_SPACEDIM);
    bcs_d.resize(1);

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

    // set up boundary conditions for Fillpatch operations
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
    {

        // lo-side bcs
        if (lo_bc[dir] == Interior) {
            // periodic uses "internal Dirichlet"
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setLo(dir, BCType::int_dir);
            }
                bcs_s[Rho ].setLo(dir, BCType::int_dir);  
                bcs_s[RhoH].setLo(dir, BCType::int_dir);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setLo(dir, BCType::int_dir);
            }
                bcs_s[Temp].setLo(dir, BCType::int_dir);
                bcs_s[Pi  ].setLo(dir, BCType::int_dir);
                bcs_S[0   ].setLo(dir, BCType::int_dir);
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_g[comp].setLo(dir, BCType::int_dir);
            }
                bcs_d[0   ].setLo(dir, BCType::int_dir);
        }
        else if (lo_bc[dir] == Inflow) {
            // inflow
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setLo(dir, BCType::ext_dir);
            }
                bcs_s[Rho ].setLo(dir, BCType::ext_dir);  
                bcs_s[RhoH].setLo(dir, BCType::ext_dir);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setLo(dir, BCType::ext_dir);
            }
                bcs_s[Temp].setLo(dir, BCType::ext_dir);
                bcs_s[Pi  ].setLo(dir, BCType::foextrap);
                bcs_S[0   ].setLo(dir, BCType::foextrap);
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_g[comp].setLo(dir, BCType::foextrap);
            }
                bcs_d[0   ].setLo(dir, BCType::foextrap);
        }
        else if (lo_bc[dir] == Outflow) {
            // outflow
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setLo(dir, BCType::foextrap);
            }
                bcs_s[Rho ].setLo(dir, BCType::foextrap);  
                bcs_s[RhoH].setLo(dir, BCType::foextrap);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setLo(dir, BCType::foextrap);
            }
                bcs_s[Temp].setLo(dir, BCType::foextrap);
                bcs_s[Pi  ].setLo(dir, BCType::ext_dir);
                bcs_S[0   ].setLo(dir, BCType::foextrap);
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_g[comp].setLo(dir, BCType::foextrap);
            }
                bcs_d[0   ].setLo(dir, BCType::foextrap);
        }
        else if (lo_bc[dir] == Symmetry) {
            // symmetry
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setLo(dir, BCType::reflect_even);
            }
                bcs_u[dir].setLo(dir, BCType::reflect_odd);
                bcs_s[Rho ].setLo(dir, BCType::reflect_even);  
                bcs_s[RhoH].setLo(dir, BCType::reflect_even);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setLo(dir, BCType::reflect_even);
            }
                bcs_s[Temp].setLo(dir, BCType::reflect_even);
                bcs_s[Pi  ].setLo(dir, BCType::reflect_even);
                bcs_S[0   ].setLo(dir, BCType::reflect_even);
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_g[comp].setLo(dir, BCType::reflect_even);
            }
                bcs_d[0   ].setLo(dir, BCType::reflect_even);
        }
        else if (lo_bc[dir] == SlipWall) {
            // slip wall
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setLo(dir, BCType::foextrap);
            }
                bcs_u[dir].setLo(dir, BCType::ext_dir);
                bcs_s[Rho ].setLo(dir, BCType::foextrap);  
                bcs_s[RhoH].setLo(dir, BCType::foextrap);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setLo(dir, BCType::foextrap);
            }
                bcs_s[Temp].setLo(dir, BCType::foextrap);
                bcs_s[Pi  ].setLo(dir, BCType::foextrap);
                bcs_S[0   ].setLo(dir, BCType::foextrap);
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_g[comp].setLo(dir, BCType::foextrap);
            }
                bcs_d[0   ].setLo(dir, BCType::foextrap);
        }
        else if (lo_bc[dir] == NoSlipWall) {
            // no-slip wall
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setLo(dir, BCType::ext_dir);
            }
                bcs_s[Rho ].setLo(dir, BCType::foextrap);  
                bcs_s[RhoH].setLo(dir, BCType::foextrap);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setLo(dir, BCType::foextrap);
            }
                bcs_s[Temp].setLo(dir, BCType::foextrap);
                bcs_s[Pi  ].setLo(dir, BCType::foextrap);
                bcs_S[0   ].setLo(dir, BCType::foextrap);
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_g[comp].setLo(dir, BCType::foextrap);
            }
                bcs_d[0   ].setLo(dir, BCType::foextrap);
        }
        else {
            Abort("Invalid lo_bc");
        }

        // hi-side bcs
        if (hi_bc[dir] == Interior) {
            // periodic uses "internal Dirichlet"
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setHi(dir, BCType::int_dir);
            }
                bcs_s[Rho ].setHi(dir, BCType::int_dir);  
                bcs_s[RhoH].setHi(dir, BCType::int_dir);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setHi(dir, BCType::int_dir);
            }
                bcs_s[Temp].setHi(dir, BCType::int_dir);
                bcs_s[Pi  ].setHi(dir, BCType::int_dir);
                bcs_S[0   ].setHi(dir, BCType::int_dir);
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_g[comp].setHi(dir, BCType::int_dir);
            }
                bcs_d[0   ].setHi(dir, BCType::int_dir);
        }
        else if (hi_bc[dir] == Inflow) {
            // inflow
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setHi(dir, BCType::ext_dir);
            }
                bcs_s[Rho ].setHi(dir, BCType::ext_dir);  
                bcs_s[RhoH].setHi(dir, BCType::ext_dir);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setHi(dir, BCType::ext_dir);
            }
                bcs_s[Temp].setHi(dir, BCType::ext_dir);
                bcs_s[Pi  ].setHi(dir, BCType::foextrap);
                bcs_S[0   ].setHi(dir, BCType::foextrap);
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_g[comp].setHi(dir, BCType::foextrap);
            }
                bcs_d[0   ].setHi(dir, BCType::foextrap);
        }
        else if (hi_bc[dir] == Outflow) {
            // outflow
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setHi(dir, BCType::foextrap);
            }
                bcs_s[Rho ].setHi(dir, BCType::foextrap);  
                bcs_s[RhoH].setHi(dir, BCType::foextrap);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setHi(dir, BCType::foextrap);
            }
                bcs_s[Temp].setHi(dir, BCType::foextrap);
                bcs_s[Pi  ].setHi(dir, BCType::ext_dir);
                bcs_S[0   ].setHi(dir, BCType::foextrap);
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_g[comp].setHi(dir, BCType::foextrap);
            }
                bcs_d[0   ].setHi(dir, BCType::foextrap);
        }
        else if (hi_bc[dir] == Symmetry) {
            // symmetry
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setHi(dir, BCType::reflect_even);
            }
                bcs_u[dir].setHi(dir, BCType::reflect_odd);
                bcs_s[Rho ].setHi(dir, BCType::reflect_even);  
                bcs_s[RhoH].setHi(dir, BCType::reflect_even);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setHi(dir, BCType::reflect_even);
            }
                bcs_s[Temp].setHi(dir, BCType::reflect_even);
                bcs_s[Pi  ].setHi(dir, BCType::reflect_even);
                bcs_S[0   ].setHi(dir, BCType::reflect_even);
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_g[comp].setHi(dir, BCType::reflect_even);
            }
                bcs_d[0   ].setHi(dir, BCType::reflect_even);
        }
        else if (hi_bc[dir] == SlipWall) {
            // slip wall
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setHi(dir, BCType::foextrap);
            }
                bcs_u[dir].setHi(dir, BCType::ext_dir);
                bcs_s[Rho ].setHi(dir, BCType::foextrap);  
                bcs_s[RhoH].setHi(dir, BCType::foextrap);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setHi(dir, BCType::foextrap);
            }
                bcs_s[Temp].setHi(dir, BCType::foextrap);
                bcs_s[Pi  ].setHi(dir, BCType::foextrap);
                bcs_S[0   ].setHi(dir, BCType::foextrap);
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_g[comp].setHi(dir, BCType::foextrap);
            }
                bcs_d[0   ].setHi(dir, BCType::foextrap);
        }
        else if (hi_bc[dir] == NoSlipWall) {
            // no-slip wall
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setHi(dir, BCType::ext_dir);
            }
                bcs_s[Rho ].setHi(dir, BCType::foextrap);  
                bcs_s[RhoH].setHi(dir, BCType::foextrap);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setHi(dir, BCType::foextrap);
            }
                bcs_s[Temp].setHi(dir, BCType::foextrap);
                bcs_s[Pi  ].setHi(dir, BCType::foextrap);
                bcs_S[0   ].setHi(dir, BCType::foextrap);
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_g[comp].setHi(dir, BCType::foextrap);
            }
                bcs_d[0   ].setHi(dir, BCType::foextrap);
        }
        else {
            Abort("Invalid hi_bc");
        }

    } // end loop over directions

    // set up boundary conditions for linear solves
    // order is xlo, xhi, ylo, yhi, zlo, zhi
    // MGT_BC_INT (Interior)
    // MGT_BC_DIR (Dirichlet)
    // MGT_BC_NEU (Neumann)
    for ( int i = 0; i < BL_SPACEDIM; ++i ) {
        if ( Geom()[0].isPeriodic(i) ) {
            mg_bcs_p[i*2 + 0] = MGT_BC_INT;
            mg_bcs_p[i*2 + 1] = MGT_BC_INT;
            mg_bcs_h[i*2 + 0] = MGT_BC_INT;
            mg_bcs_h[i*2 + 1] = MGT_BC_INT;
        }
        else {
            mg_bcs_p[i*2 + 0] = lo_bc[i]==Outflow ? MGT_BC_DIR : MGT_BC_NEU;
            mg_bcs_p[i*2 + 1] = hi_bc[i]==Outflow ? MGT_BC_DIR : MGT_BC_NEU;
            mg_bcs_h[i*2 + 0] = lo_bc[i]==Inflow  ? MGT_BC_DIR : MGT_BC_NEU;
            mg_bcs_h[i*2 + 1] = hi_bc[i]==Inflow  ? MGT_BC_DIR : MGT_BC_NEU;
        }
    }

    // lo/hi_inflow is needed for linear solves
    for (int dir=0; dir<AMREX_SPACEDIM; dir++) {
        if (lo_bc[dir] == Inflow) {
            lo_inflow[dir] = 1;
        }
        if (hi_bc[dir] == Inflow) {
            hi_inflow[dir] = 1;
        }
    }

}

// During initialization of a simulation, Maestro::InitData() calls 
// AmrCore::InitFromScratch(), which calls 
// a MakeNewGrids() function that repeatedly calls this function to build 
// and initialize finer levels.  This function creates a new fine
// level that did not exist before by interpolating from the coarser level
// overrides the pure virtual function in AmrCore
void Maestro::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
				       const DistributionMapping& dm)
{
    const int ng_s = 0;
    const int ng_u = 0;
    const int ng_S = 0;
    const int ng_g = 0;
    const int ng_d = 0;

    sold[lev]    .reset(new MultiFab(ba, dm,          NSCAL, ng_s));
    snew[lev]    .reset(new MultiFab(ba, dm,          NSCAL, ng_s));
    uold[lev]    .reset(new MultiFab(ba, dm, AMREX_SPACEDIM, ng_u));
    unew[lev]    .reset(new MultiFab(ba, dm, AMREX_SPACEDIM, ng_u));
    S_cc_old[lev].reset(new MultiFab(ba, dm,              1, ng_S));
    S_cc_new[lev].reset(new MultiFab(ba, dm,              1, ng_S));
    gpi[lev]     .reset(new MultiFab(ba, dm, AMREX_SPACEDIM, ng_g));
    dSdt[lev]    .reset(new MultiFab(ba, dm,              1, ng_d));

    t_new = time;
    t_old = time - 1.e200;

    if (lev > 0 && do_reflux) {
        flux_reg_s[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, NSCAL));
        flux_reg_u[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, AMREX_SPACEDIM));
    }

    const Real* dx = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();
    Real cur_time = t_new;

    MultiFab& scal = *snew[lev];
    MultiFab& vel = *unew[lev];

    for (MFIter mfi(scal); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        const int* lo  = box.loVect();
        const int* hi  = box.hiVect();

        initdata(lev, cur_time, ARLIM_3D(lo), ARLIM_3D(hi),
                 BL_TO_FORTRAN_3D(scal[mfi]), 
                 BL_TO_FORTRAN_3D(vel[mfi]), 
                 ZFILL(dx), ZFILL(prob_lo), NSCAL);
    }

    // now we copy data from s0_init and p0_init into s0_new and p0_new





}

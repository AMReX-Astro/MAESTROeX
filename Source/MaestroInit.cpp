

#include <Maestro.H>

using namespace amrex;

// initializes multilevel data
void
Maestro::InitData ()
{
    const Real time = 0.0;
    InitFromScratch(time);
    AverageDown(snew);
    AverageDown(unew);

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

}

// define variable mappings (Rho, RhoH, ..., NSCAL, etc.)
void Maestro::VariableSetup ()
{

    int cnt = 0;
    Rho = cnt++;
    RhoH = cnt++;

    FirstSpec = cnt;
    ca_get_num_spec(&NumSpec);
    cnt += NumSpec;

    Temp = cnt++;
    Pi = cnt++;

    NSCAL = cnt;

    ca_network_init();

}

// set up BCRec definitions for BC types
void
Maestro::BCSetup()
{
    // Get boundary conditions
    ParmParse pp("maestro");
    Array<int> lo_bc(AMREX_SPACEDIM), hi_bc(AMREX_SPACEDIM);
    pp.getarr("lo_bc",lo_bc,0,AMREX_SPACEDIM);
    pp.getarr("hi_bc",hi_bc,0,AMREX_SPACEDIM);
    
    bcs_s.resize(NSCAL);
    bcs_u.resize(AMREX_SPACEDIM);

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

    // velocity and scalars
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {

        // lo-side bcs
        if (lo_bc[idim] == Interior) {
            // periodic uses "internal Dirichlet"
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setLo(idim, BCType::int_dir);
            }
                bcs_s[Rho ].setLo(idim, BCType::int_dir);  
                bcs_s[RhoH].setLo(idim, BCType::int_dir);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setLo(idim, BCType::int_dir);
            }
                bcs_s[Temp].setLo(idim, BCType::int_dir);
                bcs_s[Pi  ].setLo(idim, BCType::int_dir);
        }
        else if (lo_bc[idim] == Inflow) {
            // inflow
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setLo(idim, BCType::ext_dir);
            }
                bcs_s[Rho ].setLo(idim, BCType::ext_dir);  
                bcs_s[RhoH].setLo(idim, BCType::ext_dir);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setLo(idim, BCType::ext_dir);
            }
                bcs_s[Temp].setLo(idim, BCType::ext_dir);
                bcs_s[Pi  ].setLo(idim, BCType::foextrap);
        }
        else if (lo_bc[idim] == Outflow) {
            // outflow
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setLo(idim, BCType::foextrap);
            }
                bcs_s[Rho ].setLo(idim, BCType::foextrap);  
                bcs_s[RhoH].setLo(idim, BCType::foextrap);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setLo(idim, BCType::foextrap);
            }
                bcs_s[Temp].setLo(idim, BCType::foextrap);
                bcs_s[Pi  ].setLo(idim, BCType::ext_dir);
        }
        else if (lo_bc[idim] == Symmetry) {
            // symmetry
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setLo(idim, BCType::reflect_even);
            }
                bcs_u[idim].setLo(idim, BCType::reflect_odd);
                bcs_s[Rho ].setLo(idim, BCType::reflect_even);  
                bcs_s[RhoH].setLo(idim, BCType::reflect_even);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setLo(idim, BCType::reflect_even);
            }
                bcs_s[Temp].setLo(idim, BCType::reflect_even);
                bcs_s[Pi  ].setLo(idim, BCType::reflect_even);
        }
        else if (lo_bc[idim] == SlipWall) {
            // slip wall
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setLo(idim, BCType::foextrap);
            }
                bcs_u[idim].setLo(idim, BCType::ext_dir);
                bcs_s[Rho ].setLo(idim, BCType::foextrap);  
                bcs_s[RhoH].setLo(idim, BCType::foextrap);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setLo(idim, BCType::foextrap);
            }
                bcs_s[Temp].setLo(idim, BCType::foextrap);
                bcs_s[Pi  ].setLo(idim, BCType::foextrap);
        }
        else if (lo_bc[idim] == NoSlipWall) {
            // no-slip wall
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setLo(idim, BCType::ext_dir);
            }
                bcs_s[Rho ].setLo(idim, BCType::foextrap);  
                bcs_s[RhoH].setLo(idim, BCType::foextrap);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setLo(idim, BCType::foextrap);
            }
                bcs_s[Temp].setLo(idim, BCType::foextrap);
                bcs_s[Pi  ].setLo(idim, BCType::foextrap);
        }
        else {
            Abort("Invalid lo_bc");
        }

        // hi-side bcs
        if (hi_bc[idim] == Interior) {
            // periodic uses "internal Dirichlet"
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setHi(idim, BCType::int_dir);
            }
                bcs_s[Rho ].setHi(idim, BCType::int_dir);  
                bcs_s[RhoH].setHi(idim, BCType::int_dir);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setHi(idim, BCType::int_dir);
            }
                bcs_s[Temp].setHi(idim, BCType::int_dir);
                bcs_s[Pi  ].setHi(idim, BCType::int_dir);
        }
        else if (hi_bc[idim] == Inflow) {
            // inflow
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setHi(idim, BCType::ext_dir);
            }
                bcs_s[Rho ].setHi(idim, BCType::ext_dir);  
                bcs_s[RhoH].setHi(idim, BCType::ext_dir);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setHi(idim, BCType::ext_dir);
            }
                bcs_s[Temp].setHi(idim, BCType::ext_dir);
                bcs_s[Pi  ].setHi(idim, BCType::foextrap);
        }
        else if (hi_bc[idim] == Outflow) {
            // outflow
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setHi(idim, BCType::foextrap);
            }
                bcs_s[Rho ].setHi(idim, BCType::foextrap);  
                bcs_s[RhoH].setHi(idim, BCType::foextrap);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setHi(idim, BCType::foextrap);
            }
                bcs_s[Temp].setHi(idim, BCType::foextrap);
                bcs_s[Pi  ].setHi(idim, BCType::ext_dir);
        }
        else if (hi_bc[idim] == Symmetry) {
            // symmetry
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setHi(idim, BCType::reflect_even);
            }
                bcs_u[idim].setHi(idim, BCType::reflect_odd);
                bcs_s[Rho ].setHi(idim, BCType::reflect_even);  
                bcs_s[RhoH].setHi(idim, BCType::reflect_even);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setHi(idim, BCType::reflect_even);
            }
                bcs_s[Temp].setHi(idim, BCType::reflect_even);
                bcs_s[Pi  ].setHi(idim, BCType::reflect_even);
        }
        else if (hi_bc[idim] == SlipWall) {
            // slip wall
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setHi(idim, BCType::foextrap);
            }
                bcs_u[idim].setHi(idim, BCType::ext_dir);
                bcs_s[Rho ].setHi(idim, BCType::foextrap);  
                bcs_s[RhoH].setHi(idim, BCType::foextrap);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setHi(idim, BCType::foextrap);
            }
                bcs_s[Temp].setHi(idim, BCType::foextrap);
                bcs_s[Pi  ].setHi(idim, BCType::foextrap);
        }
        else if (hi_bc[idim] == NoSlipWall) {
            // no-slip wall
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setHi(idim, BCType::ext_dir);
            }
                bcs_s[Rho ].setHi(idim, BCType::foextrap);  
                bcs_s[RhoH].setHi(idim, BCType::foextrap);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setHi(idim, BCType::foextrap);
            }
                bcs_s[Temp].setHi(idim, BCType::foextrap);
                bcs_s[Pi  ].setHi(idim, BCType::foextrap);
        }
        else {
            Abort("Invalid hi_bc");
        }

    } // end loop over directions
}

// Make a new level from scratch using provided BoxArray and DistributionMapping.
// Only used during initialization.
// overrides the pure virtual function in AmrCore
void Maestro::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
				       const DistributionMapping& dm)
{
    const int nghost = 0;

    snew[lev].reset(new MultiFab(ba, dm, NSCAL     , nghost));
    sold[lev].reset(new MultiFab(ba, dm, NSCAL     , nghost));
    unew[lev].reset(new MultiFab(ba, dm, AMREX_SPACEDIM, nghost));
    uold[lev].reset(new MultiFab(ba, dm, AMREX_SPACEDIM, nghost));

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
}

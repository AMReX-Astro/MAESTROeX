
#include <Maestro.H>
using namespace amrex;

// read in C++/F90 parameters
// define global C++/F90 variables and initialize network
// set up boundary conditions
// initialize base state geometry parameters
// set istep, t_new, t_old
// allocate MultiFabs and base state arrays
void
Maestro::Setup ()
{
    Print() << "Calling Setup()" << endl;

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
    set_method_params(Rho,RhoH,FirstSpec,Temp,Pi,Nscal);

    // set up BCRec definitions for BC types
    BCSetup();

    const Box domainBoxFine = geom[max_level].Domain();
    const Real* dxFine = geom[max_level].CellSize();
    const Real* probLo = geom[0].ProbLo();
    const Real* probHi = geom[0].ProbHi();

    if (spherical == 1) {

        // compute max_radial_level
        max_radial_level = 0;

        // compute dr_fine
        dr_fine = dxFine[0] / drdxfac;
       
        // compute nr_fine
        double lenx, leny, lenz, max_dist;
        if (octant) {
            lenx = probHi[0] - probLo[0];
            leny = probHi[1] - probLo[1];
            lenz = probHi[2] - probLo[2];
        }
        else {
            lenx = 0.5*(probHi[0] - probLo[0]);
            leny = 0.5*(probHi[1] - probLo[1]);
            lenz = 0.5*(probHi[2] - probLo[2]);
        }       
        max_dist = sqrt(lenx*lenx + leny*leny + lenz*lenz);
        nr_fine = int(max_dist / dr_fine) + 1;
    }
    else {
        // compute max_radial_level
        max_radial_level = max_level;

        // compute dr_fine
        dr_fine = dxFine[AMREX_SPACEDIM-1];

        // compute nr_fine
        nr_fine = domainBoxFine.bigEnd()[AMREX_SPACEDIM-1] + 1;
    }

    // vectors store the multilevel 1D states as one very long array
    // these are cell-centered
    s0_init       .resize( (max_radial_level+1)*nr_fine*Nscal );
    p0_init       .resize( (max_radial_level+1)*nr_fine );
    rho0_old      .resize( (max_radial_level+1)*nr_fine );
    rho0_new      .resize( (max_radial_level+1)*nr_fine );
    rhoh0_old     .resize( (max_radial_level+1)*nr_fine );
    rhoh0_new     .resize( (max_radial_level+1)*nr_fine );
    p0_old        .resize( (max_radial_level+1)*nr_fine );
    p0_new        .resize( (max_radial_level+1)*nr_fine );
    tempbar       .resize( (max_radial_level+1)*nr_fine );
    tempbar_init  .resize( (max_radial_level+1)*nr_fine );
    beta0_old     .resize( (max_radial_level+1)*nr_fine );
    beta0_new     .resize( (max_radial_level+1)*nr_fine );
    gamma1bar_new .resize( (max_radial_level+1)*nr_fine );
    gamma1bar_init.resize( (max_radial_level+1)*nr_fine );
    etarho_cc     .resize( (max_radial_level+1)*nr_fine );
    psi           .resize( (max_radial_level+1)*nr_fine );
    grav_cell     .resize( (max_radial_level+1)*nr_fine );
    r_cc_loc      .resize( (max_radial_level+1)*nr_fine );

    // vectors store the multilevel 1D states as one very long array
    // these are edge-centered
    w0        .resize( (max_radial_level+1)*(nr_fine+1) );
    etarho_ec .resize( (max_radial_level+1)*(nr_fine+1) );
    r_edge_loc.resize( (max_radial_level+1)*(nr_fine+1) );

    init_base_state_geometry(max_radial_level,nr_fine,dr_fine,
                             r_cc_loc.dataPtr(),
                             r_edge_loc.dataPtr(),
                             geom[max_level].CellSize(),
                             domainBoxFine.hiVect(),
                             geom[0].ProbLo(),
                             geom[0].ProbHi());

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    istep = 0;

    t_new = 0.0;
    t_old = -1.e100;

    // set this to a large number so change_max doesn't affect the first time step
    dt = 1.e100;
    dtold = 1.e100;

    sold    .resize(max_level+1);
    snew    .resize(max_level+1);
    uold    .resize(max_level+1);
    unew    .resize(max_level+1);
    S_cc_old.resize(max_level+1);
    S_cc_new.resize(max_level+1);
    gpi     .resize(max_level+1);
    dSdt    .resize(max_level+1);

    // stores fluxes at coarse-fine interface for synchronization
    // this will be sized "max_level+1+1"
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
    // therefore flux_reg[0] is never actually used in the reflux operation
    flux_reg_s.resize(max_level+1);
    flux_reg_u.resize(max_level+1);

}

// read in some parameters from inputs file
void
Maestro::ReadParameters ()
{
    Print() << "Calling ReadParameters()" << endl;

    ParmParse pp("maestro");

#include <maestro_queries.H>

    // read in boundary conditions
    lo_bc.resize(AMREX_SPACEDIM);
    hi_bc.resize(AMREX_SPACEDIM);
    pp.getarr("lo_bc",lo_bc,0,AMREX_SPACEDIM);
    pp.getarr("hi_bc",hi_bc,0,AMREX_SPACEDIM);

    // read in tagging criteria
    int n = pp.countval("temperr");
    if (n > 0) {
        pp.getarr("temperr", temperr, 0, n);
    }

}

// define variable mappings (Rho, RhoH, ..., Nscal, etc.)
void Maestro::VariableSetup ()
{
    Print() << "Calling VariableSetup()" << endl;

    int cnt = 0;
    Rho = cnt++;
    RhoH = cnt++;

    FirstSpec = cnt;
    get_num_spec(&NumSpec);
    cnt += NumSpec;

    Temp = cnt++;
    Pi = cnt++;

    Nscal = cnt;  // NumSpec + 4 (Rho, RhoH, Temp, Pi)

    maestro_network_init();
}

// set up BCRec definitions for BC types
void
Maestro::BCSetup()
{
    Print() << "Calling BCSetup()" << endl;

    bcs_s.resize(Nscal);          // scalars
    bcs_u.resize(AMREX_SPACEDIM); // velocitiy
    bcs_h.resize(1);              // rhoh
    bcs_t.resize(1);              // temperature
    bcs_f.resize(AMREX_SPACEDIM); // a vector of "first-order extrap"

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
                bcs_h[0   ].setLo(dir, BCType::int_dir);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setLo(dir, BCType::int_dir);
            }
                bcs_s[Temp].setLo(dir, BCType::int_dir);
                bcs_t[0   ].setLo(dir, BCType::int_dir);
                bcs_s[Pi  ].setLo(dir, BCType::int_dir);
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_f[comp].setLo(dir, BCType::int_dir);
            }
        }
        else if (lo_bc[dir] == Inflow) {
            // inflow
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setLo(dir, BCType::ext_dir);
            }
                bcs_s[Rho ].setLo(dir, BCType::ext_dir);  
                bcs_s[RhoH].setLo(dir, BCType::ext_dir);
                bcs_h[0   ].setLo(dir, BCType::ext_dir);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setLo(dir, BCType::ext_dir);
            }
                bcs_s[Temp].setLo(dir, BCType::ext_dir);
                bcs_t[0   ].setLo(dir, BCType::ext_dir);
                bcs_s[Pi  ].setLo(dir, BCType::foextrap);
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_f[comp].setLo(dir, BCType::foextrap);
            }
        }
        else if (lo_bc[dir] == Outflow) {
            // outflow
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setLo(dir, BCType::foextrap);
            }
                bcs_s[Rho ].setLo(dir, BCType::foextrap);  
                bcs_s[RhoH].setLo(dir, BCType::foextrap);
                bcs_h[0   ].setLo(dir, BCType::foextrap);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setLo(dir, BCType::foextrap);
            }
                bcs_s[Temp].setLo(dir, BCType::foextrap);
                bcs_t[0   ].setLo(dir, BCType::foextrap);
                bcs_s[Pi  ].setLo(dir, BCType::ext_dir);
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_f[comp].setLo(dir, BCType::foextrap);
            }
        }
        else if (lo_bc[dir] == Symmetry) {
            // symmetry
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setLo(dir, BCType::reflect_even);
            }
                bcs_u[dir ].setLo(dir, BCType::reflect_odd);
                bcs_s[Rho ].setLo(dir, BCType::reflect_even);  
                bcs_s[RhoH].setLo(dir, BCType::reflect_even);
                bcs_h[0   ].setLo(dir, BCType::reflect_even);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setLo(dir, BCType::reflect_even);
            }
                bcs_s[Temp].setLo(dir, BCType::reflect_even);
                bcs_t[0   ].setLo(dir, BCType::reflect_even);
                bcs_s[Pi  ].setLo(dir, BCType::reflect_even);
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_f[comp].setLo(dir, BCType::reflect_even);
            }
        }
        else if (lo_bc[dir] == SlipWall) {
            // slip wall
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setLo(dir, BCType::foextrap);
            }
                bcs_u[dir ].setLo(dir, BCType::ext_dir);
                bcs_s[Rho ].setLo(dir, BCType::foextrap);  
                bcs_s[RhoH].setLo(dir, BCType::foextrap);
                bcs_h[0   ].setLo(dir, BCType::foextrap);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setLo(dir, BCType::foextrap);
            }
                bcs_s[Temp].setLo(dir, BCType::foextrap);
                bcs_t[0   ].setLo(dir, BCType::foextrap);
                bcs_s[Pi  ].setLo(dir, BCType::foextrap);
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_f[comp].setLo(dir, BCType::foextrap);
            }
        }
        else if (lo_bc[dir] == NoSlipWall) {
            // no-slip wall
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setLo(dir, BCType::ext_dir);
            }
                bcs_s[Rho ].setLo(dir, BCType::foextrap);  
                bcs_s[RhoH].setLo(dir, BCType::foextrap);
                bcs_h[0   ].setLo(dir, BCType::foextrap);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setLo(dir, BCType::foextrap);
            }
                bcs_s[Temp].setLo(dir, BCType::foextrap);
                bcs_t[0   ].setLo(dir, BCType::foextrap);
                bcs_s[Pi  ].setLo(dir, BCType::foextrap);
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_f[comp].setLo(dir, BCType::foextrap);
            }
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
                bcs_h[0   ].setHi(dir, BCType::int_dir);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setHi(dir, BCType::int_dir);
            }
                bcs_s[Temp].setHi(dir, BCType::int_dir);
                bcs_t[0   ].setHi(dir, BCType::int_dir);
                bcs_s[Pi  ].setHi(dir, BCType::int_dir);
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_f[comp].setHi(dir, BCType::int_dir);
            }
        }
        else if (hi_bc[dir] == Inflow) {
            // inflow
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setHi(dir, BCType::ext_dir);
            }
                bcs_s[Rho ].setHi(dir, BCType::ext_dir);  
                bcs_s[RhoH].setHi(dir, BCType::ext_dir);
                bcs_h[0   ].setHi(dir, BCType::ext_dir);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setHi(dir, BCType::ext_dir);
            }
                bcs_s[Temp].setHi(dir, BCType::ext_dir);
                bcs_t[0   ].setHi(dir, BCType::ext_dir);
                bcs_s[Pi  ].setHi(dir, BCType::foextrap);
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_f[comp].setHi(dir, BCType::foextrap);
            }
        }
        else if (hi_bc[dir] == Outflow) {
            // outflow
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setHi(dir, BCType::foextrap);
            }
                bcs_s[Rho ].setHi(dir, BCType::foextrap);  
                bcs_s[RhoH].setHi(dir, BCType::foextrap);
                bcs_h[0   ].setHi(dir, BCType::foextrap);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setHi(dir, BCType::foextrap);
            }
                bcs_s[Temp].setHi(dir, BCType::foextrap);
                bcs_t[0   ].setHi(dir, BCType::foextrap);
                bcs_s[Pi  ].setHi(dir, BCType::ext_dir);
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_f[comp].setHi(dir, BCType::foextrap);
            }
        }
        else if (hi_bc[dir] == Symmetry) {
            // symmetry
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setHi(dir, BCType::reflect_even);
            }
                bcs_u[dir].setHi(dir, BCType::reflect_odd);
                bcs_s[Rho ].setHi(dir, BCType::reflect_even);  
                bcs_s[RhoH].setHi(dir, BCType::reflect_even);
                bcs_h[0   ].setHi(dir, BCType::reflect_even);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setHi(dir, BCType::reflect_even);
            }
                bcs_s[Temp].setHi(dir, BCType::reflect_even);
                bcs_t[0   ].setHi(dir, BCType::reflect_even);
                bcs_s[Pi  ].setHi(dir, BCType::reflect_even);
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_f[comp].setHi(dir, BCType::reflect_even);
            }
        }
        else if (hi_bc[dir] == SlipWall) {
            // slip wall
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setHi(dir, BCType::foextrap);
            }
                bcs_u[dir].setHi(dir, BCType::ext_dir);
                bcs_s[Rho ].setHi(dir, BCType::foextrap);  
                bcs_s[RhoH].setHi(dir, BCType::foextrap);
                bcs_h[0   ].setHi(dir, BCType::foextrap);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setHi(dir, BCType::foextrap);
            }
                bcs_s[Temp].setHi(dir, BCType::foextrap);
                bcs_t[0   ].setHi(dir, BCType::foextrap);
                bcs_s[Pi  ].setHi(dir, BCType::foextrap);
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_f[comp].setHi(dir, BCType::foextrap);
            }
        }
        else if (hi_bc[dir] == NoSlipWall) {
            // no-slip wall
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setHi(dir, BCType::ext_dir);
            }
                bcs_s[Rho ].setHi(dir, BCType::foextrap);  
                bcs_s[RhoH].setHi(dir, BCType::foextrap);
                bcs_h[0   ].setHi(dir, BCType::foextrap);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setHi(dir, BCType::foextrap);
            }
                bcs_s[Temp].setHi(dir, BCType::foextrap);
                bcs_t[0   ].setHi(dir, BCType::foextrap);
                bcs_s[Pi  ].setHi(dir, BCType::foextrap);
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_f[comp].setHi(dir, BCType::foextrap);
            }
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
    for ( int i = 0; i < AMREX_SPACEDIM; ++i ) {
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

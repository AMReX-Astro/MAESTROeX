
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
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Setup()",Setup);

    Print() << "Calling Setup()" << std::endl;

    // initialize the start time for our CPU-time tracker
    startCPUTime = ParallelDescriptor::second();

    ///////
    // Geometry on all levels has been defined already.
    ///////

    // read in C++ parameters in maestro_queries.H using ParmParse pp("maestro");
    ReadParameters();

    // read in F90 parameters in meth_params.F90 that are defined
    // in _cpp_parameters
    read_method_params();

    // Initialize the runtime parameters for any of the external microphysics
    // (in extern.f90)
    ExternInit();

    // define (Rho, RhoH, etc.)
    // calls network_init
    VariableSetup();

    maestro_network_init();

    burner_init();

    maestro_conductivity_init();

    const Real* probLo = geom[0].ProbLo();
    const Real* probHi = geom[0].ProbHi();

    // define additional module variables in meth_params.F90 that are defined
    // at the top of meth_params.template
    set_method_params(&Rho,&RhoH,&FirstSpec,&Temp,&Pi,&Nscal,
                      ZFILL(probLo),ZFILL(probHi));

    // set up BCRec definitions for BC types
    BCSetup();

    const Box domainBoxFine = geom[max_level].Domain();
    const Real* dxFine = geom[max_level].CellSize();

    if (spherical == 1) {

        // compute max_radial_level
        max_radial_level = 0;

        // compute dr_fine
        dr_fine = dxFine[0] / drdxfac;

	// compute nr_irreg
	int domhi = domainBoxFine.bigEnd(0)+1;
	if (!octant) {
	    nr_irreg = (3*(domhi/2-0.5)*(domhi/2-0.5)-0.75)/2.0;
	} else {
	    nr_irreg = (3*(domhi-0.5)*(domhi-0.5)-0.75)/2.0;
	}

        // compute nr_fine
	if (use_exact_base_state) {
	    nr_fine = nr_irreg + 1;
	} else {
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
    // base states are stored
    s0_init      .resize( (max_radial_level+1)*nr_fine*Nscal );
    p0_init      .resize( (max_radial_level+1)*nr_fine );
    rho0_old     .resize( (max_radial_level+1)*nr_fine );
    rho0_new     .resize( (max_radial_level+1)*nr_fine );
    rhoh0_old    .resize( (max_radial_level+1)*nr_fine );
    rhoh0_new    .resize( (max_radial_level+1)*nr_fine );
    p0_old       .resize( (max_radial_level+1)*nr_fine );
    p0_new       .resize( (max_radial_level+1)*nr_fine );
    p0_nm1       .resize( (max_radial_level+1)*nr_fine );
    tempbar      .resize( (max_radial_level+1)*nr_fine );
    tempbar_init .resize( (max_radial_level+1)*nr_fine );
    beta0_old    .resize( (max_radial_level+1)*nr_fine );
    beta0_new    .resize( (max_radial_level+1)*nr_fine );
    beta0_nm1    .resize( (max_radial_level+1)*nr_fine );
    gamma1bar_old.resize( (max_radial_level+1)*nr_fine );
    gamma1bar_new.resize( (max_radial_level+1)*nr_fine );
    grav_cell_old.resize( (max_radial_level+1)*nr_fine );
    grav_cell_new.resize( (max_radial_level+1)*nr_fine );
    r_cc_loc     .resize( (max_radial_level+1)*nr_fine );
    etarho_cc    .resize( (max_radial_level+1)*nr_fine );
    psi          .resize( (max_radial_level+1)*nr_fine );

    // vectors store the multilevel 1D states as one very long array
    // these are edge-centered
    r_edge_loc.resize( (max_radial_level+1)*(nr_fine+1) );
    w0        .resize( (max_radial_level+1)*(nr_fine+1) );
    etarho_ec .resize( (max_radial_level+1)*(nr_fine+1) );

    // tagged box array for multilevel (planar)
    tag_array .resize( (max_radial_level+1)*nr_fine );

    // diag file data arrays
    diagfile1_data.resize(diag_buf_size*11);
    diagfile2_data.resize(diag_buf_size*11);
    diagfile3_data.resize(diag_buf_size*10);

    // make sure C++ is as efficient as possible with memory usage
    s0_init      .shrink_to_fit();
    p0_init      .shrink_to_fit();
    rho0_old     .shrink_to_fit();
    rho0_new     .shrink_to_fit();
    rhoh0_old    .shrink_to_fit();
    rhoh0_new    .shrink_to_fit();
    p0_old       .shrink_to_fit();
    p0_new       .shrink_to_fit();
    p0_nm1       .shrink_to_fit();
    tempbar      .shrink_to_fit();
    tempbar_init .shrink_to_fit();
    beta0_old    .shrink_to_fit();
    beta0_new    .shrink_to_fit();
    beta0_nm1    .shrink_to_fit();
    gamma1bar_old.shrink_to_fit();
    gamma1bar_new.shrink_to_fit();
    grav_cell_old.shrink_to_fit();
    grav_cell_new.shrink_to_fit();
    w0           .shrink_to_fit();
    etarho_cc    .shrink_to_fit();
    psi          .shrink_to_fit();
    etarho_ec    .shrink_to_fit();
    r_cc_loc     .shrink_to_fit();
    r_edge_loc   .shrink_to_fit();
    tag_array    .shrink_to_fit();
    diagfile1_data.shrink_to_fit();
    diagfile2_data.shrink_to_fit();
    diagfile3_data.shrink_to_fit();

    init_base_state_geometry(&max_radial_level,&nr_fine,&dr_fine,
			     r_cc_loc.dataPtr(),
			     r_edge_loc.dataPtr(),
			     geom[max_level].CellSize(),
			     &nr_irreg);

    if (use_exact_base_state) average_base_state = 1;

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    istep = 0;

    t_new = 1.e99;
    t_old = 0.0;

    // set this to a large number so change_max doesn't affect the first time step
    dt = 1.e100;
    dtold = 1.e100;

    sold              .resize(max_level+1);
    snew              .resize(max_level+1);
    uold              .resize(max_level+1);
    unew              .resize(max_level+1);
    S_cc_old          .resize(max_level+1);
    S_cc_new          .resize(max_level+1);
    gpi               .resize(max_level+1);
    dSdt              .resize(max_level+1);
    pi                .resize(max_level+1);
    rhcc_for_nodalproj.resize(max_level+1);
    normal            .resize(max_level+1);
    cell_cc_to_r      .resize(max_level+1);

    // stores fluxes at coarse-fine interface for synchronization
    // this will be sized "max_level+2"
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
    // therefore flux_reg[0] is never actually used in the reflux operation
    flux_reg_s.resize(max_level+2);

    // number of ghost cells needed for hyperbolic step
    if (ppm_type == 2 || bds_type == 1) {
        ng_adv = 4;
    }
    else {
        ng_adv = 3;
    }

    std::fill(tag_array.begin(), tag_array.end(), 0);

    // if do_smallscale check other parameters for consistency

    if (do_smallscale && (beta0_type != 3 || evolve_base_state)) {
      std::cerr << "Error: do_smallscale = T requires beta0_type = 3 and evolve_base_state = F" << std::endl;
      std::cerr << "    do_smallscale = " << do_smallscale << std::endl;
      std::cerr << "    beta0_type = " << beta0_type << std::endl;
      std::cerr << "    evolve_base_state = " << evolve_base_state << std::endl;
      Error();
    }

}

// read in some parameters from inputs file
void
Maestro::ReadParameters ()
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::ReadParameters()",ReadParameters);

    Print() << "Calling ReadParameters()" << std::endl;

    ParmParse pp("maestro");

#include <maestro_queries.H>

    // now read in vectors for ParmParse

    // read in boundary conditions
    Vector<int> lo_bc(AMREX_SPACEDIM);
    Vector<int> hi_bc(AMREX_SPACEDIM);
    pp.getarr("lo_bc",lo_bc,0,AMREX_SPACEDIM);
    pp.getarr("hi_bc",hi_bc,0,AMREX_SPACEDIM);

    // storey boundary conditions in a single array
    // order shall be
    // LO_X, LO_Y, (LO_Z), HI_X, HI_Y, (HI_Z)
    phys_bc.resize(2*AMREX_SPACEDIM);
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        phys_bc[i]                = lo_bc[i];
        phys_bc[i+AMREX_SPACEDIM] = hi_bc[i];
    }
}

// define variable mappings (Rho, RhoH, ..., Nscal, etc.)
void Maestro::VariableSetup ()
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::VariableSetup()",VariableSetup);

    Print() << "Calling VariableSetup()" << std::endl;

    int cnt = 0;
    Rho = cnt++;
    RhoH = cnt++;

    FirstSpec = cnt;
    get_num_spec(&NumSpec);
    cnt += NumSpec;

    Temp = cnt++;
    Pi = cnt++;

    Nscal = cnt;  // NumSpec + 4 (Rho, RhoH, Temp, Pi)

    // set number of ghost cells for sold/new and uold/new
    ng_s = 3;
    if (ppm_type == 2 || bds_type == 1) {
        ng_s = 4;
    }
}

void
Maestro::ExternInit ()
{
  // initialize the external runtime parameters -- these will
  // live in the probin

  if (ParallelDescriptor::IOProcessor()) {
    std::cout << "reading extern runtime parameters ..." << std::endl;
  }

  maestro_extern_init();
}

// set up BCRec definitions for BC types
void
Maestro::BCSetup()
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::BCSetup()",BCSetup);

    Print() << "Calling BCSetup()" << std::endl;

    bcs_s.resize(Nscal);          // scalars
    bcs_u.resize(AMREX_SPACEDIM); // velocitiy
    bcs_f.resize(Nscal);          // a vector of "first-order extrap"

    // Check phys_bc against possible periodic geometry
    // if periodic, must have internal BC marked.
    //
    if (Geom(0).isAnyPeriodic())
    {
        //
        // Do idiot check.  Periodic means interior in those directions.
        //
        for (int dir = 0; dir<AMREX_SPACEDIM; dir++)
        {
            if (Geom(0).isPeriodic(dir))
            {
                if (phys_bc[dir] != Interior)
                {
                    std::cerr << "Maestro::ReadParameters:periodic in direction "
                              << dir << " but low BC is not Interior\n";
                    Error();
                }
                if (phys_bc[AMREX_SPACEDIM+dir] != Interior)
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
            if (phys_bc[dir] == Interior)
            {
                std::cerr << "Maestro::ReadParameters:interior bc in direction "
                          << dir << " but not periodic\n";
                Error();
            }
            if (phys_bc[AMREX_SPACEDIM+dir] == Interior)
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
        if (phys_bc[dir] == Interior) {
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
            for (int comp=0; comp<Nscal; ++comp) {
                bcs_f[comp].setLo(dir, BCType::int_dir);
            }
        }
        else if (phys_bc[dir] == Inflow) {
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
            for (int comp=0; comp<Nscal; ++comp) {
                bcs_f[comp].setLo(dir, BCType::foextrap);
            }
        }
        else if (phys_bc[dir] == Outflow) {
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
            for (int comp=0; comp<Nscal; ++comp) {
                bcs_f[comp].setLo(dir, BCType::foextrap);
            }
        }
        else if (phys_bc[dir] == Symmetry) {
            // symmetry
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setLo(dir, BCType::reflect_even);
            }
                bcs_u[dir ].setLo(dir, BCType::reflect_odd);
                bcs_s[Rho ].setLo(dir, BCType::reflect_even);
                bcs_s[RhoH].setLo(dir, BCType::reflect_even);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setLo(dir, BCType::reflect_even);
            }
                bcs_s[Temp].setLo(dir, BCType::reflect_even);
                bcs_s[Pi  ].setLo(dir, BCType::reflect_even);
            for (int comp=0; comp<Nscal; ++comp) {
                bcs_f[comp].setLo(dir, BCType::reflect_even);
            }
        }
        else if (phys_bc[dir] == SlipWall) {
            // slip wall
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setLo(dir, BCType::hoextrap);
            }
                bcs_u[dir ].setLo(dir, BCType::ext_dir);
                bcs_s[Rho ].setLo(dir, BCType::hoextrap);
                bcs_s[RhoH].setLo(dir, BCType::hoextrap);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setLo(dir, BCType::hoextrap);
            }
                bcs_s[Temp].setLo(dir, BCType::hoextrap);
                bcs_s[Pi  ].setLo(dir, BCType::foextrap);
            for (int comp=0; comp<Nscal; ++comp) {
                bcs_f[comp].setLo(dir, BCType::foextrap);
            }
        }
        else if (phys_bc[dir] == NoSlipWall) {
            // no-slip wall
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setLo(dir, BCType::ext_dir);
            }
                bcs_s[Rho ].setLo(dir, BCType::hoextrap);
                bcs_s[RhoH].setLo(dir, BCType::hoextrap);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setLo(dir, BCType::hoextrap);
            }
                bcs_s[Temp].setLo(dir, BCType::hoextrap);
                bcs_s[Pi  ].setLo(dir, BCType::foextrap);
            for (int comp=0; comp<Nscal; ++comp) {
                bcs_f[comp].setLo(dir, BCType::foextrap);
            }
        }
        else {
            Abort("Invalid lo_bc");
        }

        // hi-side bcs
        if (phys_bc[AMREX_SPACEDIM+dir] == Interior) {
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
            for (int comp=0; comp<Nscal; ++comp) {
                bcs_f[comp].setHi(dir, BCType::int_dir);
            }
        }
        else if (phys_bc[AMREX_SPACEDIM+dir] == Inflow) {
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
            for (int comp=0; comp<Nscal; ++comp) {
                bcs_f[comp].setHi(dir, BCType::foextrap);
            }
        }
        else if (phys_bc[AMREX_SPACEDIM+dir] == Outflow) {
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
            for (int comp=0; comp<Nscal; ++comp) {
                bcs_f[comp].setHi(dir, BCType::foextrap);
            }
        }
        else if (phys_bc[AMREX_SPACEDIM+dir] == Symmetry) {
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
            for (int comp=0; comp<Nscal; ++comp) {
                bcs_f[comp].setHi(dir, BCType::reflect_even);
            }
        }
        else if (phys_bc[AMREX_SPACEDIM+dir] == SlipWall) {
            // slip wall
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setHi(dir, BCType::hoextrap);
            }
                bcs_u[dir].setHi(dir, BCType::ext_dir);
                bcs_s[Rho ].setHi(dir, BCType::hoextrap);
                bcs_s[RhoH].setHi(dir, BCType::hoextrap);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setHi(dir, BCType::hoextrap);
            }
                bcs_s[Temp].setHi(dir, BCType::hoextrap);
                bcs_s[Pi  ].setHi(dir, BCType::foextrap);
            for (int comp=0; comp<Nscal; ++comp) {
                bcs_f[comp].setHi(dir, BCType::foextrap);
            }
        }
        else if (phys_bc[AMREX_SPACEDIM+dir] == NoSlipWall) {
            // no-slip wall
            for (int comp=0; comp<AMREX_SPACEDIM; ++comp) {
                bcs_u[comp].setHi(dir, BCType::ext_dir);
            }
                bcs_s[Rho ].setHi(dir, BCType::hoextrap);
                bcs_s[RhoH].setHi(dir, BCType::hoextrap);
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                bcs_s[comp].setHi(dir, BCType::hoextrap);
            }
                bcs_s[Temp].setHi(dir, BCType::hoextrap);
                bcs_s[Pi  ].setHi(dir, BCType::foextrap);
            for (int comp=0; comp<Nscal; ++comp) {
                bcs_f[comp].setHi(dir, BCType::foextrap);
            }
        }
        else {
            Abort("Invalid hi_bc");
        }

    } // end loop over directions

}

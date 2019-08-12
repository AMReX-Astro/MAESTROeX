
#include <Maestro.H>

using namespace amrex;

// advance a single level for a single time step, updates flux registers
void
Maestro::AdvanceTimeStepSDC (bool is_initIter) {

    // timer for profiling
    BL_PROFILE_VAR("Maestro::AdvanceTimeStepSDC()",AdvanceTimeStepSDC);

    // cell-centered MultiFabs needed within the AdvanceTimeStep routine
    Vector<MultiFab>           shat(finest_level+1);
    Vector<MultiFab>        rhohalf(finest_level+1);
    Vector<MultiFab>         cphalf(finest_level+1);
    Vector<MultiFab>         xihalf(finest_level+1);
    Vector<MultiFab>         macrhs(finest_level+1);
    Vector<MultiFab>         macphi(finest_level+1);
    Vector<MultiFab>       S_cc_nph(finest_level+1);
    Vector<MultiFab>   rho_omegadot(finest_level+1);
    Vector<MultiFab>       diff_old(finest_level+1);
    Vector<MultiFab>       diff_hat(finest_level+1);
    Vector<MultiFab> diff_hterm_old(finest_level+1);
    Vector<MultiFab> diff_hterm_hat(finest_level+1);
    Vector<MultiFab>       rho_Hnuc(finest_level+1);
    Vector<MultiFab>       rho_Hext(finest_level+1);
    Vector<MultiFab>             s1(finest_level+1);
    Vector<MultiFab>             s2(finest_level+1);
    Vector<MultiFab>          intra(finest_level+1);
    Vector<MultiFab>     sdc_source(finest_level+1);
    Vector<MultiFab>           aofs(finest_level+1);
    
    Vector<MultiFab> delta_gamma1_term(finest_level+1);
    Vector<MultiFab>      delta_gamma1(finest_level+1);
    Vector<MultiFab>      peosbar_cart(finest_level+1);
    Vector<MultiFab>           p0_cart(finest_level+1);
    Vector<MultiFab>      delta_p_term(finest_level+1);
    
    Vector<MultiFab>      Tcoeff1(finest_level+1);
    Vector<MultiFab>      hcoeff1(finest_level+1);
    Vector<MultiFab>     Xkcoeff1(finest_level+1);
    Vector<MultiFab>      pcoeff1(finest_level+1);
    Vector<MultiFab>      Tcoeff2(finest_level+1);
    Vector<MultiFab>      hcoeff2(finest_level+1);
    Vector<MultiFab>     Xkcoeff2(finest_level+1);
    Vector<MultiFab>      pcoeff2(finest_level+1);
    Vector<MultiFab>   scal_force(finest_level+1);
    Vector<MultiFab>    delta_chi(finest_level+1);
    Vector<MultiFab>       sponge(finest_level+1);
    Vector<MultiFab>         w0cc(finest_level+1);

    // face-centered in the dm-direction (planar only)
    Vector<MultiFab> etarhoflux_dummy(finest_level+1);

    // face-centered
    Vector<std::array< MultiFab, AMREX_SPACEDIM > >  umac(finest_level+1);
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > sedge(finest_level+1);
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > sflux(finest_level+1);

    ////////////////////////
    // needed for spherical routines only

    // cell-centered
    Vector<MultiFab> w0_force_cart_dummy(finest_level+1);

    // face-centered
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > w0mac(finest_level+1);
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > w0mac_dummy(finest_level+1);

    // end spherical-only MultiFabs
    ////////////////////////

    // vectors store the multilevel 1D states as one very long array
    // these are cell-centered
    RealVector     grav_cell_nph    ( (max_radial_level+1)*nr_fine );
    RealVector     rho0_nph         ( (max_radial_level+1)*nr_fine );
    RealVector     p0_nph           ( (max_radial_level+1)*nr_fine );
    RealVector     p0_minus_peosbar ( (max_radial_level+1)*nr_fine );
    RealVector     peosbar          ( (max_radial_level+1)*nr_fine );
    RealVector     w0_force_dummy   ( (max_radial_level+1)*nr_fine );
    RealVector     Sbar             ( (max_radial_level+1)*nr_fine );
    RealVector     beta0_nph        ( (max_radial_level+1)*nr_fine );
    RealVector     gamma1bar_nph    ( (max_radial_level+1)*nr_fine );
    RealVector delta_gamma1_termbar ( (max_radial_level+1)*nr_fine );
    RealVector delta_chi_w0_dummy   ( (max_radial_level+1)*nr_fine );

    // vectors store the multilevel 1D states as one very long array
    // these are edge-centered
    RealVector   w0_old             ( (max_radial_level+1)*(nr_fine+1) );
    RealVector rho0_pred_edge_dummy ( (max_radial_level+1)*(nr_fine+1) );

    // make sure C++ is as efficient as possible with memory usage
    grav_cell_nph.shrink_to_fit();
    rho0_nph.shrink_to_fit();
    p0_nph.shrink_to_fit();
    p0_minus_peosbar.shrink_to_fit();
    peosbar.shrink_to_fit();
    w0_force_dummy.shrink_to_fit();
    Sbar.shrink_to_fit();
    beta0_nph.shrink_to_fit();
    gamma1bar_nph.shrink_to_fit();
    rho0_pred_edge_dummy.shrink_to_fit();
    delta_gamma1_termbar.shrink_to_fit();
    w0_old.shrink_to_fit();
    delta_chi_w0_dummy.shrink_to_fit();

    int is_predictor;

    bool split_projection = true;

    Print() << "\nTimestep " << istep << " starts with TIME = " << t_old
	    << " DT = " << dt << std::endl << std::endl;

    if (evolve_base_state) {
	Abort("evolve_base_state not supported with SDC");
    } else if (enthalpy_pred_type == 1) {
	Abort("enthalpy_pred_type == 1 not supported with SDC");
    }

    if (maestro_verbose > 0) {
    	Print() << "Cell Count:" << std::endl;
    	for (int lev=0; lev<=finest_level; ++lev) {
    	    Print() << "Level " << lev << ", " << CountCells(lev) << " cells" << std::endl;
    	}
    }

    for (int lev=0; lev<=finest_level; ++lev) {
	// cell-centered MultiFabs
	shat        [lev].define(grids[lev], dmap[lev],   Nscal, ng_s);
	rhohalf     [lev].define(grids[lev], dmap[lev],       1,    1);
	cphalf      [lev].define(grids[lev], dmap[lev],       1,    1);
	xihalf      [lev].define(grids[lev], dmap[lev], NumSpec,    1);
	macrhs      [lev].define(grids[lev], dmap[lev],       1,    0);
	macphi      [lev].define(grids[lev], dmap[lev],       1,    1);
	S_cc_nph    [lev].define(grids[lev], dmap[lev],       1,    0);
	rho_omegadot[lev].define(grids[lev], dmap[lev], NumSpec,    0);
	diff_old    [lev].define(grids[lev], dmap[lev],       1,    0);
	diff_hat    [lev].define(grids[lev], dmap[lev],       1,    0);
	diff_hterm_old[lev].define(grids[lev], dmap[lev],     1,    0);
	diff_hterm_hat[lev].define(grids[lev], dmap[lev],     1,    0);
	rho_Hnuc    [lev].define(grids[lev], dmap[lev],       1,    0);
	rho_Hext    [lev].define(grids[lev], dmap[lev],       1,    0);
	s1          [lev].define(grids[lev], dmap[lev],   Nscal, ng_s);
	s2          [lev].define(grids[lev], dmap[lev],   Nscal, ng_s);
	intra       [lev].define(grids[lev], dmap[lev],   Nscal,    0);
	sdc_source  [lev].define(grids[lev], dmap[lev],   Nscal,    0);
	aofs        [lev].define(grids[lev], dmap[lev],   Nscal,    0);
	delta_gamma1_term[lev].define(grids[lev], dmap[lev],  1,    0);
	delta_gamma1[lev].define(grids[lev], dmap[lev],       1,    0);
        peosbar_cart[lev].define(grids[lev], dmap[lev],       1,    0);
	p0_cart     [lev].define(grids[lev], dmap[lev],       1,    0);
	delta_p_term[lev].define(grids[lev], dmap[lev],       1,    0);
	Tcoeff1     [lev].define(grids[lev], dmap[lev],       1,    1);
	hcoeff1     [lev].define(grids[lev], dmap[lev],       1,    1);
	Xkcoeff1    [lev].define(grids[lev], dmap[lev], NumSpec,    1);
	pcoeff1     [lev].define(grids[lev], dmap[lev],       1,    1);
	Tcoeff2     [lev].define(grids[lev], dmap[lev],       1,    1);
	hcoeff2     [lev].define(grids[lev], dmap[lev],       1,    1);
	Xkcoeff2    [lev].define(grids[lev], dmap[lev], NumSpec,    1);
	pcoeff2     [lev].define(grids[lev], dmap[lev],       1,    1);
	if (ppm_trace_forces == 0) {
	    scal_force  [lev].define(grids[lev], dmap[lev],   Nscal,    1);
	} else {
	    // we need more ghostcells if we are tracing the forces
	    scal_force  [lev].define(grids[lev], dmap[lev],   Nscal, ng_s);
	}
	delta_chi   [lev].define(grids[lev], dmap[lev],       1,    0);
	sponge      [lev].define(grids[lev], dmap[lev],       1,    0);
	w0cc     [lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 0);

	// face-centered in the dm-direction (planar only)
	AMREX_D_TERM(etarhoflux_dummy[lev].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1); ,
		     etarhoflux_dummy[lev].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1); ,
		     etarhoflux_dummy[lev].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1); );

	// face-centered arrays of MultiFabs
	AMREX_D_TERM(umac [lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1,     1); ,
		     umac [lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1,     1); ,
		     umac [lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1,     1); );
	AMREX_D_TERM(sedge[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], Nscal, 0); ,
		     sedge[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], Nscal, 0); ,
		     sedge[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], Nscal, 0); );
	AMREX_D_TERM(sflux[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], Nscal, 0); ,
		     sflux[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], Nscal, 0); ,
		     sflux[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], Nscal, 0); );

	// initialize umac
	for (int d=0; d < AMREX_SPACEDIM; ++d)
	    umac[lev][d].setVal(0.);
    }

#if (AMREX_SPACEDIM == 3)
    for (int lev=0; lev<=finest_level; ++lev) {
	w0mac[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);
	w0mac[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
	w0mac[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1);
	w0mac_dummy[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);
	w0mac_dummy[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
	w0mac_dummy[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1);
    }
    if (spherical == 1) {
	for (int lev=0; lev<=finest_level; ++lev) {
	    w0_force_cart_dummy[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
	}
    }
#endif


    // set etarhoflux_dummy to zero
    for (int lev=0; lev<=finest_level; ++lev) {
	etarhoflux_dummy[lev].setVal(0.);
    }

#if (AMREX_SPACEDIM == 3)
    // initialize MultiFabs and Vectors to ZERO
    for (int lev=0; lev<=finest_level; ++lev) {
	for (int d=0; d<AMREX_SPACEDIM; ++d) {
	    w0mac[lev][d].setVal(0.);
	    w0mac_dummy[lev][d].setVal(0.);
	}
    }
    if (spherical == 1) {
	for (int lev=0; lev<=finest_level; ++lev) {
	    w0_force_cart_dummy[lev].setVal(0.);
	}
    }
#endif

    // initialize to zero
    std::fill(Sbar.begin(), Sbar.end(), 0.);
    std::fill(w0.begin()  , w0.end()  , 0.);

    // set dummy variables to zero
    std::fill(w0_force_dummy.begin()      , w0_force_dummy.end()      , 0.);
    std::fill(rho0_pred_edge_dummy.begin(), rho0_pred_edge_dummy.end(), 0.);

    // make the sponge for all levels
    if (do_sponge) {
	init_sponge(rho0_old.dataPtr());
	MakeSponge(sponge);
    }

    //////////////////////////////////////////////////////////////////////////////
    // STEP 1 -- Compute advection velocities
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
	Print() << "<<< STEP 1 : Compute advection velocities >>>" << std::endl;
    }


    //////////////////////////////////////////////////////////////////////////////
    // STEP 2 -- Predictor
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
	Print() << "<<< STEP 2 : Predictor >>>" << std::endl;
    }

    //////////////////////////////////////////////////////////////////////////////
    // STEP 2A -- compute advective flux divergences
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
	Print() << "<<< STEP 2A : compute advective flux div >>>" << std::endl;
    }


    //////////////////////////////////////////////////////////////////////////////
    // STEP 2B (optional) -- compute diffusive flux divergence
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
	Print() << "<<< STEP 2B : compute diffusive flux div >>>" << std::endl;
    }

    
    //////////////////////////////////////////////////////////////////////////////
    // STEP 2C -- advance thermodynamic variables
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
	Print() << "<<< STEP 2C: advance thermo variables >>>" << std::endl;
    }


    //////////////////////////////////////////////////////////////////////////////
    // STEP 3 -- Update advection velocities
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
	Print() << "<<< STEP 3 : Update advection velocities >>>" << std::endl;
    }


    //////////////////////////////////////////////////////////////////////////////
    // STEP 4 -- Corrector loop
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
	Print() << "<<< STEP 4 : Corrector loop >>>" << std::endl;
    }


    //////////////////////////////////////////////////////////////////////////////
    // STEP 4A -- compute advective 
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
	Print() << "<<< STEP 4A : compute advective flux div >>>" << std::endl;
    }


    //////////////////////////////////////////////////////////////////////////////
    // STEP 4B (optional) -- compute diffusive flux divergences
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
	Print() << "<<< STEP 4B : compute diffusive flux div >>>" << std::endl;
    }


    //////////////////////////////////////////////////////////////////////////////
    // STEP 4C -- advance thermodynamic variables
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
	Print() << "<<< STEP 4C: advance thermo variables >>>" << std::endl;
    }


    //////////////////////////////////////////////////////////////////////////////
    // STEP 5 -- Advance velocity and dynamic pressure
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
	Print() << "<<< STEP 5 : Advance velocity and dynamic pressure >>>" << std::endl;
    }

    // Define rho at half time using the new rho from Step 4
    FillPatch(0.5*(t_old+t_new), rhohalf, sold, snew, Rho, 0, 1, Rho, bcs_s);

    VelocityAdvance(rhohalf,umac,w0mac_dummy,w0_force_dummy,w0_force_cart_dummy,
		    rho0_nph,grav_cell_nph,sponge);

    if (evolve_base_state && is_initIter) {
        // throw away w0 by setting w0 = w0_old
        w0 = w0_old;
    }

    if (spherical == 1 && evolve_base_state && split_projection) {
	// subtract w0 from uold and unew for nodal projection
	for (int lev = 0; lev <= finest_level; ++lev) {
	    MultiFab::Subtract(uold[lev],w0cc[lev],0,0,AMREX_SPACEDIM,0);
	    MultiFab::Subtract(unew[lev],w0cc[lev],0,0,AMREX_SPACEDIM,0);
	}
    }
    if (evolve_base_state && !split_projection) {
        for (int i=0; i<Sbar.size(); ++i) {
            Sbar[i] = (p0_new[i] - p0_old[i])/(dt*gamma1bar_new[i]*p0_new[i]);
        }
    }

    int proj_type;

    // Project the new velocity field
    if (is_initIter) {

	proj_type = pressure_iters_comp;

	// rhcc_for_nodalproj needs to contain
	// (beta0^nph S^1 - beta0^n S^0 ) / dt

	Vector<MultiFab> rhcc_for_nodalproj_old(finest_level+1);
	for (int lev=0; lev<=finest_level; ++lev) {
	    rhcc_for_nodalproj_old[lev].define(grids[lev], dmap[lev], 1, 1);
	    MultiFab::Copy(rhcc_for_nodalproj_old[lev], rhcc_for_nodalproj[lev], 0, 0, 1, 1);
	}

	MakeRHCCforNodalProj(rhcc_for_nodalproj,S_cc_new,Sbar,beta0_nph,delta_gamma1_term);

	for (int lev=0; lev<=finest_level; ++lev) {
	    MultiFab::Subtract(rhcc_for_nodalproj[lev], rhcc_for_nodalproj_old[lev], 0, 0, 1, 1);
	    rhcc_for_nodalproj[lev].mult(1./dt,0,1,1);
	}

    }
    else {

	proj_type = regular_timestep_comp;

	MakeRHCCforNodalProj(rhcc_for_nodalproj,S_cc_new,Sbar,beta0_nph,delta_gamma1_term);

	// compute delta_p_term = peos_new - p0_new (for RHS of projection)
	if (dpdt_factor > 0.) {
	    // peos_new now holds the thermodynamic p computed from snew(rho h X)
	    PfromRhoH(snew,snew,delta_p_term);

	    // compute peosbar = Avg(peos_new)
            Average(delta_p_term,peosbar,0);

	    // no need to compute p0_minus_peosbar since make_w0 is not called after here

	    // compute peosbar_cart from peosbar
            Put1dArrayOnCart(peosbar, peosbar_cart, 0, 0, bcs_f, 0);

            // compute delta_p_term = peos_new - peosbar_cart
            for (int lev=0; lev<=finest_level; ++lev) {
                MultiFab::Subtract(delta_p_term[lev],peosbar_cart[lev],0,0,1,0);
            }

	    CorrectRHCCforNodalProj(rhcc_for_nodalproj,rho0_new,beta0_nph,gamma1bar_new,
				    p0_new,delta_p_term);
	}
    }

    // wallclock time
    const Real start_total_nodalproj = ParallelDescriptor::second();

    // call nodal projection
    NodalProj(proj_type,rhcc_for_nodalproj);

    // wallclock time
    Real end_total_nodalproj = ParallelDescriptor::second() - start_total_nodalproj;
    ParallelDescriptor::ReduceRealMax(end_total_nodalproj,ParallelDescriptor::IOProcessorNumber());

    if (spherical == 1 && evolve_base_state && split_projection) {
    	// add w0 back to unew
    	for (int lev = 0; lev <= finest_level; ++lev) {
    	    MultiFab::Add(unew[lev],w0cc[lev],0,0,AMREX_SPACEDIM,0);
    	}
    	AverageDown(unew,0,AMREX_SPACEDIM);
    	FillPatch(t_new, unew, unew, unew, 0, 0, AMREX_SPACEDIM, 0, bcs_u, 1);
    }

    for(int i=0; i<beta0_nm1.size(); ++i) {
        beta0_nm1[i] = 0.5*(beta0_old[i]+beta0_new[i]);
    }

    if (!is_initIter) {
	if (!fix_base_state) {
	    // compute tempbar by "averaging"
	    Average(snew,tempbar,Temp);
	}
    }

    Print() << "\nTimestep " << istep << " ends with TIME = " << t_new
	    << " DT = " << dt << std::endl;

    // print wallclock time
    if (maestro_verbose > 0) {
        Print() << "Time to solve mac proj   : " << end_total_macproj << '\n';
        Print() << "Time to solve nodal proj : " << end_total_nodalproj << '\n';
        Print() << "Time to solve reactions  : " << end_total_react << '\n';
    }

}

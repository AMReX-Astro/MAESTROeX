
#include <Maestro.H>

using namespace amrex;

// advance a single level for a single time step, updates flux registers
void
Maestro::AdvanceTimeStepAverage (bool is_initIter) {

    // timer for profiling
    BL_PROFILE_VAR("Maestro::AdvanceTimeStepAverage()",AdvanceTimeStepAverage);

    // cell-centered MultiFabs needed within the AdvanceTimeStep routine
    Vector<MultiFab>      rhohalf(finest_level+1);
    Vector<MultiFab>       macrhs(finest_level+1);
    Vector<MultiFab>       macphi(finest_level+1);
    Vector<MultiFab>     S_cc_nph(finest_level+1);
    Vector<MultiFab> rho_omegadot(finest_level+1);
    Vector<MultiFab>     thermal1(finest_level+1);
    Vector<MultiFab>     thermal2(finest_level+1);
    Vector<MultiFab>     rho_Hnuc(finest_level+1);
    Vector<MultiFab>     rho_Hext(finest_level+1);
    Vector<MultiFab>           s1(finest_level+1);
    Vector<MultiFab>           s2(finest_level+1);
    Vector<MultiFab>       s2star(finest_level+1);
    Vector<MultiFab> delta_gamma1_term(finest_level+1);
    Vector<MultiFab> delta_gamma1(finest_level+1);
    Vector<MultiFab> peosbar_cart(finest_level+1);
    Vector<MultiFab>      p0_cart(finest_level+1);
    Vector<MultiFab> delta_p_term(finest_level+1);
    Vector<MultiFab>       Tcoeff(finest_level+1);
    Vector<MultiFab>      hcoeff1(finest_level+1);
    Vector<MultiFab>     Xkcoeff1(finest_level+1);
    Vector<MultiFab>      pcoeff1(finest_level+1);
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
    RealVector grav_cell_nph   ( (max_radial_level+1)*nr_fine );
    RealVector   rho0_nph        ( (max_radial_level+1)*nr_fine );
    RealVector p0_nph          ( (max_radial_level+1)*nr_fine );
    RealVector p0_minus_peosbar( (max_radial_level+1)*nr_fine );
    RealVector   peosbar         ( (max_radial_level+1)*nr_fine );
    RealVector   w0_force_dummy  ( (max_radial_level+1)*nr_fine );
    RealVector   Sbar            ( (max_radial_level+1)*nr_fine );
    RealVector   beta0_nph       ( (max_radial_level+1)*nr_fine );
    RealVector   gamma1bar_nph   ( (max_radial_level+1)*nr_fine );
    RealVector   delta_gamma1_termbar ( (max_radial_level+1)*nr_fine );
    RealVector delta_chi_w0_dummy   ( (max_radial_level+1)*nr_fine );

    // vectors store the multilevel 1D states as one very long array
    // these are edge-centered
    RealVector   w0_old             ( (max_radial_level+1)*(nr_fine+1) );
    RealVector rho0_pred_edge_dummy( (max_radial_level+1)*(nr_fine+1) );

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

    if (maestro_verbose > 0) {
    	Print() << "Cell Count:" << std::endl;
    	for (int lev=0; lev<=finest_level; ++lev) {
    	    Print() << "Level " << lev << ", " << CountCells(lev) << " cells" << std::endl;
    	}
    }

    for (int lev=0; lev<=finest_level; ++lev) {
	// cell-centered MultiFabs
	rhohalf     [lev].define(grids[lev], dmap[lev],       1,    1);
	macrhs      [lev].define(grids[lev], dmap[lev],       1,    0);
	macphi      [lev].define(grids[lev], dmap[lev],       1,    1);
	S_cc_nph    [lev].define(grids[lev], dmap[lev],       1,    0);
	rho_omegadot[lev].define(grids[lev], dmap[lev], NumSpec,    0);
	thermal1    [lev].define(grids[lev], dmap[lev],       1,    0);
	thermal2    [lev].define(grids[lev], dmap[lev],       1,    0);
	rho_Hnuc    [lev].define(grids[lev], dmap[lev],       1,    0);
	rho_Hext    [lev].define(grids[lev], dmap[lev],       1,    0);
	s1          [lev].define(grids[lev], dmap[lev],   Nscal, ng_s);
	s2          [lev].define(grids[lev], dmap[lev],   Nscal, ng_s);
	s2star      [lev].define(grids[lev], dmap[lev],   Nscal, ng_s);
	delta_gamma1_term[lev].define(grids[lev], dmap[lev],  1,    0);
	delta_gamma1[lev].define(grids[lev], dmap[lev],       1,    0);
	p0_cart     [lev].define(grids[lev], dmap[lev],       1,    0);
        peosbar_cart[lev].define(grids[lev], dmap[lev],       1,    0);
	delta_p_term[lev].define(grids[lev], dmap[lev],       1,    0);
	Tcoeff      [lev].define(grids[lev], dmap[lev],       1,    1);
	hcoeff1     [lev].define(grids[lev], dmap[lev],       1,    1);
	Xkcoeff1    [lev].define(grids[lev], dmap[lev], NumSpec,    1);
	pcoeff1     [lev].define(grids[lev], dmap[lev],       1,    1);
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
    // STEP 1 -- react the full state and then base state through dt/2
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
	Print() << "<<< STEP 1 : react state >>>" << std::endl;
    }

    // wallclock time
    Real start_total_react = ParallelDescriptor::second();

    React(sold,s1,rho_Hext,rho_omegadot,rho_Hnuc,p0_old,0.5*dt);

    // wallclock time
    Real end_total_react = ParallelDescriptor::second() - start_total_react;
    ParallelDescriptor::ReduceRealMax(end_total_react,ParallelDescriptor::IOProcessorNumber());

    //////////////////////////////////////////////////////////////////////////////
    // STEP 2 -- define average expansion at time n+1/2
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
	Print() << "<<< STEP 2 : compute provisional S >>>" << std::endl;
    }

    if (t_old == 0.) {
	// this is either a pressure iteration or the first time step
	// set S_cc_nph = (1/2) (S_cc_old + S_cc_new)
	for (int lev=0; lev<=finest_level; ++lev) {
	    MultiFab::LinComb(S_cc_nph[lev],0.5,S_cc_old[lev],0,0.5,S_cc_new[lev],0,0,1,0);
	}
    }
    else {
	// set S_cc_nph = S_cc_old + (dt/2) * dSdt
	for (int lev=0; lev<=finest_level; ++lev) {
	    MultiFab::LinComb(S_cc_nph[lev],1.0,S_cc_old[lev],0,0.5*dt,dSdt[lev],0,0,1,0);
	}
    }
    // no ghost cells for S_cc_nph
    AverageDown(S_cc_nph,0,1);

    // compute p0_minus_peosbar = p0_old - peosbar (for making w0) and
    // compute delta_p_term = peos_old - p0_old (for RHS of projections)
    if (dpdt_factor > 0.0) {
	// peos_old (delta_p_term) now holds the thermodynamic p computed from sold(rho,h,X)
	PfromRhoH(sold,sold,delta_p_term);

	// compute peosbar = Avg(peos_old)
        Average(delta_p_term,peosbar,0);

        // compute p0_minus_peosbar = p0_old - peosbar
        for (int i=0; i<p0_minus_peosbar.size(); ++i) {
            p0_minus_peosbar[i] = p0_old[i] - peosbar[i];
        }

	// compute peosbar_cart from peosbar
        Put1dArrayOnCart(peosbar, peosbar_cart, 0, 0, bcs_f, 0);

        // compute delta_p_term = peos_old - peosbar_cart
        for (int lev=0; lev<=finest_level; ++lev) {
            MultiFab::Subtract(delta_p_term[lev],peosbar_cart[lev],0,0,1,0);
        }
    }
    else {
	// these should have no effect if dpdt_factor <= 0
        std::fill(p0_minus_peosbar.begin(), p0_minus_peosbar.end(), 0.);
	for (int lev=0; lev<=finest_level; ++lev) {
	    delta_p_term[lev].setVal(0.);
	}
    }

    if (evolve_base_state) {

        if (split_projection) {

            // compute Sbar = average(S_cc_nph)
            Average(S_cc_nph,Sbar,0);

            // save old-time value
            w0_old = w0;

            // compute w0, w0_force, and delta_chi_w0
            is_predictor = 1;
            make_w0(w0.dataPtr(),w0_old.dataPtr(),w0_force_dummy.dataPtr(),Sbar.dataPtr(),
                    rho0_old.dataPtr(),rho0_old.dataPtr(),p0_old.dataPtr(),p0_old.dataPtr(),
                    gamma1bar_old.dataPtr(),gamma1bar_old.dataPtr(),p0_minus_peosbar.dataPtr(),
                    etarho_ec.dataPtr(),etarho_cc.dataPtr(),delta_chi_w0_dummy.dataPtr(),
                    r_cc_loc.dataPtr(),r_edge_loc.dataPtr(),&dt,&dtold,&is_predictor);

            if (spherical == 1) {
                // put w0 on Cartesian edges
                MakeW0mac(w0mac);
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////////
    // STEP 3 -- construct the advective velocity
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
	Print() << "<<< STEP 3 : create MAC velocities >>>" << std::endl;
    }

    // compute unprojected MAC velocities
    is_predictor = 1;
    AdvancePremac(umac,w0mac_dummy,w0_force_dummy,w0_force_cart_dummy);

    for (int lev=0; lev<=finest_level; ++lev) {
	delta_chi[lev].setVal(0.);
	macphi   [lev].setVal(0.);
	delta_gamma1_term[lev].setVal(0.);
    }

    if (evolve_base_state && !split_projection) {
	for (int i=0; i<Sbar.size(); ++i) {
            Sbar[i] = 1.0/(gamma1bar_old[i]*p0_old[i]) * (p0_old[i] - p0_nm1[i])/dtold;
	}
    }

    // compute RHS for MAC projection, beta0*(S_cc-Sbar) + beta0*delta_chi
    MakeRHCCforMacProj(macrhs,rho0_old,S_cc_nph,Sbar,beta0_old,delta_gamma1_term,
		       gamma1bar_old,p0_old,delta_p_term,delta_chi,is_predictor);

    if (spherical == 1 && evolve_base_state && split_projection) {
	// subtract w0mac from umac
	for (int lev = 0; lev <= finest_level; ++lev) {
	    for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
		MultiFab::Subtract(umac[lev][dim],w0mac[lev][dim],0,0,1,1);
	    }
	}
    }

    // wallclock time
    Real start_total_macproj = ParallelDescriptor::second();

    // MAC projection
    // includes spherical option in C++ function
    MacProj(umac,macphi,macrhs,beta0_old,is_predictor);

    // wallclock time
    Real end_total_macproj = ParallelDescriptor::second() - start_total_macproj;
    ParallelDescriptor::ReduceRealMax(end_total_macproj,ParallelDescriptor::IOProcessorNumber());

    if (spherical == 1 && evolve_base_state && split_projection) {
	// add w0mac back to umac
	for (int lev = 0; lev <= finest_level; ++lev) {
	    for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
		MultiFab::Add(umac[lev][dim],w0mac[lev][dim],0,0,1,1);
	    }
	}
    }

    //////////////////////////////////////////////////////////////////////////////
    // STEP 4 -- advect the full state through dt
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
	Print() << "<<< STEP 4 : advect base >>>" << std::endl;
    }

    // no need to advect the base state density
    rho0_new = rho0_old;

    // thermal is the forcing for rhoh or temperature
    if (use_thermal_diffusion) {
	MakeThermalCoeffs(s1,Tcoeff,hcoeff1,Xkcoeff1,pcoeff1);

	MakeExplicitThermal(thermal1,s1,Tcoeff,hcoeff1,Xkcoeff1,pcoeff1,p0_old,
			    temp_diffusion_formulation);
    }
    else {
	for (int lev=0; lev<=finest_level; ++lev) {
	    thermal1[lev].setVal(0.);
	}
    }

    // copy temperature from s1 into s2 for seeding eos calls
    // temperature will be overwritten later after enthalpy advance
    for (int lev=0; lev<=finest_level; ++lev) {
	s2[lev].setVal(0.);
	MultiFab::Copy(s2[lev],s1[lev],Temp,Temp,1,ng_s);
    }

    if (maestro_verbose >= 1) {
	Print() << "            :  density_advance >>>" << std::endl;
	Print() << "            :   tracer_advance >>>" << std::endl;
    }

    // set sedge and sflux to zero
    for (int lev=0; lev<=finest_level; ++lev) {
	for (int idim=0; idim<AMREX_SPACEDIM; ++idim) {
	    sedge[lev][idim].setVal(0.);
	    sflux[lev][idim].setVal(0.);
	}
    }

    // advect rhoX, rho, and tracers
    DensityAdvance(1,s1,s2,sedge,sflux,scal_force,etarhoflux_dummy,umac,w0mac_dummy,rho0_pred_edge_dummy);

    // no need to compute etarho
    if (evolve_base_state) {
	// correct the base state density by "averaging"
	Average(s2, rho0_new, Rho);
	compute_cutoff_coords(rho0_new.dataPtr());

	if (use_etarho) {
            // compute the new etarho
            if (spherical == 0) {
                MakeEtarho(etarho_ec,etarho_cc,etarhoflux_dummy);
            } else {
                MakeEtarhoSphr(s1,s2,umac,w0mac_dummy,etarho_ec,etarho_cc);
            }
	}
    }

    // update grav_cell_new
    if (evolve_base_state) {
	make_grav_cell(grav_cell_new.dataPtr(),
		       rho0_new.dataPtr(),
		       r_cc_loc.dataPtr(),
		       r_edge_loc.dataPtr());
    }
    else {
	grav_cell_new = grav_cell_old;
    }

    // base state pressure update
    if (evolve_base_state) {

	// set new p0 through HSE
	p0_new = p0_old;

	enforce_HSE(rho0_new.dataPtr(),
		    p0_new.dataPtr(),
		    grav_cell_new.dataPtr(),
		    r_cc_loc.dataPtr(),
		    r_edge_loc.dataPtr());

	// compute p0_nph
	for (int i=0; i<p0_nph.size(); ++i) {
	    p0_nph[i] = 0.5*(p0_old[i] + p0_new[i]);
	}

	// hold dp0/dt in psi for enthalpy advance
	for (int i=0; i<p0_old.size(); ++i) {
	    psi[i] = (p0_new[i] - p0_old[i])/dt;
	}

    }
    else {
	p0_new = p0_old;
    }

    // base state enthalpy update
    if (evolve_base_state) {
	// compute rhoh0_old by "averaging"
	Average(s1, rhoh0_old, RhoH);
	Average(s2, rhoh0_new, RhoH); // -> rhoh0_new = rhoh0_old (bad?)
    }
    else {
	rhoh0_new = rhoh0_old;
    }

    if (maestro_verbose >= 1) {
	Print() << "            : enthalpy_advance >>>" << std::endl;
    }

    EnthalpyAdvance(1,s1,s2,sedge,sflux,scal_force,umac,w0mac_dummy,thermal1);

    //////////////////////////////////////////////////////////////////////////////
    // STEP 4a (Option I) -- Add thermal conduction (only enthalpy terms)
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
	Print() << "<<< STEP 4a: thermal conduct >>>" << std::endl;
    }

    if (use_thermal_diffusion) {
	ThermalConduct(s1,s2,hcoeff1,Xkcoeff1,pcoeff1,hcoeff1,Xkcoeff1,pcoeff1);
    }

    // pass temperature through for seeding the temperature update eos call
    // pi goes along for the ride
    for (int lev=0; lev<=finest_level; ++lev) {
	MultiFab::Copy(s2[lev],s1[lev],Temp,Temp,1,ng_s);
	MultiFab::Copy(s2[lev],s1[lev],  Pi,  Pi,1,ng_s);
    }

    // now update temperature
    if (use_tfromp) {
	TfromRhoP(s2,p0_new);
    }
    else {
	TfromRhoH(s2,p0_new);
    }

    if (use_thermal_diffusion) {
	// make a copy of s2star since these are needed to compute
	// coefficients in the call to thermal_conduct_full_alg
	for (int lev=0; lev<=finest_level; ++lev) {
	    MultiFab::Copy(s2star[lev],s2[lev],0,0,Nscal,ng_s);
	}
    }

    //////////////////////////////////////////////////////////////////////////////
    // STEP 5 -- react the full state and then base state through dt/2
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
	Print() << "<<< STEP 5 : react state >>>" << std::endl;
    }

    // wallclock time
    start_total_react = ParallelDescriptor::second();

    React(s2,snew,rho_Hext,rho_omegadot,rho_Hnuc,p0_new,0.5*dt);

    // wallclock time
    end_total_react += ParallelDescriptor::second() - start_total_react;
    ParallelDescriptor::ReduceRealMax(end_total_react,ParallelDescriptor::IOProcessorNumber());

    if (evolve_base_state) {
	// compute beta0 and gamma1bar
	MakeGamma1bar(snew,gamma1bar_new,p0_new);
	make_beta0(beta0_new.dataPtr(), rho0_new.dataPtr(), p0_new.dataPtr(),
                   gamma1bar_new.dataPtr(), grav_cell_new.dataPtr());
    }
    else {
	// Just pass beta0 and gamma1bar through if not evolving base state
	beta0_new = beta0_old;
	gamma1bar_new = gamma1bar_old;
    }

    for(int i=0; i<beta0_nph.size(); ++i) {
	beta0_nph[i] = 0.5*(beta0_old[i]+beta0_new[i]);
	gamma1bar_nph[i] = 0.5*(gamma1bar_old[i]+gamma1bar_new[i]);
    }

    //////////////////////////////////////////////////////////////////////////////
    // STEP 6 -- define a new average expansion rate at n+1/2
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
	Print() << "<<< STEP 6 : make new S >>>" << std::endl;
    }

    if (evolve_base_state) {
	// reset cutoff coordinates to old time value
	compute_cutoff_coords(rho0_old.dataPtr());
    }

    if (use_thermal_diffusion) {
	MakeThermalCoeffs(snew,Tcoeff,hcoeff2,Xkcoeff2,pcoeff2);

	MakeExplicitThermal(thermal2,snew,Tcoeff,hcoeff2,Xkcoeff2,pcoeff2,p0_new,
			    temp_diffusion_formulation);
    }
    else {
	for (int lev=0; lev<=finest_level; ++lev) {
	    thermal2[lev].setVal(0.);
	}
    }

    // compute S at cell-centers
    Make_S_cc(S_cc_new,delta_gamma1_term,delta_gamma1,snew,uold,rho_omegadot,rho_Hnuc,
	      rho_Hext,thermal2,p0_new,gamma1bar_new,delta_gamma1_termbar,psi);

    // set S_cc_nph = (1/2) (S_cc_old + S_cc_new)
    for (int lev=0; lev<=finest_level; ++lev) {
	MultiFab::LinComb(S_cc_nph[lev],0.5,S_cc_old[lev],0,0.5,S_cc_new[lev],0,0,1,0);
    }
    AverageDown(S_cc_nph,0,1);

    // and delta_p_term = peos_new - p0_new (for RHS of projection)
    if (dpdt_factor > 0.) {
	// peos_new now holds the thermodynamic p computed from snew(rho,h,X)
	PfromRhoH(snew,snew,delta_p_term);

	// compute peosbar = Avg(peos_new)
        Average(delta_p_term,peosbar,0);

	// compute p0_minus_peosbar = p0_new - peosbar
        for (int i=0; i<p0_minus_peosbar.size(); ++i) {
            p0_minus_peosbar[i] = p0_new[i] - peosbar[i];
        }

        // compute peosbar_cart from peosbar
        Put1dArrayOnCart(peosbar, peosbar_cart, 0, 0, bcs_f, 0);

        // compute delta_p_term = peos_new - peosbar_cart
        for (int lev=0; lev<=finest_level; ++lev) {
            MultiFab::Subtract(delta_p_term[lev],peosbar_cart[lev],0,0,1,0);
        }
    }
    else {
	// these should have no effect if dpdt_factor <= 0
        std::fill(p0_minus_peosbar.begin(), p0_minus_peosbar.end(), 0.);
	for (int lev=0; lev<=finest_level; ++lev) {
	    delta_p_term[lev].setVal(0.);
	}
    }

    if (evolve_base_state) {

        if (split_projection) {

            // compute Sbar = average(S_cc_nph)
            Average(S_cc_nph,Sbar,0);

            // compute Sbar = Sbar + delta_gamma1_termbar
            if (use_delta_gamma1_term) {
                for(int i=0; i<Sbar.size(); ++i) {
                    Sbar[i] += delta_gamma1_termbar[i];
                }
            }

            // compute w0, w0_force, and delta_chi_w0
            is_predictor = 0;
            make_w0(w0.dataPtr(),w0_old.dataPtr(),w0_force_dummy.dataPtr(),Sbar.dataPtr(),
                    rho0_old.dataPtr(),rho0_new.dataPtr(),p0_old.dataPtr(),p0_new.dataPtr(),
                    gamma1bar_old.dataPtr(),gamma1bar_new.dataPtr(),p0_minus_peosbar.dataPtr(),
                    etarho_ec.dataPtr(),etarho_cc.dataPtr(),delta_chi_w0_dummy.dataPtr(),
                    r_cc_loc.dataPtr(),r_edge_loc.dataPtr(),&dt,&dtold,&is_predictor);

            if (spherical == 1) {
                // put w0 on Cartesian edges
                MakeW0mac(w0mac);
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////////
    // STEP 7 -- redo the construction of the advective velocity
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
	Print() << "<<< STEP 7 : create MAC velocities >>>" << std::endl;
    }

    // compute unprojected MAC velocities
    is_predictor = 0;
    AdvancePremac(umac,w0mac_dummy,w0_force_dummy,w0_force_cart_dummy);

    if (evolve_base_state && !split_projection) {
	for (int i=0; i<Sbar.size(); ++i) {
	    Sbar[i] = (1.0/(gamma1bar_nph[i]*p0_nph[i]))*(p0_new[i] - p0_old[i])/dt;
	}
    }

    // compute RHS for MAC projection, beta0*(S_cc-Sbar) + beta0*delta_chi
    MakeRHCCforMacProj(macrhs,rho0_new,S_cc_nph,Sbar,beta0_nph,delta_gamma1_term,
		       gamma1bar_new,p0_new,delta_p_term,delta_chi,is_predictor);

    if (spherical == 1 && evolve_base_state && split_projection) {
	// subtract w0mac from umac
	for (int lev = 0; lev <= finest_level; ++lev) {
	    for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
		MultiFab::Subtract(umac[lev][dim],w0mac[lev][dim],0,0,1,1);
	    }
	}
    }

    // wallclock time
    start_total_macproj = ParallelDescriptor::second();

    // MAC projection
    // includes spherical option in C++ function
    MacProj(umac,macphi,macrhs,beta0_nph,is_predictor);

    // wallclock time
    end_total_macproj += ParallelDescriptor::second() - start_total_macproj;
    ParallelDescriptor::ReduceRealMax(end_total_macproj,ParallelDescriptor::IOProcessorNumber());

    if (spherical == 1 && evolve_base_state && split_projection) {
	// add w0mac back to umac
	for (int lev = 0; lev <= finest_level; ++lev) {
	    for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
		MultiFab::Add(umac[lev][dim],w0mac[lev][dim],0,0,1,1);
	    }
	}
    }

    //////////////////////////////////////////////////////////////////////////////
    // STEP 8 -- advect the full state through dt
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
	Print() << "<<< STEP 8 : advect base >>>" << std::endl;
    }

    // no need to advect the base state density
    rho0_new = rho0_old;

    // copy temperature from s1 into s2 for seeding eos calls
    // temperature will be overwritten later after enthalpy advance
    for (int lev=0; lev<=finest_level; ++lev) {
	MultiFab::Copy(s2[lev],s1[lev],Temp,Temp,1,ng_s);
    }

    if (maestro_verbose >= 1) {
	Print() << "            :  density_advance >>>" << std::endl;
	Print() << "            :   tracer_advance >>>" << std::endl;
    }

    // advect rhoX, rho, and tracers
    DensityAdvance(2,s1,s2,sedge,sflux,scal_force,etarhoflux_dummy,umac,w0mac_dummy,rho0_pred_edge_dummy);

    if (evolve_base_state) {
	// correct the base state density by "averaging"
	Average(s2, rho0_new, Rho);
	compute_cutoff_coords(rho0_new.dataPtr());

	if (use_etarho) {
	    // compute the new etarho
            if (spherical == 0) {
                MakeEtarho(etarho_ec,etarho_cc,etarhoflux_dummy);
            } else {
                MakeEtarhoSphr(s1,s2,umac,w0mac_dummy,etarho_ec,etarho_cc);
            }
	}
    }

    // update grav_cell_new, rho0_nph, grav_cell_nph
    if (evolve_base_state) {
	make_grav_cell(grav_cell_new.dataPtr(),
		       rho0_new.dataPtr(),
		       r_cc_loc.dataPtr(),
		       r_edge_loc.dataPtr());

	for(int i=0; i<beta0_nph.size(); ++i) {
	    rho0_nph[i] = 0.5*(rho0_old[i]+rho0_new[i]);
	}

	make_grav_cell(grav_cell_nph.dataPtr(),
		       rho0_nph.dataPtr(),
		       r_cc_loc.dataPtr(),
		       r_edge_loc.dataPtr());
    } else {
	rho0_nph = rho0_old;
	grav_cell_nph = grav_cell_old;
    }

    // base state pressure update
    if (evolve_base_state) {

	// set new p0 through HSE
	p0_new = p0_old;

	enforce_HSE(rho0_new.dataPtr(),
		    p0_new.dataPtr(),
		    grav_cell_new.dataPtr(),
		    r_cc_loc.dataPtr(),
		    r_edge_loc.dataPtr());

	for (int i=0; i<p0_nph.size(); ++i) {
	    p0_nph[i] = 0.5*(p0_old[i] + p0_new[i]);
	}

	// hold dp0/dt in psi for enthalpy advance
	for (int i=0; i<p0_old.size(); ++i) {
            psi[i] = (p0_new[i] - p0_old[i])/dt;
        }
    }

    // base state enthalpy averaging
    if (evolve_base_state) {
    	Average(s2, rhoh0_new, RhoH);
    }

    // base state enthalpy update
    if (maestro_verbose >= 1) {
	Print() << "            : enthalpy_advance >>>" << std::endl;
    }

    EnthalpyAdvance(2,s1,s2,sedge,sflux,scal_force,umac,w0mac_dummy,thermal1);

    //////////////////////////////////////////////////////////////////////////////
    // STEP 8a (Option I) -- Add thermal conduction (only enthalpy terms)
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
	Print() << "<<< STEP 8a: thermal conduct >>>" << std::endl;
    }

    if (use_thermal_diffusion) {
	MakeThermalCoeffs(s2star,Tcoeff,hcoeff2,Xkcoeff2,pcoeff2);

	ThermalConduct(s1,s2,hcoeff1,Xkcoeff1,pcoeff1,hcoeff2,Xkcoeff2,pcoeff2);
    }

    // pass temperature through for seeding the temperature update eos call
    // pi goes along for the ride
    for (int lev=0; lev<=finest_level; ++lev) {
	MultiFab::Copy(s2[lev],s1[lev],Temp,Temp,1,ng_s);
	MultiFab::Copy(s2[lev],s1[lev],  Pi,  Pi,1,ng_s);
    }

    // now update temperature
    if (use_tfromp) {
	TfromRhoP(s2,p0_new);
    }
    else {
	TfromRhoH(s2,p0_new);
    }

    //////////////////////////////////////////////////////////////////////////////
    // STEP 9 -- react the full state and then base state through dt/2
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
	Print() << "<<< STEP 9 : react state >>>" << std::endl;
    }

    // wallclock time
    start_total_react = ParallelDescriptor::second();

    React(s2,snew,rho_Hext,rho_omegadot,rho_Hnuc,p0_new,0.5*dt);

    // wallclock time
    end_total_react += ParallelDescriptor::second() - start_total_react;
    ParallelDescriptor::ReduceRealMax(end_total_react,ParallelDescriptor::IOProcessorNumber());

    if (evolve_base_state) {
	//compute beta0 and gamma1bar
	MakeGamma1bar(snew,gamma1bar_new,p0_new);
	make_beta0(beta0_new.dataPtr(), rho0_new.dataPtr(), p0_new.dataPtr(),
                   gamma1bar_new.dataPtr(), grav_cell_new.dataPtr());
    }

    for(int i=0; i<beta0_nph.size(); ++i) {
	beta0_nph[i] = 0.5*(beta0_old[i]+beta0_new[i]);
    }

    //////////////////////////////////////////////////////////////////////////////
    // STEP 10 -- compute S^{n+1} for the final projection
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
	Print() << "<<< STEP 10: make new S >>>" << std::endl;
    }

    if (use_thermal_diffusion) {
	MakeThermalCoeffs(snew,Tcoeff,hcoeff2,Xkcoeff2,pcoeff2);

	MakeExplicitThermal(thermal2,snew,Tcoeff,hcoeff2,Xkcoeff2,pcoeff2,p0_new,
			    temp_diffusion_formulation);
    }

    Make_S_cc(S_cc_new,delta_gamma1_term,delta_gamma1,snew,uold,rho_omegadot,rho_Hnuc,
	      rho_Hext,thermal2,p0_new,gamma1bar_new,delta_gamma1_termbar,psi);

    // define dSdt = (S_cc_new - S_cc_old) / dt
    for (int lev=0; lev<=finest_level; ++lev) {
	MultiFab::LinComb(dSdt[lev],-1./dt,S_cc_old[lev],0,1./dt,S_cc_new[lev],0,0,1,0);
    }

    if (evolve_base_state) {

        if (split_projection) {

            // compute Sbar = average(S_cc_new)
            Average(S_cc_new,Sbar,0);

            // compute Sbar = Sbar + delta_gamma1_termbar
            if (use_delta_gamma1_term) {
                for(int i=0; i<Sbar.size(); ++i) {
                    Sbar[i] += delta_gamma1_termbar[i];
                }
            }

            // compute w0, w0_force, and delta_chi_w0
            is_predictor = 0;
            make_w0(w0.dataPtr(),w0_old.dataPtr(),w0_force_dummy.dataPtr(),Sbar.dataPtr(),
                    rho0_new.dataPtr(),rho0_new.dataPtr(),p0_new.dataPtr(),p0_new.dataPtr(),
                    gamma1bar_new.dataPtr(),gamma1bar_new.dataPtr(),p0_minus_peosbar.dataPtr(),
                    etarho_ec.dataPtr(),etarho_cc.dataPtr(),delta_chi_w0_dummy.dataPtr(),
                    r_cc_loc.dataPtr(),r_edge_loc.dataPtr(),&dt,&dtold,&is_predictor);

            if (spherical == 1) {
                // put w0 on Cartesian cell-centers
                Put1dArrayOnCart(w0, w0cc, 1, 1, bcs_u, 0, 1);
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////////
    // STEP 11 -- update the velocity
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
	Print() << "<<< STEP 11: update and project new velocity >>>" << std::endl;
    }

    // Define rho at half time using the new rho from Step 8
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

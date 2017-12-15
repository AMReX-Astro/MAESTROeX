
#include <Maestro.H>

using namespace amrex;

// advance a single level for a single time step, updates flux registers
void
Maestro::AdvanceTimeStep (bool is_initIter)
{

    // features to be added later:
    // -delta_gamma1_term

    // MultiFabs needed within the AdvanceTimeStep routine
    Vector<MultiFab>           rhohalf(finest_level+1);
    Vector<MultiFab>            macrhs(finest_level+1);
    Vector<MultiFab>            macphi(finest_level+1);
    Vector<MultiFab>          S_cc_nph(finest_level+1);
    Vector<MultiFab>      rho_omegadot(finest_level+1);
    Vector<MultiFab>          thermal1(finest_level+1);
    Vector<MultiFab>          thermal2(finest_level+1);
    Vector<MultiFab>          rho_Hnuc(finest_level+1);
    Vector<MultiFab>          rho_Hext(finest_level+1);
    Vector<MultiFab>                s1(finest_level+1);
    Vector<MultiFab>                s2(finest_level+1);
    Vector<MultiFab>            s2star(finest_level+1);
    Vector<MultiFab>        beta0_cart(finest_level+1);
    Vector<MultiFab>      peosbar_cart(finest_level+1);
    Vector<MultiFab>      delta_p_term(finest_level+1);
    Vector<MultiFab>            Tcoeff(finest_level+1);
    Vector<MultiFab>           hcoeff1(finest_level+1);
    Vector<MultiFab>          Xkcoeff1(finest_level+1);
    Vector<MultiFab>           pcoeff1(finest_level+1);
    Vector<MultiFab>           hcoeff2(finest_level+1);
    Vector<MultiFab>          Xkcoeff2(finest_level+1);
    Vector<MultiFab>           pcoeff2(finest_level+1);
    Vector<MultiFab>        scal_force(finest_level+1);
    Vector<MultiFab>         delta_chi(finest_level+1);
    Vector<MultiFab>            sponge(finest_level+1);
    Vector<MultiFab> delta_gamma1_term(finest_level+1);

    // face-centered in the dm-direction (planar only)
    Vector<MultiFab> etarhoflux(finest_level+1);

    // nodal
    Vector<MultiFab>     nodalrhs(finest_level+1);
    Vector<MultiFab> nodalrhs_old(finest_level+1);
    Vector<MultiFab>           pi(finest_level+1);

    // face-centered
    Vector<std::array< MultiFab, AMREX_SPACEDIM > >            umac(finest_level+1);
    Vector<std::array< MultiFab, AMREX_SPACEDIM > >           sedge(finest_level+1);
    Vector<std::array< MultiFab, AMREX_SPACEDIM > >           sflux(finest_level+1);
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > beta0_cart_edge(finest_level+1);

    ////////////////////////
    // needed for spherical routines only

    Vector<MultiFab> w0_force_cart(finest_level+1);

    // face-centered
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > w0mac(finest_level+1);

    // end spherical-only MultiFabs
    ////////////////////////

    // vectors store the multilevel 1D states as one very long array
    // these are cell-centered
    Vector<Real> grav_cell_nph   ( (max_radial_level+1)*nr_fine );
    Vector<Real> rho0_nph        ( (max_radial_level+1)*nr_fine );
    Vector<Real> p0_nph          ( (max_radial_level+1)*nr_fine );
    Vector<Real> p0_minus_peosbar( (max_radial_level+1)*nr_fine );
    Vector<Real> peosbar         ( (max_radial_level+1)*nr_fine );
    Vector<Real> w0_force        ( (max_radial_level+1)*nr_fine );
    Vector<Real> Sbar            ( (max_radial_level+1)*nr_fine );
    Vector<Real> beta0_nph       ( (max_radial_level+1)*nr_fine );
    Vector<Real> gamma1bar_temp1 ( (max_radial_level+1)*nr_fine );
    Vector<Real> gamma1bar_temp2 ( (max_radial_level+1)*nr_fine );
    Vector<Real> delta_chi_w0    ( (max_radial_level+1)*nr_fine );

    // vectors store the multilevel 1D states as one very long array
    // these are edge-centered
    Vector<Real> w0_old             ( (max_radial_level+1)*(nr_fine+1) );
    Vector<Real> beta0_edge         ( (max_radial_level+1)*(nr_fine+1) );
    Vector<Real> rho0_predicted_edge( (max_radial_level+1)*(nr_fine+1) );

    // wallclock time
    const Real strt_total = ParallelDescriptor::second();

    Print() << "\nTimestep " << istep << " starts with TIME = " << t_old
            << " DT = " << dt << endl << endl;

    if (maestro_verbose > 0) {
        Print() << "Cell Count:" << endl;
        for (int lev=0; lev<=finest_level; ++lev) {
            Print() << "Level " << lev << ", " << CountCells(lev) << " cells" << endl;
        }
    }

    for (int lev=0; lev<=finest_level; ++lev) {
        // cell-centered MultiFabs
        rhohalf          [lev].define(grids[lev], dmap[lev],       1, 1);
        macrhs           [lev].define(grids[lev], dmap[lev],       1, 0);
        macphi           [lev].define(grids[lev], dmap[lev],       1, 1);
        S_cc_nph         [lev].define(grids[lev], dmap[lev],       1, 0);
        rho_omegadot     [lev].define(grids[lev], dmap[lev], NumSpec, 0);
        thermal1         [lev].define(grids[lev], dmap[lev],       1, 0);
        thermal2         [lev].define(grids[lev], dmap[lev],       1, 0);
        rho_Hnuc         [lev].define(grids[lev], dmap[lev],       1, 0);
        rho_Hext         [lev].define(grids[lev], dmap[lev],       1, 0);
        s1               [lev].define(grids[lev], dmap[lev],   Nscal, 0);
        s2               [lev].define(grids[lev], dmap[lev],   Nscal, 0);
        s2star           [lev].define(grids[lev], dmap[lev],   Nscal, 0);
        beta0_cart       [lev].define(grids[lev], dmap[lev],       1, 1);
        peosbar_cart     [lev].define(grids[lev], dmap[lev],       1, 0);
        delta_p_term     [lev].define(grids[lev], dmap[lev],       1, 0);
        Tcoeff           [lev].define(grids[lev], dmap[lev],       1, 1);
        hcoeff1          [lev].define(grids[lev], dmap[lev],       1, 1);
        Xkcoeff1         [lev].define(grids[lev], dmap[lev], NumSpec, 1);
        pcoeff1          [lev].define(grids[lev], dmap[lev],       1, 1);
        hcoeff2          [lev].define(grids[lev], dmap[lev],       1, 1);
        Xkcoeff2         [lev].define(grids[lev], dmap[lev],       1, 1);
        pcoeff2          [lev].define(grids[lev], dmap[lev], NumSpec, 1);
        scal_force       [lev].define(grids[lev], dmap[lev],   Nscal, 1);
        delta_chi        [lev].define(grids[lev], dmap[lev],       1, 0);
        sponge           [lev].define(grids[lev], dmap[lev],       1, 0);
        delta_gamma1_term[lev].define(grids[lev], dmap[lev],       1, 0);

        // nodal MultiFabs
        nodalrhs[lev].define    (convert(grids[lev],nodal_flag), dmap[lev], 1, 0);
        nodalrhs_old[lev].define(convert(grids[lev],nodal_flag), dmap[lev], 1, 0);
        pi[lev].define          (convert(grids[lev],nodal_flag), dmap[lev], 1, 0);

        // face-centered in the dm-direction (planar only)
#if (AMREX_SPACEDIM == 2)
        etarhoflux[lev].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
#elif (AMREX_SPACEDIM == 3)
        etarhoflux[lev].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1);
#endif

        // face-centered arrays of MultiFabs
        umac           [lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);
        umac           [lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
        sedge          [lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);
        sedge          [lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
        sflux          [lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);
        sflux          [lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
        beta0_cart_edge[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);
        beta0_cart_edge[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
#if (AMREX_SPACEDIM == 3)
        umac           [lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1);
        sedge          [lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1);
        sflux          [lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1);
        beta0_cart_edge[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1);
#endif
    }

#if (AMREX_SPACEDIM == 3)
    if (spherical == 1) {
        for (int lev=0; lev<=finest_level; ++lev) {
            w0_force_cart[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
            w0mac[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);
            w0mac[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
            w0mac[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1);
        }
    }
#endif

    // initialize MultiFabs and Vectors to ZERO
    for (int lev=0; lev<=finest_level; ++lev) {
        delta_p_term     [lev].setVal(0.);
        delta_gamma1_term[lev].setVal(0.);
        delta_chi        [lev].setVal(0.);
        macphi           [lev].setVal(0.);
        thermal1         [lev].setVal(0.);
        thermal2         [lev].setVal(0.);
        etarhoflux       [lev].setVal(0.);
    }
    std::fill(p0_minus_peosbar.begin(), p0_minus_peosbar.end(), 0.);
    std::fill(w0_force        .begin(), w0_force        .end(), 0.);
    std::fill(Sbar            .begin(), Sbar            .end(), 0.);

#if (AMREX_SPACEDIM == 3)
    if (spherical == 1) {
        for (int lev=0; lev<=finest_level; ++lev) {
            w0_force_cart[lev].setVal(0.);
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                w0mac[lev][d].setVal(0.);
            }
        }
    }
#endif

    // make the sponge for all levels
    if (do_sponge) {
        init_sponge(rho0_old.dataPtr());
        MakeSponge(sponge);
    }

    //////////////////////////////////////////////////////////////////////////////
    // STEP 1 -- react the full state and then base state through dt/2
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // STEP 2 -- define average expansion at time n+1/2
    //////////////////////////////////////////////////////////////////////////////

    if (t_old == 0.0) {
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

    if (dpdt_factor > 0.0) {
        Abort("MaestroAdvance.cpp: dpdt_factor not implemented");
    }

    //////////////////////////////////////////////////////////////////////////////
    // STEP 3 -- construct the advective velocity
    //////////////////////////////////////////////////////////////////////////////

    AdvancePremac(umac);

    //////////////////////////////////////////////////////////////////////////////
    // STEP 4 -- advect the base state and full state through dt
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // STEP 4a (Option I) -- Add thermal conduction (only enthalpy terms)
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // STEP 5 -- react the full state and then base state through dt/2
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // STEP 6 -- define a new average expansion rate at n+1/2
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // STEP 7 -- redo the construction of the advective velocity using the current w0
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // STEP 8 -- advect the base state and full state through dt
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // STEP 8a (Option I) -- Add thermal conduction (only enthalpy terms)
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // STEP 9 -- react the full state and then base state through dt/2
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // STEP 10 -- compute S^{n+1} for the final projection
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // STEP 11 -- update the velocity
    //////////////////////////////////////////////////////////////////////////////


        


    for (int lev=0; lev<=finest_level; ++lev) 
    {

        MultiFab& S_new = snew[lev];

        MultiFab fluxes[AMREX_SPACEDIM];
        if (do_reflux)
        {
            for (int i = 0; i < AMREX_SPACEDIM; ++i)
            {
                BoxArray ba = grids[lev];
                ba.surroundingNodes(i);
                fluxes[i].define(ba, dmap[lev], S_new.nComp(), 0);
            }
        }

        // advection

        // increment or decrement the flux registers by area and time-weighted fluxes
        // Note that the fluxes have already been scaled by dt and area
        // In this example we are solving s_t = -div(+F)
        // The fluxes contain, e.g., F_{i+1/2,j} = (s*u)_{i+1/2,j}
        // Keep this in mind when considering the different sign convention for updating
        // the flux registers from the coarse or fine grid perspective
        // NOTE: the flux register associated with flux_reg_s[lev] is associated
        // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
        if (do_reflux) { 
            if (lev < finest_level)
            {
                for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                    // update the lev+1/lev flux register (index lev+1)   
                    flux_reg_s[lev+1]->CrseInit(fluxes[i],i,0,0,fluxes[i].nComp(), -1.0);
                }	
            }
            if (lev > 0)
            {
                for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                    // update the lev/lev-1 flux register (index lev) 
                    flux_reg_s[lev]->FineAdd(fluxes[i],i,0,0,fluxes[i].nComp(), 1.0);
                }
            }
        }


    } // end loop over levels

    // synchronize by refluxing and averaging down, starting from the finest_level-1/finest_level pair
    for (int lev=finest_level-1; lev>=0; --lev)
    {
        if (do_reflux) {
            // update lev based on coarse-fine flux mismatch
            flux_reg_s[lev+1]->Reflux(snew[lev], 1.0, 0, 0, snew[lev].nComp(), geom[lev]);
        }

        AverageDownTo(lev,snew,0,Nscal); // average lev+1 down to lev
    }

    Print() << "\nTimestep " << istep << " ends with TIME = " << t_new
            << " DT = " << dt << endl;

    // wallclock time
    Real end_total = ParallelDescriptor::second() - strt_total;
	
    // print wallclock time
    ParallelDescriptor::ReduceRealMax(end_total ,ParallelDescriptor::IOProcessorNumber());
    if (maestro_verbose > 0) {
        Print() << "Time to advance time step: " << end_total << '\n';
    }

}

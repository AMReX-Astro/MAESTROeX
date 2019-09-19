
#include <Maestro.H>

#ifdef AMREX_USE_CUDA
#include <cuda_runtime_api.h>
#include <AMReX_Arena.H>
#endif

using namespace amrex;

void
Maestro::EstDt ()
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::EstDt()",EstDt);

#ifdef AMREX_USE_CUDA
    auto not_launched = Gpu::notInLaunchRegion();
    // turn on GPU
    if (not_launched) Gpu::setLaunchRegion(true);
#endif

    dt = 1.e20;

    // allocate a dummy w0_force and set equal to zero
    RealVector w0_force_dummy( (max_radial_level+1)*nr_fine );
    w0_force_dummy.shrink_to_fit();
    std::fill(w0_force_dummy.begin(),w0_force_dummy.end(), 0.);

    // build dummy w0_force_cart and set equal to zero
    Vector<MultiFab> w0_force_cart_dummy(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        w0_force_cart_dummy[lev].define(grids[lev], dmap[lev], 3, 1);
        w0_force_cart_dummy[lev].setVal(0.);
    }

    // build a dummy umac and set equal to zero
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > umac_dummy(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        umac_dummy[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);
        umac_dummy[lev][0].setVal(0.);
        umac_dummy[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
        umac_dummy[lev][1].setVal(0.);
#if (AMREX_SPACEDIM == 3)
        umac_dummy[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1);
        umac_dummy[lev][2].setVal(0.);
#endif
    }

    // face-centered
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > w0mac(finest_level+1);

    if (spherical == 1) {
        // initialize
        for (int lev=0; lev<=finest_level; ++lev) {
            w0mac[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);
            w0mac[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
#if (AMREX_SPACEDIM == 3)
            w0mac[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1);
#endif
        }

        for (int lev=0; lev<=finest_level; ++lev) {
            for (int idim=0; idim<AMREX_SPACEDIM; ++idim) {
                w0mac[lev][idim].setVal(0.);
            }
        }

        if (evolve_base_state && (use_exact_base_state == 0 && average_base_state == 0)) {
            MakeW0mac(w0mac);
        }
    }

    // build and compute vel_force
    Vector<MultiFab> vel_force(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        vel_force[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        // needed to avoid NaNs in filling corner ghost cells with 2 physical boundaries
        vel_force[lev].setVal(0.);
    }

    int do_add_utilde_force = 0;
    MakeVelForce(vel_force,umac_dummy,sold,rho0_old,grav_cell_old,
                 w0_force_cart_dummy,do_add_utilde_force);

#if (AMREX_SPACEDIM == 3)
    // build and initialize grad_p0 for spherical case
    Vector<MultiFab> gp0_cart(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        gp0_cart[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
	gp0_cart[lev].setVal(0.);
    }
    RealVector gp0( (max_radial_level+1)*(nr_fine+1) );
    gp0.shrink_to_fit();
    std::fill(gp0.begin(),gp0.end(), 0.);

    // divU constraint
    estdt_divu(gp0.dataPtr(), p0_old.dataPtr(), gamma1bar_old.dataPtr(),
	       r_cc_loc.dataPtr(), r_edge_loc.dataPtr());

    Put1dArrayOnCart (gp0,gp0_cart,1,1,bcs_f,0);
#endif

    Vector<MultiFab> p0_cart(finest_level+1);
    Vector<MultiFab> gamma1bar_cart(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        p0_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        gamma1bar_cart[lev].define(grids[lev], dmap[lev], 1, 1);
    }

    Put1dArrayOnCart(p0_old,p0_cart,0,0,bcs_f,0);
    Put1dArrayOnCart(gamma1bar_old,gamma1bar_cart,0,0,bcs_f,0);

    Real umax = 0.;

    Real dt_lev = 1.e99;
    Real umax_lev = 0.;

    for (int lev = 0; lev <= finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& uold_mf = uold[lev];
        MultiFab& sold_mf = sold[lev];
        MultiFab& vel_force_mf = vel_force[lev];
        MultiFab& S_cc_old_mf = S_cc_old[lev];
        MultiFab& dSdt_mf = dSdt[lev];
#if (AMREX_SPACEDIM == 3)
        MultiFab& w0macx_mf = w0mac[lev][0];
        MultiFab& w0macy_mf = w0mac[lev][1];
        MultiFab& w0macz_mf = w0mac[lev][2];
        const MultiFab& gp0_cart_mf = gp0_cart[lev];
#endif
        const MultiFab& w0_mf = w0_cart[lev];
        const MultiFab& p0_mf = p0_cart[lev];
        const MultiFab& gamma1bar_mf = gamma1bar_cart[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel reduction(min:dt_lev) reduction(max:umax_lev)
#endif
	{

        Real dt_grid = 1.e99;
        Real umax_grid = 0.;

        for ( MFIter mfi(uold_mf, true); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const Real* dx = geom[lev].CellSize();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            if (spherical == 0) {
#pragma gpu box(tileBox)
                estdt(lev, AMREX_MFITER_REDUCE_MIN(&dt_grid),
                      AMREX_MFITER_REDUCE_MAX(&umax_grid),
                      AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),
                      AMREX_REAL_ANYD(dx),
                      BL_TO_FORTRAN_ANYD(sold_mf[mfi]), sold_mf[mfi].nCompPtr(),
                      BL_TO_FORTRAN_ANYD(uold_mf[mfi]), uold_mf[mfi].nCompPtr(),
                      BL_TO_FORTRAN_ANYD(vel_force_mf[mfi]), vel_force_mf[mfi].nCompPtr(),
                      BL_TO_FORTRAN_ANYD(S_cc_old_mf[mfi]),
                      BL_TO_FORTRAN_ANYD(dSdt_mf[mfi]),
                      BL_TO_FORTRAN_ANYD(w0_mf[mfi]),
                      BL_TO_FORTRAN_ANYD(p0_mf[mfi]),
                      BL_TO_FORTRAN_ANYD(gamma1bar_mf[mfi]));
            } else {
#if (AMREX_SPACEDIM == 3)

#pragma gpu box(tileBox)
		estdt_sphr(AMREX_MFITER_REDUCE_MIN(&dt_grid),
			   AMREX_MFITER_REDUCE_MAX(&umax_grid),
			   AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),
			   AMREX_REAL_ANYD(dx),
			   BL_TO_FORTRAN_ANYD(sold_mf[mfi]), sold_mf[mfi].nCompPtr(),
			   BL_TO_FORTRAN_ANYD(uold_mf[mfi]), uold_mf[mfi].nCompPtr(),
			   BL_TO_FORTRAN_ANYD(vel_force_mf[mfi]), vel_force_mf[mfi].nCompPtr(),
			   BL_TO_FORTRAN_ANYD(S_cc_old_mf[mfi]),
			   BL_TO_FORTRAN_ANYD(dSdt_mf[mfi]),
               BL_TO_FORTRAN_ANYD(w0_mf[mfi]),
			   BL_TO_FORTRAN_ANYD(w0macx_mf[mfi]),
			   BL_TO_FORTRAN_ANYD(w0macy_mf[mfi]),
			   BL_TO_FORTRAN_ANYD(w0macz_mf[mfi]),
			   BL_TO_FORTRAN_ANYD(gp0_cart_mf[mfi]));

#else
                Abort("EstDt: Spherical is not valid for DIM < 3");
#endif
            }
        }

	dt_lev = std::min(dt_lev,dt_grid);
	umax_lev = std::max(umax_lev,umax_grid);

	} //end openmp

        // find the smallest dt over all processors
        ParallelDescriptor::ReduceRealMin(dt_lev);

        // find the largest umax over all processors
        ParallelDescriptor::ReduceRealMax(umax_lev);

        // update umax over all levels
        umax = std::max(umax,umax_lev);

        if (maestro_verbose > 0) {
            Print() << "Call to estdt for level " << lev << " gives dt_lev = " << dt_lev << std::endl;
        }

        // update dt over all levels
        dt = std::min(dt,dt_lev);

    }     // end loop over levels

#ifdef AMREX_USE_CUDA
    // turn off GPU
    if (not_launched) Gpu::setLaunchRegion(false);
#endif

    if (maestro_verbose > 0) {
        Print() << "Minimum estdt over all levels = " << dt << std::endl;
    }

    if (dt < small_dt) {
        Abort("EstDt: dt < small_dt");
    }

    if (dt > max_dt) {
        if (maestro_verbose > 0) {
            Print() << "max_dt limits the new dt = " << max_dt << std::endl;
        }
        dt = max_dt;
    }

    if (fixed_dt != -1.0) {
        // fixed dt
        dt = fixed_dt;
        if (maestro_verbose > 0) {
            Print() << "Setting fixed dt = " << dt << std::endl;
        }
    }

    // set rel_eps in fortran module
    umax *= 1.e-8;
    set_rel_eps(&umax);
}


void
Maestro::FirstDt ()
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::FirstDt()",FirstDt);

#ifdef AMREX_USE_CUDA
    auto not_launched = Gpu::notInLaunchRegion();
    // turn on GPU
    if (not_launched) Gpu::setLaunchRegion(true);
#endif

    dt = 1.e20;

    // allocate a dummy w0_force and set equal to zero
    RealVector w0_force_dummy( (max_radial_level+1)*nr_fine );
    w0_force_dummy.shrink_to_fit();
    std::fill(w0_force_dummy.begin(),w0_force_dummy.end(), 0.);

    // build dummy w0_force_cart and set equal to zero
    Vector<MultiFab> w0_force_cart_dummy(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        w0_force_cart_dummy[lev].define(grids[lev], dmap[lev], 3, 1);
        w0_force_cart_dummy[lev].setVal(0.);
    }

    // build a dummy umac and set equal to zero
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > umac_dummy(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        umac_dummy[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);
        umac_dummy[lev][0].setVal(0.);
        umac_dummy[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
        umac_dummy[lev][1].setVal(0.);
#if (AMREX_SPACEDIM == 3)
        umac_dummy[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1);
        umac_dummy[lev][2].setVal(0.);
#endif
    }

    // build and compute vel_force
    Vector<MultiFab> vel_force(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        vel_force[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        // needed to avoid NaNs in filling corner ghost cells with 2 physical boundaries
        vel_force[lev].setVal(0.);
    }

    int do_add_utilde_force = 0;
    MakeVelForce(vel_force,umac_dummy,sold,rho0_old,grav_cell_old,
                 w0_force_cart_dummy,do_add_utilde_force);

#if (AMREX_SPACEDIM == 3)
    // build and initialize grad_p0 for spherical case
    Vector<MultiFab> gp0_cart(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        gp0_cart[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        gp0_cart[lev].setVal(0.);
    }
    RealVector gp0( (max_radial_level+1)*(nr_fine+1) );
    gp0.shrink_to_fit();
    std::fill(gp0.begin(),gp0.end(), 0.);

    // divU constraint
    if (use_divu_firstdt) {
	estdt_divu(gp0.dataPtr(), p0_old.dataPtr(), gamma1bar_old.dataPtr(),
		   r_cc_loc.dataPtr(), r_edge_loc.dataPtr());
    }

    Put1dArrayOnCart (gp0,gp0_cart,1,1,bcs_f,0);
#endif

    Vector<MultiFab> p0_cart(finest_level+1);
    Vector<MultiFab> gamma1bar_cart(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        p0_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        gamma1bar_cart[lev].define(grids[lev], dmap[lev], 1, 1);
    }

    Put1dArrayOnCart(p0_old,p0_cart,0,0,bcs_f,0);
    Put1dArrayOnCart(gamma1bar_old,gamma1bar_cart,0,0,bcs_f,0);

    Real umax = 0.;

    Real dt_lev = 1.e99;
    Real umax_lev = 0.;

    for (int lev = 0; lev <= finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& uold_mf = uold[lev];
        const MultiFab& sold_mf = sold[lev];
        const MultiFab& vel_force_mf = vel_force[lev];
        const MultiFab& S_cc_old_mf = S_cc_old[lev];
#if (AMREX_SPACEDIM == 3)
        const MultiFab& gp0_cart_mf = gp0_cart[lev];
#endif
        const MultiFab& p0_mf = p0_cart[lev];
        const MultiFab& gamma1bar_mf = gamma1bar_cart[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel reduction(min:dt_lev) reduction(max:umax_lev)
#endif
	{

        Real dt_grid = 1.e99;
        Real umax_grid = 0.;

        for ( MFIter mfi(sold_mf,true); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const Real* dx = geom[lev].CellSize();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            if (spherical == 0 ) {

#pragma gpu box(tileBox)
                firstdt(lev,
                    AMREX_MFITER_REDUCE_MIN(&dt_grid),
                    AMREX_MFITER_REDUCE_MAX(&umax_grid),
                    AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),
                    AMREX_REAL_ANYD(dx),
                    BL_TO_FORTRAN_ANYD(sold_mf[mfi]), sold_mf[mfi].nCompPtr(),
                    BL_TO_FORTRAN_ANYD(uold_mf[mfi]), uold_mf[mfi].nCompPtr(),
                    BL_TO_FORTRAN_ANYD(vel_force_mf[mfi]), vel_force_mf[mfi].nCompPtr(),
                    BL_TO_FORTRAN_ANYD(S_cc_old_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(p0_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(gamma1bar_mf[mfi]));

            } else {
#if (AMREX_SPACEDIM == 3)

#pragma gpu box(tileBox)
                firstdt_sphr(AMREX_MFITER_REDUCE_MIN(&dt_grid),
			     AMREX_MFITER_REDUCE_MAX(&umax_grid),
                             AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),
                             AMREX_REAL_ANYD(dx),
                             BL_TO_FORTRAN_ANYD(sold_mf[mfi]), sold_mf[mfi].nCompPtr(),
                             BL_TO_FORTRAN_ANYD(uold_mf[mfi]), uold_mf[mfi].nCompPtr(),
                             BL_TO_FORTRAN_ANYD(vel_force_mf[mfi]), vel_force_mf[mfi].nCompPtr(),
                             BL_TO_FORTRAN_ANYD(S_cc_old_mf[mfi]),
                             BL_TO_FORTRAN_ANYD(gp0_cart_mf[mfi]));
#else
                Abort("FirstDt: Spherical is not valid for DIM < 3");
#endif
            }
        }

	dt_lev = std::min(dt_lev,dt_grid);
	umax_lev = std::max(umax_lev,umax_grid);

	} //end openmp

        // find the smallest dt over all processors
        ParallelDescriptor::ReduceRealMin(dt_lev);

        // find the largest umax over all processors
        ParallelDescriptor::ReduceRealMax(umax_lev);

        // update umax over all levels
        umax = std::max(umax,umax_lev);

        if (maestro_verbose > 0) {
            Print() << "Call to firstdt for level " << lev << " gives dt_lev = " << dt_lev << std::endl;
        }

        // multiply by init_shrink
        dt_lev *= init_shrink;

        if (maestro_verbose > 0) {
            Print() << "Multiplying dt_lev by init_shrink; dt_lev = " << dt_lev << std::endl;
        }

        // update dt over all levels
        dt = std::min(dt,dt_lev);

    }     // end loop over levels

#ifdef AMREX_USE_CUDA
    // turn off GPU
    if (not_launched) Gpu::setLaunchRegion(false);
#endif

    if (maestro_verbose > 0) {
        Print() << "Minimum firstdt over all levels = " << dt << std::endl;
    }

    if (dt < small_dt) {
        Abort("FirstDt: dt < small_dt");
    }

    if (dt > max_dt) {
        if (maestro_verbose > 0) {
            Print() << "max_dt limits the new dt = " << max_dt << std::endl;
        }
        dt = max_dt;
    }

    if (fixed_dt != -1.0) {
        // fixed dt
        dt = fixed_dt;
        if (maestro_verbose > 0) {
            Print() << "Setting fixed dt = " << dt << std::endl;
        }
    }

    // set rel_eps in fortran module
    umax *= 1.e-8;
    set_rel_eps(&umax);
}

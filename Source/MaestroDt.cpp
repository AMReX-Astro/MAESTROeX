
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
#if (AMREX_SPACEDIM >= 2)
        umac_dummy[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
        umac_dummy[lev][1].setVal(0.);
#endif
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
#if (AMREX_SPACEDIM >= 2)
            w0mac[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
#endif
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
    }

    int do_add_utilde_force = 0;
    MakeVelForce(vel_force,umac_dummy,sold,rho0_old,grav_cell_old,
                 w0_force_dummy,w0_force_cart_dummy,do_add_utilde_force);

    Real umax = 0.;

    for (int lev = 0; lev <= finest_level; ++lev) {
        Real dt_lev = 1.e99;
        Real umax_lev = 0.;

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
        const MultiFab& cc_to_r = cell_cc_to_r[lev];
#endif

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel reduction(min:dt_lev) reduction(max:umax_lev)
#endif
        for ( MFIter mfi(uold_mf, true); mfi.isValid(); ++mfi ) {

            Real dt_grid = 1.e99;
            Real umax_grid = 0.;

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const Real* dx = geom[lev].CellSize();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            if (spherical == 0) {
                estdt(&lev,&dt_grid,&umax_grid,
                      ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
                      ZFILL(dx),
                      BL_TO_FORTRAN_FAB(sold_mf[mfi]),
                      BL_TO_FORTRAN_FAB(uold_mf[mfi]),
                      BL_TO_FORTRAN_FAB(vel_force_mf[mfi]),
                      BL_TO_FORTRAN_3D(S_cc_old_mf[mfi]),
                      BL_TO_FORTRAN_3D(dSdt_mf[mfi]),
                      w0.dataPtr(),
                      p0_old.dataPtr(),
                      gamma1bar_old.dataPtr());
            } else {
#if (AMREX_SPACEDIM == 3)
                    estdt_sphr(&dt_grid,&umax_grid,
                               ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
                               ZFILL(dx),
                               BL_TO_FORTRAN_FAB(sold_mf[mfi]),
                               BL_TO_FORTRAN_FAB(uold_mf[mfi]),
                               BL_TO_FORTRAN_FAB(vel_force_mf[mfi]),
                               BL_TO_FORTRAN_3D(S_cc_old_mf[mfi]),
                               BL_TO_FORTRAN_3D(dSdt_mf[mfi]),
                               w0.dataPtr(),
                               BL_TO_FORTRAN_3D(w0macx_mf[mfi]),
                               BL_TO_FORTRAN_3D(w0macy_mf[mfi]),
                               BL_TO_FORTRAN_3D(w0macz_mf[mfi]),
                               p0_old.dataPtr(),
                               gamma1bar_old.dataPtr(),
                               r_cc_loc.dataPtr(),
                               r_edge_loc.dataPtr(),
                               BL_TO_FORTRAN_3D(cc_to_r[mfi]));
#else
                Abort("EstDt: Spherical is not valid for DIM < 3");
#endif
            }

            dt_lev = std::min(dt_lev,dt_grid);
            umax_lev = std::max(umax_lev,umax_grid);
        }

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

// #ifdef AMREX_USE_CUDA
//     // turn off GPU
//     Gpu::setLaunchRegion(false);
// #endif

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
#if (AMREX_SPACEDIM >= 2)
        umac_dummy[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
        umac_dummy[lev][1].setVal(0.);
#endif
#if (AMREX_SPACEDIM == 3)
        umac_dummy[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1);
        umac_dummy[lev][2].setVal(0.);
#endif
    }

    // build and compute vel_force
    Vector<MultiFab> vel_force(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        vel_force[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
    }

    int do_add_utilde_force = 0;
    MakeVelForce(vel_force,umac_dummy,sold,rho0_old,grav_cell_old,
                 w0_force_dummy,w0_force_cart_dummy,do_add_utilde_force);

    Real umax = 0.;

    for (int lev = 0; lev <= finest_level; ++lev) {
        Real dt_lev = 1.e99;
        Real umax_lev = 0.;

        // get references to the MultiFabs at level lev
        MultiFab& uold_mf = uold[lev];
        MultiFab& sold_mf = sold[lev];
        MultiFab& vel_force_mf = vel_force[lev];
        MultiFab& S_cc_old_mf = S_cc_old[lev];
        const MultiFab& cc_to_r = cell_cc_to_r[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel reduction(min:dt_lev) reduction(max:umax_lev)
#endif
        for ( MFIter mfi(sold_mf,true); mfi.isValid(); ++mfi ) {

            Real dt_grid = 1.e99;
            Real umax_grid = 0.;

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const Real* dx = geom[lev].CellSize();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            if (spherical == 0 ) {
                firstdt(&lev,&dt_grid,&umax_grid,
                        ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
                        ZFILL(dx),
                        BL_TO_FORTRAN_FAB(sold_mf[mfi]),
                        BL_TO_FORTRAN_FAB(uold_mf[mfi]),
                        BL_TO_FORTRAN_FAB(vel_force_mf[mfi]),
                        BL_TO_FORTRAN_3D(S_cc_old_mf[mfi]),
                        p0_old.dataPtr(),
                        gamma1bar_old.dataPtr());
            } else {
#if (AMREX_SPACEDIM == 3)
                firstdt_sphr(&dt_grid,&umax_grid,
                             ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
                             ZFILL(dx),
                             BL_TO_FORTRAN_FAB(sold_mf[mfi]),
                             BL_TO_FORTRAN_FAB(uold_mf[mfi]),
                             BL_TO_FORTRAN_FAB(vel_force_mf[mfi]),
                             BL_TO_FORTRAN_3D(S_cc_old_mf[mfi]),
                             p0_old.dataPtr(),
                             gamma1bar_old.dataPtr(),
                             r_cc_loc.dataPtr(), r_edge_loc.dataPtr(),
                             BL_TO_FORTRAN_3D(cc_to_r[mfi]));
#else
                Abort("FirstDt: Spherical is not valid for DIM < 3");
#endif
            }

            dt_lev = std::min(dt_lev,dt_grid);
            umax_lev = std::max(umax_lev,umax_grid);
        }

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

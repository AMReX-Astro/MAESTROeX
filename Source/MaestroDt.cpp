
#include <Maestro.H>

using namespace amrex;

void
Maestro::EstDt ()
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::EstDt()",EstDt);

    dt = 1.e20;

    // allocate a dummy w0_force and set equal to zero
    Vector<Real> w0_force_dummy( (max_radial_level+1)*nr_fine );
    w0_force_dummy.shrink_to_fit();
    std::fill(w0_force_dummy.begin(),w0_force_dummy.end(), 0.);

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
    }

    int do_add_utilde_force = 0;
    MakeVelForce(vel_force,umac_dummy,sold,rho0_old,grav_cell_old,w0_force_dummy,do_add_utilde_force);
    
    Real dt_lev = 1.e99;
    Real umax_lev = 0.;
    Real umax = 0.;

    for (int lev = 0; lev <= finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& uold_mf = uold[lev];
        MultiFab& sold_mf = sold[lev];
        MultiFab& vel_force_mf = vel_force[lev];
        MultiFab& S_cc_old_mf = S_cc_old[lev];
        MultiFab& dSdt_mf = dSdt[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(uold_mf); mfi.isValid(); ++mfi ) {

            Real dt_grid = 1.e99;
            Real umax_grid = 0.;

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();

            const Real* dx = geom[lev].CellSize();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            estdt(&lev,&dt_grid,&umax_grid,
                  ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                  ZFILL(dx),
                  BL_TO_FORTRAN_FAB(sold_mf[mfi]),
                  BL_TO_FORTRAN_FAB(uold_mf[mfi]),
                  BL_TO_FORTRAN_FAB(vel_force_mf[mfi]),
                  BL_TO_FORTRAN_3D(S_cc_old_mf[mfi]),
                  BL_TO_FORTRAN_3D(dSdt_mf[mfi]),
                  w0.dataPtr(),
                  p0_old.dataPtr(),
                  gamma1bar_old.dataPtr());

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
            Print() << "Call to estdt for level " << lev << " gives dt_lev = " << dt_lev << endl;
        }

        // multiply by init_shrink
        dt_lev *= init_shrink;

        if (maestro_verbose > 0) {
            Print() << "Multiplying dt_lev by init_shrink; dt_lev = " << dt_lev << endl;
        }

        // update dt over all levels
        dt = std::min(dt,dt_lev);

    }  // end loop over levels
   
    if (maestro_verbose > 0) {
        Print() << "Minimum estdt over all levels = " << dt << endl;
    }

    if (dt < small_dt) {
        Abort("EstDt: dt < small_dt");
    }

    if (dt > max_dt) {
        if (maestro_verbose > 0) {
            Print() << "max_dt limits the new dt = " << max_dt << endl;
        }
        dt = max_dt;
    }

    if (fixed_dt != -1.0) {
        // fixed dt
        dt = fixed_dt;
        if (maestro_verbose > 0) {
            Print() << "Setting fixed dt = " << dt << endl;
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
    Vector<Real> w0_force_dummy( (max_radial_level+1)*nr_fine );
    w0_force_dummy.shrink_to_fit();
    std::fill(w0_force_dummy.begin(),w0_force_dummy.end(), 0.);

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
    }

    int do_add_utilde_force = 0;
    MakeVelForce(vel_force,umac_dummy,sold,rho0_old,grav_cell_old,w0_force_dummy,do_add_utilde_force);
    
    Real dt_lev = 1.e99;
    Real umax_lev = 0.;
    Real umax = 0.;

    for (int lev = 0; lev <= finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& uold_mf = uold[lev];
        MultiFab& sold_mf = sold[lev];
        MultiFab& vel_force_mf = vel_force[lev];
        MultiFab& S_cc_old_mf = S_cc_old[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(uold_mf); mfi.isValid(); ++mfi ) {

            Real dt_grid = 1.e99;
            Real umax_grid = 0.;

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();

            const Real* dx = geom[lev].CellSize();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            firstdt(&lev,&dt_grid,&umax_grid,
                    ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                    ZFILL(dx),
                    BL_TO_FORTRAN_FAB(sold_mf[mfi]),
                    BL_TO_FORTRAN_FAB(uold_mf[mfi]),
                    BL_TO_FORTRAN_FAB(vel_force_mf[mfi]),
                    BL_TO_FORTRAN_3D(S_cc_old_mf[mfi]),
                    p0_old.dataPtr(),
                    gamma1bar_old.dataPtr());

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
            Print() << "Call to firstdt for level " << lev << " gives dt_lev = " << dt_lev << endl;
        }

        // multiply by init_shrink
        dt_lev *= init_shrink;

        if (maestro_verbose > 0) {
            Print() << "Multiplying dt_lev by init_shrink; dt_lev = " << dt_lev << endl;
        }

        // update dt over all levels
        dt = std::min(dt,dt_lev);

    }  // end loop over levels
   
    if (maestro_verbose > 0) {
        Print() << "Minimum firstdt over all levels = " << dt << endl;
    }

    if (dt < small_dt) {
        Abort("FirstDt: dt < small_dt");
    }

    if (dt > max_dt) {
        if (maestro_verbose > 0) {
            Print() << "max_dt limits the new dt = " << max_dt << endl;
        }
        dt = max_dt;
    }

    if (fixed_dt != -1.0) {
        // fixed dt
        dt = fixed_dt;
        if (maestro_verbose > 0) {
            Print() << "Setting fixed dt = " << dt << endl;
        }
    }

    // set rel_eps in fortran module
    umax *= 1.e-8;
    set_rel_eps(&umax);
}

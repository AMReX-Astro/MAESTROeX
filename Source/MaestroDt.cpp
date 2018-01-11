
#include <Maestro.H>

using namespace amrex;

// a wrapper for EstDtLevel
void
Maestro::EstDt ()
{
    dt = 1.e99;
        
    Real umax = 0.;

    for (int lev = 0; lev <= finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& uold_mf = uold[lev];
        MultiFab& sold_mf = sold[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(uold_mf); mfi.isValid(); ++mfi ) {

            Real dt_grid = 1.e99;
            Real umax_grid = 0.;

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();

            const Real* dx = geom[lev].CellSize();

            // FIXME - call estdt instead

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
/*
            firstdt(dt_grid,umax_grid,
                    ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                    ZFILL(dx),
                    BL_TO_FORTRAN_FAB(sold_mf[mfi]),
                    BL_TO_FORTRAN_FAB(uold_mf[mfi]));
*/

            dt = std::min(dt,dt_grid);
            umax = std::max(umax,umax_grid);
        }
    }

    // find the smallest dt over all processors
    ParallelDescriptor::ReduceRealMin(dt);

    // find the largest umax over all processors
    ParallelDescriptor::ReduceRealMax(umax);

    // set rel_eps in fortran module
    set_rel_eps(umax*1.e-8);        

    if (fixed_dt == -1.) {

        if (maestro_verbose > 0) {
            Print() << "Call to EstDt at beginning of step" << istep << endl;
            Print() << "gives dt = " << dt << endl;
        }

        if (nuclear_dt_fac < 1.0) {
            dt = nuclear_dt_fac*dtold;
            Print() << "Applying nuclear_dt_fac; new dt =" << dt << endl;
        }

        if (dt > max_dt_growth*dtold)
        {
            dt = max_dt_growth*dtold;
            if (maestro_verbose > 0) {
                Print() << "dt_growth factor limits the new dt = " << dt << endl;
            }
        }

        if (dt > max_dt) {
            dt = max_dt;
            if (maestro_verbose > 0) {
                Print() << "max_dt limits the new dt = " << max_dt << endl;
            }
        }

        if (dt < small_dt) {
            Abort("EstDt: dt < small_dt");
        }

    }
    else {

        // fixed dt
        dt = fixed_dt;
        if (maestro_verbose > 0) {
            Print() << "Setting fixed dt = " << dt << endl;
        }

    }

    // limit dt by stop_time
    if (t_new + dt > stop_time) {
        dt = stop_time - t_new;
        Print() << "Stop time limits dt = " << dt << endl;
    }
}

void
Maestro::FirstDt ()
{

    // allocate a dummy w0_force and set equal to zero
    Vector<Real> w0_force_dummy( (max_radial_level+1)*nr_fine );
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
    

    dt = 1.e99;
        
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
            firstdt(lev,dt_grid,umax_grid,
                    ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                    ZFILL(dx),
                    BL_TO_FORTRAN_FAB(sold_mf[mfi]),
                    BL_TO_FORTRAN_FAB(uold_mf[mfi]),
                    BL_TO_FORTRAN_FAB(vel_force_mf[mfi]),
                    BL_TO_FORTRAN_3D(S_cc_old_mf[mfi]),
                    p0_old.dataPtr(),
                    gamma1bar_old.dataPtr());

            dt = std::min(dt,dt_grid);
            umax = std::max(umax,umax_grid);
        }
    }

    // find the smallest dt over all processors
    ParallelDescriptor::ReduceRealMin(dt);

    // find the largest umax over all processors
    ParallelDescriptor::ReduceRealMax(umax);

    // set rel_eps in fortran module
    set_rel_eps(umax*1.e-8);    

    if (fixed_dt == -1.) {
   
        // use dt obtained with FirstDt
        if (maestro_verbose > 0) {
            Print() << "Call to FirstDt gives dt = " << dt << endl;
        }

        if (init_shrink != 1.0) {
            dt *= init_shrink;
            if (maestro_verbose > 0) {
                Print() << "Multiplying dt by init_shrink gives dt = " << dt << endl;
            }
        }

        // limit dt by max_dt
        if (dt > max_dt) {
            dt = max_dt;
            if (maestro_verbose > 0) {
                Print() << "max_dt limits the new dt = " << max_dt << endl;
            }
        }

        if (dt < small_dt) {
            Abort("FirstDt: dt < small_dt");
        }

    }
    else {

        // fixed dt
        dt = fixed_dt;
        if (maestro_verbose > 0) {
            Print() << "Setting fixed dt = " << dt << endl;
        }

    }

    // limit dt by stop_time
    if (t_new + dt > stop_time) {
        dt = stop_time - t_new;
        if (maestro_verbose > 0) {
            Print() << "Stop time limits dt = " << dt << endl;
        }
    }
}

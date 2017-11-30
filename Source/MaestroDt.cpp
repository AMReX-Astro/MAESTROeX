
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
        MultiFab& u_mf = unew[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(u_mf); mfi.isValid(); ++mfi ) {

            Real dt_grid = 1.e99;
            Real umax_grid = 0.;

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.



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

        if (verbose > 0) {
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
            if (verbose > 0) {
                Print() << "dt_growth factor limits the new dt = " << dt << endl;
            }
        }

        if (dt > max_dt) {
            dt = max_dt;
            if (verbose > 0) {
                Print() << "max_dt limits the new dt = " << max_dt << endl;
            }
        }
    }
    else {

        // fixed dt
        dt = fixed_dt;
        if (verbose > 0) {
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
    dt = 1.e99;
        
    Real umax = 0.;

    for (int lev = 0; lev <= finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& u_mf = unew[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(u_mf); mfi.isValid(); ++mfi ) {

            Real dt_grid = 1.e99;
            Real umax_grid = 0.;

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.



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
        if (verbose > 0) {
            Print() << "Call to FirstDt gives dt = " << dt << endl;
        }

        if (init_shrink != 1.0) {
            dt *= init_shrink;
            if (verbose > 0) {
                Print() << "Multiplying dt by init_shrink gives dt = " << dt << endl;
            }
        }

        // limit dt by max_dt
        if (dt > max_dt) {
            dt = max_dt;
            if (verbose > 0) {
                Print() << "max_dt limits the new dt = " << max_dt << endl;
            }
        }

    }
    else {

        // fixed dt
        dt = fixed_dt;
        if (verbose > 0) {
            Print() << "Setting fixed dt = " << dt << endl;
        }

    }

    // limit dt by stop_time
    if (t_new + dt > stop_time) {
        dt = stop_time - t_new;
        if (verbose > 0) {
            Print() << "Stop time limits dt = " << dt << endl;
        }
    }
}

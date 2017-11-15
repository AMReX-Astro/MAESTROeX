
#include <Maestro.H>

using namespace amrex;

// a wrapper for EstDtLevel
void
Maestro::EstDt ()
{

    // initial timestepping is handled in FirstDt
    if (istep == 1) return;

    if (fixed_dt == -1.0) {

        // compute dt at each level from CFL considerations
        Vector<Real> dt_tmp(finest_level+1);
        for (int lev = 0; lev <= finest_level; ++lev) {
            dt_tmp[lev] = EstDtLevel(lev);
        }

        // select the smallest time step over all levels
        dt = dt_tmp[0];
        for (int lev = 1; lev <= finest_level; ++lev) {
            dt = std::min(dt, dt_tmp[lev]);
        }

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

        dt = fixed_dt;
        if (verbose > 0) {
            Print() << "Setting fixed dt = " << dt << endl;
        }

    }

    // Limit dt's by the value of stop_time.
    if (t_new + dt > stop_time) {
        dt = stop_time - t_new;
        Print() << "Stop time limits dt = " << dt << endl;
    }
}

// compute dt at a level from CFL considerations
Real
Maestro::EstDtLevel (int lev) const
{
    BL_PROFILE("Maestro::EstTimeStep()");

    // FIXME
    return fixed_dt;

}

void
Maestro::FirstDt ()
{
    // FIXME

}

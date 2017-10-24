
#include <Maestro.H>

using namespace amrex;

// a wrapper for ComputeDtLevel
void
Maestro::ComputeDt ()
{

    // compute dt at each level from CFL considerations
    Vector<Real> dt_tmp(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        dt_tmp[lev] = ComputeDtLevel(lev);
    }

    // select the smallest time step over all levels
    Real dt_0 = dt_tmp[0];
    for (int lev = 1; lev <= finest_level; ++lev) {
        dt_0 = std::min(dt_0, dt_tmp[lev]);
    }

    // do not allow time step to increase too much from previous
    constexpr Real change_max = 1.1;
    if (dt_0 > change_max*dt)
    {
        Print() << "Reducing time step to respect change_max" << std::endl;
        dt_0 = change_max*dt;
    }

    // Limit dt's by the value of stop_time.
    const Real eps = 1.e-3*dt_0;
    if (t_new + dt_0 > stop_time - eps) {
        Print() << "\nModifying time step to respect stop_time" << std::endl;
        Print() << "Original dt = " << dt_0 << std::endl;
        Print() << "New dt = " << stop_time-t_new << std::endl;
        dt_0 = stop_time - t_new;
    }

    dt = dt_0;
}

// compute dt at a level from CFL considerations
Real
Maestro::ComputeDtLevel (int lev) const
{
    BL_PROFILE("Maestro::EstTimeStep()");

    Real dt_est = std::numeric_limits<Real>::max();

    const Real* dx = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();
    const Real cur_time = t_new;
    const MultiFab& S_new = *snew[lev];

#ifdef _OPENMP
#pragma omp parallel reduction(min:dt_est)
#endif
    {
        FArrayBox uface[BL_SPACEDIM];

        for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
        {
            for (int i = 0; i < BL_SPACEDIM ; i++) {
                const Box& bx = mfi.nodaltilebox(i);
                uface[i].resize(bx,1);
            }

            get_face_velocity(lev, cur_time,
                              AMREX_D_DECL(BL_TO_FORTRAN(uface[0]),
                                           BL_TO_FORTRAN(uface[1]),
                                           BL_TO_FORTRAN(uface[2])),
                              dx, prob_lo);

            for (int i = 0; i < BL_SPACEDIM; ++i) {
                Real umax = uface[i].norm(0);
                if (umax > 1.e-100) {
                    dt_est = std::min(dt_est, dx[i] / umax);
                }
            }
        }
    }

    ParallelDescriptor::ReduceRealMin(dt_est);

    dt_est *= cflfac;

    return dt_est;
}

#include <Maestro.H>
#include <Postprocess.H>
#include "AMReX_PlotFileUtil.H"

using namespace amrex;

// ------------------------------------------
// Write 2D r-theta slice diagnostics
// ------------------------------------------
void Postprocess::Write2dSliceFile(const Vector<MultiFab>& u_in,
                                   const Vector<MultiFab>& w0_in,
                                   const int deltat, const int nfiles) {
    // timer for profiling
    BL_PROFILE_VAR("Postprocess::Write2dSliceFile()", Write2dSliceFile);

    // MakeMeridionalCirculation
    int max_grid_size = GetMaxGridSize(iFile);
    Geometry geomFine;
    Vector<MultiFab> meridion_vel(1);
    MakeMeridionalCirculation(u_in, w0_in, meridion_vel, max_grid_size,
                              geomFine);

    // time-averaged circulation need multiple files
    if (deltat > 0 && nfiles > 0) {
        amrex::PlotFileData pltfile(iFile);

        // read in base file name to get time step
        std::string basefilename =
            GetVarFromJobInfo(iFile, "maestro.plot_base_name");
        int tn = GetTimeStep(iFile, basefilename);
        // Print() << "basefilename = " << basefilename << ", timestep = " << tn << std::endl;

        // save initial time
        Real t0 = pltfile.time();

        // save old parameters and outputs
        Vector<MultiFab> u_old(finest_level + 1);
        Vector<MultiFab> w0_old(finest_level + 1);
        for (int lev = 0; lev <= finest_level; ++lev) {
            u_old[lev].define(grid[lev], dmap[lev], AMREX_SPACEDIM, 0);
            w0_old[lev].define(grid[lev], dmap[lev], AMREX_SPACEDIM, 0);

            MultiFab::Copy(u_old[lev], u_in[lev], 0, 0, AMREX_SPACEDIM, 0);
            MultiFab::Copy(w0_old[lev], w0_in[lev], 0, 0, AMREX_SPACEDIM, 0);
        }

        int ncomp = meridion_vel[0].nComp();
        Vector<MultiFab> meridion_vel_old(1);
        meridion_vel_old[0].define(meridion_vel[0].boxArray(),
                                   meridion_vel[0].DistributionMap(), ncomp, 0);
        MultiFab::Copy(meridion_vel_old[0], meridion_vel[0], 0, 0, ncomp, 0);
        meridion_vel[0].setVal(0.0);
        Real t_old = t0;

        // new plotfile name and time
        std::string plotfilename_new = basefilename;
        Real t_new;

        // define new velocities
        Vector<MultiFab> u_new(finest_level + 1);
        Vector<MultiFab> w0_new(finest_level + 1);
        for (int lev = 0; lev <= finest_level; ++lev) {
            u_new[lev].define(grid[lev], dmap[lev], AMREX_SPACEDIM, 0);
            w0_new[lev].define(grid[lev], dmap[lev], AMREX_SPACEDIM, 0);
        }

        for (int n = 0; n < nfiles; ++n) {
            // update timestep and plotfile name
            tn += deltat;
            PlotFileName(tn, &basefilename, &plotfilename_new);

            // open file
            amrex::PlotFileData pltfile_new(plotfilename_new);

            // get new time
            t_new = pltfile_new.time();

            // read in new velocities
            for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                std::string x = "vel";
                std::string w = "w0";
                x += (120 + i);
                w += (120 + i);

                for (int lev = 0; lev <= finest_level; ++lev) {
                    MultiFab::Copy(u_new[lev], pltfile_new.get(lev, x), 0, i, 1,
                                   0);
                    MultiFab::Copy(w0_new[lev], pltfile_new.get(lev, w), 0, i,
                                   1, 0);
                }
            }

            // MakeMeridionalCirculation
            Vector<MultiFab> meridion_vel_new(1);
            MakeMeridionalCirculation(u_new, w0_new, meridion_vel_new,
                                      max_grid_size, geomFine);

            // vel = vel + dt/2 * (vel_old + vel_new)
            Real dt = t_new - t_old;
            MultiFab::Add(meridion_vel_old[0], meridion_vel_new[0], 0, 0, ncomp,
                          0);
            MultiFab::Saxpy(meridion_vel[0], 0.5 * dt, meridion_vel_old[0], 0,
                            0, ncomp, 0);

            // set old states to new states values
            for (int lev = 0; lev <= finest_level; ++lev) {
                MultiFab::Copy(u_old[lev], u_new[lev], 0, 0, AMREX_SPACEDIM, 0);
                MultiFab::Copy(w0_old[lev], w0_new[lev], 0, 0, AMREX_SPACEDIM,
                               0);
            }
            MultiFab::Copy(meridion_vel_old[0], meridion_vel_new[0], 0, 0,
                           ncomp, 0);
            t_old = t_new;
        }

        // divide vel by total time elapsed
        const Real ttot = t_new - t0;
        int lev = 0;

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(meridion_vel[lev], TilingIfNotGPU()); mfi.isValid();
             ++mfi) {
            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const Array4<Real> vel_arr = meridion_vel[lev].array(mfi);

            AMREX_PARALLEL_FOR_3D(tileBox, i, j, k,
                                  { vel_arr(i, j, k) /= ttot; });
        }
    }

    // write to disk
    std::string slicefilename = "slice_" + iFile;

    WriteSingleLevelPlotfile(slicefilename, meridion_vel[0],
                             {"vel_radial", "vel_theta"}, geomFine, 0, 0);
}

// get plotfile name
void Postprocess::PlotFileName(const int timestep, std::string* basefilename,
                               std::string* plotfilename) {
    *plotfilename = Concatenate(*basefilename, timestep, 7);
}

void Postprocess::MakeMeridionalCirculation(const Vector<MultiFab>& vel,
                                            const Vector<MultiFab>& w0rcart,
                                            Vector<MultiFab>& circ_vel,
                                            const int max_grid_size,
                                            Geometry& geomFine) {
    // timer for profiling
    BL_PROFILE_VAR("Postprocess::MakeMeridionalCirculation()",
                   MakeMeridionalCirculation);

    // make circulation velocities
    Vector<MultiFab> rad_vel(finest_level + 1);
    Vector<MultiFab> ang_vel(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        rad_vel[lev].define(vel[lev].boxArray(), vel[lev].DistributionMap(), 1,
                            0);
        ang_vel[lev].define(vel[lev].boxArray(), vel[lev].DistributionMap(), 1,
                            0);
    }

    const auto& center_p = center;

    for (int lev = finest_level; lev <= finest_level; ++lev) {
        const auto dx = pgeom[lev].CellSizeArray();
        const auto prob_lo = pgeom[lev].ProbLoArray();

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(vel[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const Array4<const Real> vel_arr = vel[lev].array(mfi);
            const Array4<Real> radvel_arr = rad_vel[lev].array(mfi);
            const Array4<Real> angvel_arr = ang_vel[lev].array(mfi);
            const Array4<const Real> w0rcart_arr = w0rcart[lev].array(mfi);

            AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                radvel_arr(i, j, k) = 0.0;
                angvel_arr(i, j, k) = 0.0;

                Real x = prob_lo[0] + (Real(i) + 0.5) * dx[0] - center_p[0];
                Real y = prob_lo[1] + (Real(j) + 0.5) * dx[1] - center_p[1];
                Real z = prob_lo[2] + (Real(k) + 0.5) * dx[2] - center_p[2];
                Real inv_radius = 1.0 / sqrt(x * x + y * y + z * z);
                Real inv_xy = 1.0 / sqrt(x * x + y * y);

                Vector<Real> normal(3);
                normal[0] = x * inv_radius;
                normal[1] = y * inv_radius;
                normal[2] = z * inv_radius;

                for (auto n = 0; n < AMREX_SPACEDIM; ++n) {
                    radvel_arr(i, j, k) += vel_arr(i, j, k, n) * normal[n];
                }

                Vector<Real> theta_dir(3);
                theta_dir[0] = x * inv_radius * z * inv_xy;
                theta_dir[1] = y * inv_radius * z * inv_xy;
                theta_dir[2] = -inv_radius / inv_xy;

                for (auto n = 0; n < AMREX_SPACEDIM; ++n) {
                    angvel_arr(i, j, k) += vel_arr(i, j, k, n) * theta_dir[n];
                }

                // add base state vel to get full radial velocity
                radvel_arr(i, j, k) += w0rcart_arr(i, j, k);
            });
        }
    }

    // need single-level finest grid for whole domain
    //  TODO: domain size should be larger than simply half x-z slice
    //        since its diagonal plane has the largest area
    const Box& domain0 = pgeom[0].Domain();
    auto dom0_lo = domain0.loVect();
    auto dom0_hi = domain0.hiVect();
    IntVect ncell(0, 1, 0);
    ncell[0] = (dom0_hi[0] + 1) * (finest_level + 1) / 2 - 1;
    ncell[2] = (dom0_hi[2] + 1) * (finest_level + 1) - 1;

    IntVect domLo(AMREX_D_DECL(dom0_lo[0], dom0_lo[1], dom0_lo[2]));
    IntVect domHi(AMREX_D_DECL(ncell[0], ncell[1], ncell[2]));
    Box domain(domLo, domHi);

    // we only need half x-z plane because theta is between (0,pi)
    auto prob_lo = pgeom[0].ProbLo();
    auto prob_hi = pgeom[0].ProbHi();
    const auto dx = pgeom[finest_level].CellSizeArray();

    Real probLo[] = {prob_lo[0] + (prob_hi[0] - prob_lo[0]) / 2.0, prob_lo[1],
                     prob_lo[2]};
    Real probHi[] = {prob_hi[0], prob_lo[1] + 2.0 * dx[1], prob_hi[2]};
    RealBox real_box(probLo, probHi);

    // make BoxArray and Geometry
    BoxArray ba(domain);
    ba.maxSize(max_grid_size);
    DistributionMapping dm(ba);

    geomFine.define(domain, &real_box);

    // define velocities on new single-level grid
    int numcomp = 2;
    circ_vel[0].define(ba, dm, numcomp, 0);

    // average to 2D r-theta (x-z) plane
    Average2d(rad_vel, circ_vel, 0, {domain}, 0);
    Average2d(ang_vel, circ_vel, 1, {domain}, 0);
}


#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

// advance solution to final time
void Maestro::Evolve() {
    Vector<MultiFab> dens_orig(finest_level + 1);
    Vector<MultiFab> dens_final(finest_level + 1);
    Vector<MultiFab> error(finest_level + 1);
    Vector<MultiFab> scal_force(finest_level + 1);
    Vector<MultiFab> etarhoflux(finest_level + 1);
    Vector<std::array<MultiFab, AMREX_SPACEDIM> > umac(finest_level + 1);
    Vector<std::array<MultiFab, AMREX_SPACEDIM> > sedge(finest_level + 1);
    Vector<std::array<MultiFab, AMREX_SPACEDIM> > sflux(finest_level + 1);
    Vector<std::array<MultiFab, AMREX_SPACEDIM> > w0mac(finest_level + 1);
    BaseState<Real> rho0_predicted_edge(base_geom.max_radial_level + 1,
                                        base_geom.nr_fine + 1);

    Vector<Real> abs_norm(finest_level + 1);
    Vector<Real> rel_norm(finest_level + 1);

    for (int lev = 0; lev <= finest_level; ++lev) {
        dens_orig[lev].define(grids[lev], dmap[lev], 1, ng_s);
        dens_orig[lev].setVal(0.);
        dens_final[lev].define(grids[lev], dmap[lev], 1, ng_s);
        dens_final[lev].setVal(0.);
        error[lev].define(grids[lev], dmap[lev], 1, ng_s);
        error[lev].setVal(0.);
        scal_force[lev].define(grids[lev], dmap[lev], Nscal, 1);
        scal_force[lev].setVal(0.);

        AMREX_D_TERM(etarhoflux[lev].define(convert(grids[lev], nodal_flag_x),
                                            dmap[lev], 1, 1);
                     , etarhoflux[lev].define(convert(grids[lev], nodal_flag_y),
                                              dmap[lev], 1, 1);
                     , etarhoflux[lev].define(convert(grids[lev], nodal_flag_z),
                                              dmap[lev], 1, 1););

        // face-centered arrays of MultiFabs
        AMREX_D_TERM(umac[lev][0].define(convert(grids[lev], nodal_flag_x),
                                         dmap[lev], 1, 1);
                     , umac[lev][1].define(convert(grids[lev], nodal_flag_y),
                                           dmap[lev], 1, 1);
                     , umac[lev][2].define(convert(grids[lev], nodal_flag_z),
                                           dmap[lev], 1, 1););

        AMREX_D_TERM(sedge[lev][0].define(convert(grids[lev], nodal_flag_x),
                                          dmap[lev], Nscal, 0);
                     , sedge[lev][1].define(convert(grids[lev], nodal_flag_y),
                                            dmap[lev], Nscal, 0);
                     , sedge[lev][2].define(convert(grids[lev], nodal_flag_z),
                                            dmap[lev], Nscal, 0););

        AMREX_D_TERM(sflux[lev][0].define(convert(grids[lev], nodal_flag_x),
                                          dmap[lev], Nscal, 0);
                     , sflux[lev][1].define(convert(grids[lev], nodal_flag_y),
                                            dmap[lev], Nscal, 0);
                     , sflux[lev][2].define(convert(grids[lev], nodal_flag_z),
                                            dmap[lev], Nscal, 0););

        etarhoflux[lev].setVal(0.);
    }

#if (AMREX_SPACEDIM == 3)
    for (int lev = 0; lev <= finest_level; ++lev) {
        w0mac[lev][0].define(convert(grids[lev], nodal_flag_x), dmap[lev], 1,
                             1);
        w0mac[lev][1].define(convert(grids[lev], nodal_flag_y), dmap[lev], 1,
                             1);
        w0mac[lev][2].define(convert(grids[lev], nodal_flag_z), dmap[lev], 1,
                             1);

        for (int d = 0; d < AMREX_SPACEDIM; ++d) w0mac[lev][d].setVal(0.);
    }
#endif

    // -------------------------------------------------------------------------
    //  loop over all possible directions
    // -------------------------------------------------------------------------

    for (auto i = 0; i < AMREX_SPACEDIM; i++) {
        for (auto j = -1; j <= 1; j += 2) {
            t_old = 0.0;
            t_new = 0.0;
            istep = 1;

            // reset the density
            for (int lev = 0; lev <= finest_level; ++lev) {
                MakeNewLevelFromScratch(lev, t_old, grids[lev], dmap[lev]);
            }

            // reset tagging array to include buffer zones
            TagArray();

            // compute numdisjointchunks, r_start_coord, r_end_coord
            BaseState<int> tag_array_b(
                tag_array, base_geom.max_radial_level + 1, base_geom.nr_fine);
            base_geom.InitMultiLevel(finest_level, tag_array_b.array());

            // average down data and fill ghost cells
            AverageDown(sold, 0, Nscal);
            FillPatch(t_old, sold, sold, sold, 0, 0, Nscal, 0, bcs_s);

            // first compute cutoff coordinates using initial density profile
            ComputeCutoffCoords(rho0_old);
            base_geom.ComputeCutoffCoords(rho0_old.array());

            // set rho0 to be the average
            Average(sold, rho0_old, Rho);
            ComputeCutoffCoords(rho0_old);
            base_geom.ComputeCutoffCoords(rho0_old.array());

            // call eos with r,p as input to recompute T,h
            TfromRhoP(sold, p0_old, 1);

            // set rhoh0 to be the average
            Average(sold, rhoh0_old, RhoH);

            Print() << "\nInitial velocity = ";
            if (i == 0) {
                Print() << "(" << j << ", 0)" << std::endl;
            } else {
                Print() << "(0, " << j << ")" << std::endl;
            }

            // the base state will not carry any information in this test problem
            rho0_old.setVal(0.0);
            rho0_new.setVal(0.0);
            rhoh0_old.setVal(0.0);
            rhoh0_new.setVal(0.0);
            p0_old.setVal(0.0);
            p0_new.setVal(0.0);
            w0.setVal(0.0);
            rho0_predicted_edge.setVal(0.);

            // initialize the velocity field -- it is unity in the
            // direction of propagation, a negative itest_dir indicates
            // negative velocity
            for (int lev = 0; lev <= finest_level; ++lev) {
                for (int d = 0; d < AMREX_SPACEDIM; ++d)
                    umac[lev][d].setVal(0.0);

                umac[lev][i].setVal(double(j));
            }

            if (finest_level == 0) {
                // fill periodic ghost cells
                for (int lev = 0; lev <= finest_level; ++lev) {
                    for (int d = 0; d < AMREX_SPACEDIM; ++d)
                        umac[lev][d].FillBoundary(geom[lev].periodicity());
                }
                // fill ghost cells behind physical boundaries
                FillUmacGhost(umac);
            } else {
                // edge_restriction for velocities
                AverageDownFaces(umac);
                // fill level n ghost cells using interpolation from level n-1 data
                FillPatchUedge(umac);
            }

            // Store the initial density here
            for (int lev = 0; lev <= finest_level; ++lev)
                MultiFab::Copy(dens_orig[lev], sold[lev], Rho, 0, 1, 0);

            Print() << "original density = " << dens_orig[0].norm2()
                    << std::endl;

            // compute the initial timestep -- dt = dx / u, where u = 1
            const Real* dx = geom[finest_level].CellSize();
            dt = cfl * dx[0];

            // advance the density using the constant velocity field
            for (istep = start_step; istep <= max_step && t_old < stop_time;
                 ++istep) {
                t_old = t_new;

                DensityAdvance(1, sold, snew, sedge, sflux, scal_force,
                               etarhoflux, umac, w0mac, rho0_predicted_edge);

                // move new state into old state by swapping pointers
                for (int lev = 0; lev <= finest_level; ++lev) {
                    std::swap(sold[lev], snew[lev]);
                    std::swap(rho0_old, rho0_new);
                }

                t_new = t_old + dt;

                if (t_new + dt > stop_time) dt = stop_time - t_new;
            }

            // Store the final density here
            for (int lev = 0; lev <= finest_level; ++lev)
                MultiFab::Copy(dens_final[lev], snew[lev], Rho, 0, 1, 0);

            Print() << "final density = " << dens_final[0].norm2() << std::endl;

            // compare the initial and final density
            // compute dens_final - dens_orig
            for (int lev = 0; lev <= finest_level; ++lev) {
                MultiFab::Copy(error[lev], dens_final[lev], 0, 0, 1, 0);
                MultiFab::Subtract(error[lev], dens_orig[lev], 0, 0, 1, 0);

                abs_norm[lev] = error[lev].norm2();

                rel_norm[lev] = error[lev].norm2() / dens_orig[lev].norm2();

                Print() << "\tAbs norm = " << abs_norm[lev]
                        << "  Rel norm = " << rel_norm[lev] << std::endl;
            }
        }
    }
}

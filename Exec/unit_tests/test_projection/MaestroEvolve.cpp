
#include <Maestro.H>
#include <extern_parameters.H>

using namespace amrex;

// advance solution to final time
void Maestro::Evolve() {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Evolve()", Evolve);

    Print() << "Calling Evolve()" << std::endl;

    // -------------------------------------------------------------------------
    //  allocate arrays
    // -------------------------------------------------------------------------

    Print() << "...allocate arrays" << std::endl;

    Vector<MultiFab> umid(finest_level + 1);
    Vector<MultiFab> gphi(finest_level + 1);
    Vector<std::array<MultiFab, AMREX_SPACEDIM> > umac_old(finest_level + 1);
    Vector<std::array<MultiFab, AMREX_SPACEDIM> > umac_mid(finest_level + 1);
    Vector<std::array<MultiFab, AMREX_SPACEDIM> > umac_new(finest_level + 1);
    Vector<std::array<MultiFab, AMREX_SPACEDIM> > gphi_mac(finest_level + 1);
    Vector<MultiFab> utemp(finest_level + 1);

    if (project_type == 1) {
        // HG projection.  Velocities are cell-centered
        for (int lev = 0; lev <= finest_level; ++lev) {
            umid[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, ng_s);
            gphi[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, ng_s);

            umid[lev].setVal(0.);
            gphi[lev].setVal(0.);
        }
    } else {
        // MAC projection.  Velocities are nodal in respective dimension
        for (int lev = 0; lev <= finest_level; ++lev) {
            AMREX_D_TERM(
                umac_old[lev][0].define(convert(grids[lev], nodal_flag_x),
                                        dmap[lev], 1, ng_s);
                , umac_old[lev][1].define(convert(grids[lev], nodal_flag_y),
                                          dmap[lev], 1, ng_s);
                , umac_old[lev][2].define(convert(grids[lev], nodal_flag_z),
                                          dmap[lev], 1, ng_s););
            AMREX_D_TERM(
                umac_mid[lev][0].define(convert(grids[lev], nodal_flag_x),
                                        dmap[lev], 1, ng_s);
                , umac_mid[lev][1].define(convert(grids[lev], nodal_flag_y),
                                          dmap[lev], 1, ng_s);
                , umac_mid[lev][2].define(convert(grids[lev], nodal_flag_z),
                                          dmap[lev], 1, ng_s););
            AMREX_D_TERM(
                umac_new[lev][0].define(convert(grids[lev], nodal_flag_x),
                                        dmap[lev], 1, ng_s);
                , umac_new[lev][1].define(convert(grids[lev], nodal_flag_y),
                                          dmap[lev], 1, ng_s);
                , umac_new[lev][2].define(convert(grids[lev], nodal_flag_z),
                                          dmap[lev], 1, ng_s););
            AMREX_D_TERM(
                gphi_mac[lev][0].define(convert(grids[lev], nodal_flag_x),
                                        dmap[lev], 1, ng_s);
                , gphi_mac[lev][1].define(convert(grids[lev], nodal_flag_y),
                                          dmap[lev], 1, ng_s);
                , gphi_mac[lev][2].define(convert(grids[lev], nodal_flag_z),
                                          dmap[lev], 1, ng_s););

            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                umac_old[lev][d].setVal(0.);
                umac_mid[lev][d].setVal(0.);
                umac_new[lev][d].setVal(0.);
                gphi_mac[lev][d].setVal(0.);
            }
            utemp[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, ng_s);
            utemp[lev].setVal(0.);
        }
    }

    // -------------------------------------------------------------------------
    //  initialize velocity field
    // -------------------------------------------------------------------------

    Print() << "...initialize velocity field" << std::endl;

    if (project_type == 1) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab& vel = uold[lev];
            const Real* dx = geom[lev].CellSize();
            const auto prob_lo = geom[lev].ProbLoArray();

            // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(vel, true); mfi.isValid(); ++mfi) {
                const Box& tilebox = mfi.tilebox();

                Array4<Real> const vel_arr = vel.array(mfi);

                amrex::ParallelFor(tilebox,
                [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
                {

                    Real x = (static_cast<Real>(i)+0.5_rt) * dx[0] + prob_lo[0];
                    Real y = (static_cast<Real>(j)+0.5_rt) * dx[1] + prob_lo[1];
#if (AMREX_SPACEDIM == 3)
                    Real z = (static_cast<Real>(k)+0.5_rt) * dx[2] + prob_lo[2];
#endif

#if (AMREX_SPACEDIM==2)
                    vel_arr(i,j,k,0) = -std::pow(std::sin(M_PI*x), 2) * std::sin(2.0_rt*M_PI*y);
                    vel_arr(i,j,k,1) = std::pow(std::sin(M_PI*y), 2) * std::sin(2.0_rt*M_PI*x);

#else
                    vel_arr(i,j,k,0) =
                        2.0_rt * M_PI * std::sin(4.0_rt*M_PI*x) * std::cos(2.0_rt*M_PI*y) -
                        4.0_rt * M_PI * std::sin(2.0_rt*M_PI*x) * std::cos(4.0_rt*M_PI*z);

                    vel_arr(i,j,k,1) =
                        2.0_rt * M_PI * std::sin(4.0_rt*M_PI*y) * std::cos(2.0_rt*M_PI*z) -
                        4.0_rt * M_PI * std::cos(4.0_rt*M_PI*x) * std::sin(2.0_rt*M_PI*y);

                    vel_arr(i,j,k,2) =
                        2.0_rt * M_PI * std::cos(2.0_rt*M_PI*x) * std::sin(4.0_rt*M_PI*z) -
                        4.0_rt * M_PI * std::cos(4.0_rt*M_PI*y) * std::sin(2.0_rt*M_PI*z);
#endif
                });
            }

        }

        AverageDown(uold, 0, AMREX_SPACEDIM);
        FillPatch(t_old, uold, uold, uold, 0, 0, AMREX_SPACEDIM, 0, bcs_u);

    } else {
        // need to initialize the mac velocity
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab& umac_mf = umac_old[lev][0];
            MultiFab& vmac_mf = umac_old[lev][1];
            MultiFab& wmac_mf = umac_old[lev][2];
            const Real* dx = geom[lev].CellSize();
            const auto prob_lo = geom[lev].ProbLoArray();

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(umac_mf, true); mfi.isValid(); ++mfi) {

                Array4<Real> const umac_arr = umac_mf.array(mfi);
                Array4<Real> const vmac_arr = vmac_mf.array(mfi);
#if AMREX_SPACEDIM == 3
                Array4<Real> const wmac_arr = wmac_mf.array(mfi);
#endif

                // x-velocity  (x are edges, y and z are centers)

                amrex::ParallelFor(mfi.nodaltilebox(0),
                [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
                {

                    Real x = static_cast<Real>(i) * dx[0] + prob_lo[0];
                    Real y = (static_cast<Real>(j)+0.5_rt) * dx[1] + prob_lo[1];
#if (AMREX_SPACEDIM==3)
                    Real z = (static_cast<Real>(k)+0.5_rt) * dx[2] + prob_lo[2];
#endif

#if (AMREX_SPACEDIM==2)
                    umac_arr(i,j,k) = -std::pow(std::sin(M_PI*x), 2) * std::sin(2.0_rt*M_PI*y);
#else
                    umac_arr(i,j,k) =
                        2.0_rt * M_PI * std::sin(4.0_rt*M_PI*x) * std::cos(2.0_rt*M_PI*y) -
                        4.0_rt * M_PI * std::sin(2.0_rt*M_PI*x) * std::cos(4.0_rt*M_PI*z);
#endif
                });

                // y-velocity  (x and z are centers, y are edges)

                amrex::ParallelFor(mfi.nodaltilebox(1),
                [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
                {

                    Real x = (static_cast<Real>(i)+0.5_rt) * dx[0] + prob_lo[0];
                    Real y = static_cast<Real>(j) * dx[1] + prob_lo[1];
#if (AMREX_SPACEDIM==3)
                    Real z = (static_cast<Real>(k)+0.5_rt) * dx[2] + prob_lo[2];
#endif

#if AMREX_SPACEDIM == 2
                    vmac_arr(i,j,k) = std::pow(std::sin(M_PI*y), 2) * std::sin(2.0_rt*M_PI*x);
#else

                    vmac_arr(i,j,k) =
                        2.0_rt * M_PI * std::sin(4.0_rt*M_PI*y) * std::cos(2.0_rt*M_PI*z) -
                        4.0_rt * M_PI * std::cos(4.0_rt*M_PI*x) * std::sin(2.0_rt*M_PI*y);
#endif
                });

#if (AMREX_SPACEDIM == 3)
                // z-velocity  (x and y are centers, z are edges)

                amrex::ParallelFor(mfi.nodaltilebox(2),
                [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
                {

                    Real x = (static_cast<Real>(i)+0.5_rt) * dx[0] + prob_lo[0];
                    Real y = (static_cast<Real>(j)+0.5_rt) * dx[1] + prob_lo[1];
                    Real z = static_cast<Real>(k) * dx[2] + prob_lo[2];

                    wmac_arr(i,j,k) =
                        2.0_rt * M_PI * std::cos(2.0_rt*M_PI*x) * std::sin(4.0_rt*M_PI*z) -
                        4.0_rt * M_PI * std::cos(4.0_rt*M_PI*y) * std::sin(2.0_rt*M_PI*z);

                });
#endif

            }
        }

        if (finest_level == 0) {
            // fill periodic ghost cells
            for (int lev = 0; lev <= finest_level; ++lev) {
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    umac_old[lev][d].FillBoundary(geom[lev].periodicity());
                }
            }
            // fill ghost cells behind physical boundaries
            FillUmacGhost(umac_old);
        } else {
            // edge_restriction for velocities
            AverageDownFaces(umac_old);
            // fill level n ghost cells using interpolation from level n-1 data
            FillPatchUedge(umac_old);
        }
    }

    BaseState<Real> dummy;

    if (project_type == 1) {
        // write uold
        WritePlotFile(0, t_new, dt, dummy, dummy, dummy, dummy, uold, uold,
                      uold);

        // copy the velocity field over to the intermediate state, umid
        for (int lev = 0; lev <= finest_level; ++lev)
            MultiFab::Copy(umid[lev], uold[lev], 0, 0, AMREX_SPACEDIM, ng_s);
    } else {
        // convect MAC field to cell-centered

        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab& umac_mf = umac_old[lev][0];
            MultiFab& vmac_mf = umac_old[lev][1];
            MultiFab& wmac_mf = umac_old[lev][2];
            MultiFab& utemp_mf = utemp[lev];

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(umac_mf, true); mfi.isValid(); ++mfi) {
                const Box& tilebox = mfi.tilebox();

                auto utemp_arr = utemp_mf.array(mfi);
                auto umac_arr = umac_mf.array(mfi);
                auto vmac_arr = vmac_mf.array(mfi);
#if AMREX_SPACEDIM == 3
                auto wmac_arr = wmac_mf.array(mfi);
#endif

                amrex::ParallelFor(tilebox,
                [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
                {

                    utemp_arr(i,j,k,0) = 0.5_rt * (umac_arr(i,j,k) + umac_arr(i+1,j,k));
                    utemp_arr(i,j,k,1) = 0.5_rt * (vmac_arr(i,j,k) + vmac_arr(i,j+1,k));
#if (AMREX_SPACEDIM==3)
                    utemp_arr(i,j,k,2) = 0.5_rt * (wmac_arr(i,j,k) + wmac_arr(i,j,k+1));
#endif

                });
            }
        }

        // write umac_old
        WritePlotFile(0, t_new, dt, dummy, dummy, dummy, dummy, utemp, uold,
                      uold);

        // copy the velocity field over to the intermediate state, umid
        for (int lev = 0; lev <= finest_level; ++lev) {
            for (auto comp = 0; comp < AMREX_SPACEDIM; ++comp)
                MultiFab::Copy(umac_mid[lev][comp], umac_old[lev][comp], 0, 0,
                               1, ng_s);
        }
    }

    //--------------------------------------------------------------------------
    // 'pollute' the velocity field by adding the gradient of a scalar
    //--------------------------------------------------------------------------

    Print() << "...pollute velocity field" << std::endl;

    if (project_type == 1) {
        // add_grad_scalar

        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab& gphi_mf = gphi[lev];
            MultiFab& umid_mf = umid[lev];
            const Real* dx = geom[lev].CellSize();
            const auto prob_lo = geom[lev].ProbLoArray();

#ifdef _OPENMP
#pragma omp parallel
#endif
            FArrayBox phi;

            for (MFIter mfi(gphi_mf, true); mfi.isValid(); ++mfi) {
                const Box& tilebox = mfi.tilebox();
                const Box& gbox = amrex::grow(tilebox, 1);

                if (   phys_bc[0] == SlipWall
                    && phys_bc[1] == SlipWall
#if AMREX_SPACEDIM == 3
                    && phys_bc[2] == SlipWall
#endif
                    && phys_bc[AMREX_SPACEDIM] == SlipWall
                    && phys_bc[AMREX_SPACEDIM+1] == SlipWall
#if AMREX_SPACEDIM == 3
                    && phys_bc[AMREX_SPACEDIM+2] == SlipWall
#endif
                    ) {

                    Array4<Real> const gphi_arr = gphi_mf.array(mfi);
                    Array4<Real> const umid_arr = umid_mf.array(mfi);

                    amrex::ParallelFor(tilebox,
                    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
                    {

                        Real x = (static_cast<Real>(i)+0.5_rt) * dx[0] + prob_lo[0];
                        Real y = (static_cast<Real>(j)+0.5_rt) * dx[1] + prob_lo[1];
#if AMREX_SPACEDIM == 3
                        Real z = (static_cast<Real>(k)+0.5_rt) * dx[2] + prob_lo[2];
#endif

#if AMREX_SPACEDIM == 2
                        gphi_arr(i,j,k,0) = 4.0_rt * x * (1.0_rt - x);
                        gphi_arr(i,j,k,1) = 4.0_rt * y * (1.0_rt - y);
#else
                        gphi_arr(i,j,k,0) = 160.0_rt * x * (1.0_rt - x);
                        gphi_arr(i,j,k,1) = 160.0_rt * y * (1.0_rt - y);
                        gphi_arr(i,j,k,2) = 160.0_rt * z * (1.0_rt - z);
#endif

                        umid_arr(i,j,k,0) += gphi_arr(i,j,k,0);
                        umid_arr(i,j,k,1) += gphi_arr(i,j,k,1);
#if AMREX_SPACEDIM == 3
                        umid_arr(i,j,k,2) += gphi_arr(i,j,k,2);
#endif

                        });

                } else if (    phys_bc[0] == Interior
                            && phys_bc[1] == Interior
#if AMREX_SPACEDIM == 3
                            && phys_bc[2] == Interior
#endif
                            && phys_bc[AMREX_SPACEDIM] == Interior
                            && phys_bc[AMREX_SPACEDIM+1] == Interior
#if AMREX_SPACEDIM == 3
                            && phys_bc[AMREX_SPACEDIM+2] == Interior
#endif
                ) {

                    // first fill phi on a local FAB on a grown tile box

                    phi.resize(gbox, 1);
                    Elixir elix_phi = phi.elixir();
                    auto phi_arr = phi.array();

                    // now we can explicitly difference phi

                    amrex::ParallelFor(gbox,
                    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
                    {

                        Real x = (static_cast<Real>(i)+0.5_rt) * dx[0] + prob_lo[0];
                        Real y = (static_cast<Real>(j)+0.5_rt) * dx[1] + prob_lo[1];
#if AMREX_SPACEDIM == 3
                        Real z = (static_cast<Real>(k)+0.5_rt) * dx[2] + prob_lo[2];
#endif

#if AMREX_SPACEDIM == 2
                        //std::cout << "here: " << i << " " << j << " " << k << std::endl;
                        phi_arr(i,j,k) = 0.1_rt * std::cos(2.0_rt*M_PI*y) * std::cos(2.0_rt*M_PI*x);
#else
                        phi_arr(i,j,k) = 5.0_rt * std::cos(2.0_rt*M_PI*y) * std::cos(2.0_rt*M_PI*x) * std::cos(2.0_rt*M_PI*z);
#endif

                    });

                    Array4<Real> const gphi_arr = gphi_mf.array(mfi);
                    Array4<Real> const umid_arr = umid_mf.array(mfi);


                    amrex::ParallelFor(tilebox,
                    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
                    {

                        Real x = (static_cast<Real>(i)+0.5_rt) * dx[0] + prob_lo[0];
                        Real y = (static_cast<Real>(j)+0.5_rt) * dx[1] + prob_lo[1];
#if (AMREX_SPACEDIM==3)
                        Real z = (static_cast<Real>(k)+0.5_rt) * dx[2] + prob_lo[2];
#endif

                        gphi_arr(i,j,k,0) = (phi_arr(i+1,j,k) - phi_arr(i-1,j,k)) / (2.0_rt * dx[0]);
                        gphi_arr(i,j,k,1) = (phi_arr(i,j+1,k) - phi_arr(i,j-1,k)) / (2.0_rt * dx[1]);
#if (AMREX_SPACEDIM==3)
                        gphi_arr(i,j,k,2) = (phi_arr(i,j,k+1) - phi_arr(i,j,k-1)) / (2.0_rt * dx[2]);
#endif

                        umid_arr(i,j,k,0) += gphi_arr(i,j,k,0);
                        umid_arr(i,j,k,1) += gphi_arr(i,j,k,1);
#if (AMREX_SPACEDIM==3)
                        umid_arr(i,j,k,2) += gphi_arr(i,j,k,2);
#endif

                    });

                } else {
                    amrex::Error("Not set up for these boundary conditions");
                }
            }
        }

        // fill ghosts and boundary
        AverageDown(umid, 0, AMREX_SPACEDIM);
        FillPatch(t_old, umid, umid, umid, 0, 0, AMREX_SPACEDIM, 0, bcs_u);

        // write umid
        WritePlotFile(1, t_new, dt, dummy, dummy, dummy, dummy, umid, uold,
                      uold);
        // write gphi
        WritePlotFile(2, t_new, dt, dummy, dummy, dummy, dummy, gphi, uold,
                      uold);

        // copy the velocity field over to the final state, unew
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(unew[lev], umid[lev], 0, 0, AMREX_SPACEDIM, ng_s);

            // swap pointers so NodalProj works properly
            std::swap(uold[lev], unew[lev]);
        }

    } else {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab& gphix_mac_mf = gphi_mac[lev][0];
            MultiFab& gphiy_mac_mf = gphi_mac[lev][1];
#if (AMREX_SPACEDIM == 3)
            MultiFab& gphiz_mac_mf = gphi_mac[lev][2];
#endif
            MultiFab& umac_mid_mf = umac_mid[lev][0];
            MultiFab& vmac_mid_mf = umac_mid[lev][1];
#if (AMREX_SPACEDIM == 3)
            MultiFab& wmac_mid_mf = umac_mid[lev][2];
#endif
            const Real* dx = geom[lev].CellSize();
            const auto prob_lo = geom[lev].ProbLoArray();

#ifdef _OPENMP
#pragma omp parallel
#endif
            FArrayBox phi;

            for (MFIter mfi(gphix_mac_mf, true); mfi.isValid(); ++mfi) {
                const Box& tilebox = mfi.tilebox();
                const Box& gbox = amrex::grow(tilebox, 1);

                if (   phys_bc[0] == SlipWall
                    && phys_bc[1] == SlipWall
#if AMREX_SPACEDIM == 3
                    && phys_bc[2] == SlipWall
#endif
                    && phys_bc[AMREX_SPACEDIM] == SlipWall
                    && phys_bc[AMREX_SPACEDIM+1] == SlipWall
#if AMREX_SPACEDIM == 3
                    && phys_bc[AMREX_SPACEDIM+2] == SlipWall
#endif
                    ) {


                    // Add on the gradient of a scalar (phi) that satisfies
                    // grad(phi).n = 0.

                    // x-velocity  (x are edges, y and z are centers)

                    Array4<Real> const gphix_arr = gphix_mac_mf.array(mfi);
                    Array4<Real> const umac_arr = umac_mid_mf.array(mfi);

                    amrex::ParallelFor(mfi.nodaltilebox(0),
                    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
                    {

                        Real x = static_cast<Real>(i) * dx[0] + prob_lo[0];
                        Real y = (static_cast<Real>(j)+0.5_rt) * dx[1] + prob_lo[1];
#if (AMREX_SPACEDIM==3)
                        Real z = (static_cast<Real>(k)+0.5_rt) * dx[2] + prob_lo[2];
#endif

#if AMREX_SPACEDIM == 2
                        gphix_arr(i,j,k) = 4.0_rt * x * (1.0_rt - x);
#else
                        gphix_arr(i,j,k) = 160.0_rt * x * (1.0_rt - x);
#endif
                        umac_arr(i,j,k) += gphix_arr(i,j,k);

                    });

                    // y-velocity  (x and z are centers, y are edges)

                    Array4<Real> const gphiy_arr = gphiy_mac_mf.array(mfi);
                    Array4<Real> const vmac_arr = vmac_mid_mf.array(mfi);

                    amrex::ParallelFor(mfi.nodaltilebox(1),
                    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
                    {

                        Real x = (static_cast<Real>(i)+0.5_rt) * dx[0] + prob_lo[0];
                        Real y = static_cast<Real>(j) * dx[1] + prob_lo[1];
#if (AMREX_SPACEDIM==3)
                        Real z = (static_cast<Real>(k)+0.5_rt) * dx[2] + prob_lo[2];
#endif

#if AMREX_SPACEDIM == 2
                        gphiy_arr(i,j,k) = 4.0_rt * y * (1.0_rt - y);
#else
                        gphiy_arr(i,j,k) = 160.0_rt * y * (1.0_rt - y);
#endif
                        vmac_arr(i,j,k) += gphiy_arr(i,j,k);

                    });

#if (AMREX_SPACEDIM==3)

                    // z-velocity  (x and y are centers, z are edges)

                    Array4<Real> const gphiz_arr = gphiz_mac_mf.array(mfi);
                    Array4<Real> const wmac_arr = wmac_mid_mf.array(mfi);

                    amrex::ParallelFor(mfi.nodaltilebox(2),
                    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
                    {

                        Real x = (static_cast<Real>(i)+0.5_rt) * dx[0] + prob_lo[0];
                        Real y = (static_cast<Real>(j)+0.5_rt) * dx[1] + prob_lo[1];
                        Real z = static_cast<Real>(k) * dx[2] + prob_lo[2];

                        gphiz_arr(i,j,k) = 160.0_rt * z * (1.0_rt - z);
                        wmac_arr(i,j,k) += gphiz_arr(i,j,k);

                    });
#endif

                } else if (    phys_bc[0] == Interior
                            && phys_bc[1] == Interior
#if AMREX_SPACEDIM == 3
                            && phys_bc[2] == Interior
#endif
                            && phys_bc[AMREX_SPACEDIM] == Interior
                            && phys_bc[AMREX_SPACEDIM+1] == Interior
#if AMREX_SPACEDIM == 3
                            && phys_bc[AMREX_SPACEDIM+2] == Interior
#endif
                ) {

                    // first fill phi on a local FAB on a grown tile box

                    phi.resize(gbox, 1);
                    Elixir elix_phi = phi.elixir();
                    auto phi_arr = phi.array();

                    // now we can explicitly difference phi

                    amrex::ParallelFor(gbox,
                    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
                    {

                        Real x = (static_cast<Real>(i)+0.5_rt) * dx[0] + prob_lo[0];
                        Real y = (static_cast<Real>(j)+0.5_rt) * dx[1] + prob_lo[1];
#if (AMREX_SPACEDIM==3)
                        Real z = (static_cast<Real>(k)+0.5_rt) * dx[2] + prob_lo[2];
#endif


#if AMREX_SPACEDIM == 2
                        phi_arr(i,j,k) = 0.1_rt * std::cos(2.0_rt*M_PI*y) * std::cos(2.0_rt*M_PI*x);
#else
                        phi_arr(i,j,k) = 5.0_rt * std::cos(2.0_rt*M_PI*y) * std::cos(2.0_rt*M_PI*x) * std::cos(2.0_rt*M_PI*z);
#endif

                    });

                    // lo and hi for the valid box for imposing BCs

                    const Box& tilebox = mfi.tilebox();
                    const int* lo = tilebox.loVect();
                    const int* hi = tilebox.hiVect();

                    // x-velocity  (x are edges, y and z are centers)

                    Array4<Real> const gphix_arr = gphix_mac_mf.array(mfi);
                    Array4<Real> const umac_arr = umac_mid_mf.array(mfi);

                    amrex::ParallelFor(mfi.nodaltilebox(0),
                    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
                    {

                        Real x = static_cast<Real>(i) * dx[0] + prob_lo[0];
                        Real y = (static_cast<Real>(j)+0.5_rt) * dx[1] + prob_lo[1];
#if (AMREX_SPACEDIM==3)
                        Real z = (static_cast<Real>(k)+0.5_rt) * dx[2] + prob_lo[2];
#endif

                        gphix_arr(i,j,k) = (phi_arr(i,j,k) - phi_arr(i-1,j,k)) / dx[0];
                        umac_arr(i,j,k) += gphix_arr(i,j,k);

                        // impose BCs
                        if (phys_bc[0] == SlipWall && i == lo[0]) {
                            umac_arr(i,j,k) = 0.0_rt;
                        }

                        if (phys_bc[AMREX_SPACEDIM] == SlipWall && i == hi[0]+1) {
                            umac_arr(i,j,k) = 0.0_rt;
                        }

                    });

                    // y-velocity  (x and z are centers, y are edges)

                    Array4<Real> const gphiy_arr = gphiy_mac_mf.array(mfi);
                    Array4<Real> const vmac_arr = vmac_mid_mf.array(mfi);

                    amrex::ParallelFor(mfi.nodaltilebox(1),
                    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
                    {

                        gphiy_arr(i,j,k) = (phi_arr(i,j,k) - phi_arr(i,j-1,k)) / dx[1];
                        vmac_arr(i,j,k) += gphiy_arr(i,j,k);

                        // impose BCs
                        if (phys_bc[1] == SlipWall && j == lo[1]) {
                            vmac_arr(i,j,k) = 0.0_rt;
                        }

                        if (phys_bc[AMREX_SPACEDIM+1] == SlipWall && j == hi[1]+1) {
                            vmac_arr(i,j,k) = 0.0_rt;
                        }

                    });

#if (AMREX_SPACEDIM==3)

                    // z-velocity  (x and y are centers, z are edges)

                    Array4<Real> const gphiz_arr = gphiz_mac_mf.array(mfi);
                    Array4<Real> const wmac_arr = wmac_mid_mf.array(mfi);

                    amrex::ParallelFor(mfi.nodaltilebox(2),
                    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
                    {

                        gphiz_arr(i,j,k) = (phi_arr(i,j,k) - phi_arr(i,j,k-1)) / dx[2];
                        wmac_arr(i,j,k) += gphiz_arr(i,j,k);

                        // impose BCs
                        if (phys_bc[2] == SlipWall && k == lo[2]) {
                            wmac_arr(i,j,k) = 0.0_rt;
                        }

                        if (phys_bc[AMREX_SPACEDIM+2] == SlipWall && k == hi[2]+1) {
                            wmac_arr(i,j,k) = 0.0_rt;
                        }

                    });
#endif
                } else {
                    amrex::Error("Not set up for these boundary conditions");
                }

            }
        }

        // average down and do boundaries
        if (finest_level == 0) {
            // fill periodic ghost cells
            for (int lev = 0; lev <= finest_level; ++lev) {
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    umac_mid[lev][d].FillBoundary(geom[lev].periodicity());
                }
            }
            // fill ghost cells behind physical boundaries
            FillUmacGhost(umac_mid);
        } else {
            // edge_restriction for velocities
            AverageDownFaces(umac_mid);
            // fill level n ghost cells using interpolation from level n-1 data
            FillPatchUedge(umac_mid);
        }

        // write umac_mid
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab& umac_mf = umac_mid[lev][0];
            MultiFab& vmac_mf = umac_mid[lev][1];
            MultiFab& wmac_mf = umac_mid[lev][2];
            MultiFab& utemp_mf = utemp[lev];

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(umac_mf, true); mfi.isValid(); ++mfi) {
                const Box& tilebox = mfi.tilebox();

                auto utemp_arr = utemp_mf.array(mfi);
                auto umac_arr = umac_mf.array(mfi);
                auto vmac_arr = vmac_mf.array(mfi);
#if AMREX_SPACEDIM == 3
                auto wmac_arr = wmac_mf.array(mfi);
#endif

                amrex::ParallelFor(tilebox,
                [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
                {

                    utemp_arr(i,j,k,0) = 0.5_rt * (umac_arr(i,j,k) + umac_arr(i+1,j,k));
                    utemp_arr(i,j,k,1) = 0.5_rt * (vmac_arr(i,j,k) + vmac_arr(i,j+1,k));
#if (AMREX_SPACEDIM==3)
                    utemp_arr(i,j,k,2) = 0.5_rt * (wmac_arr(i,j,k) + wmac_arr(i,j,k+1));
#endif

                });

            }
        }

        WritePlotFile(1, t_new, dt, dummy, dummy, dummy, dummy, utemp, uold,
                      uold);

        // write gphi_mac
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab& umac_mf = gphi_mac[lev][0];
            MultiFab& vmac_mf = gphi_mac[lev][1];
            MultiFab& wmac_mf = gphi_mac[lev][2];
            MultiFab& utemp_mf = utemp[lev];

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(umac_mf, true); mfi.isValid(); ++mfi) {
                const Box& tilebox = mfi.tilebox();

                auto utemp_arr = utemp_mf.array(mfi);
                auto umac_arr = umac_mf.array(mfi);
                auto vmac_arr = vmac_mf.array(mfi);
#if AMREX_SPACEDIM == 3
                auto wmac_arr = wmac_mf.array(mfi);
#endif

                amrex::ParallelFor(tilebox,
                [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
                {

                    utemp_arr(i,j,k,0) = 0.5_rt * (umac_arr(i,j,k) + umac_arr(i+1,j,k));
                    utemp_arr(i,j,k,1) = 0.5_rt * (vmac_arr(i,j,k) + vmac_arr(i,j+1,k));
#if (AMREX_SPACEDIM==3)
                    utemp_arr(i,j,k,2) = 0.5_rt * (wmac_arr(i,j,k) + wmac_arr(i,j,k+1));
#endif

                });
            }
        }

        WritePlotFile(2, t_new, dt, dummy, dummy, dummy, dummy, utemp, uold,
                      uold);

        // copy the velocity field over to the final state, unew
        for (int lev = 0; lev <= finest_level; ++lev) {
            for (auto comp = 0; comp < AMREX_SPACEDIM; ++comp)
                MultiFab::Copy(umac_new[lev][comp], umac_mid[lev][comp], 0, 0,
                               1, umac_new[lev][comp].nGrow());
        }
    }

    //--------------------------------------------------------------------------
    // project out the divergent portion of the velocity field
    //--------------------------------------------------------------------------

    Print() << "...projection" << std::endl;

    if (project_type == 1) {
        // hgprojection -- here pi is nodal and u is cell-centered

        for (int lev = 0; lev <= finest_level; ++lev) {
            // build the density used in the projection -- we are just doing
            // constant density, so set it to 1
            sold[lev].setVal(0.);
            snew[lev].setVal(0.);
            sold[lev].setVal(1., Rho, 1, 1);
            snew[lev].setVal(1., Rho, 1, 1);
            pi[lev].setVal(0.);
            gpi[lev].setVal(0.);
            rhcc_for_nodalproj[lev].setVal(0.);
        }

        // build the coefficient in the divergence.  We are doing
        // divergence-free (incompressible), so set beta0 = 1
        beta0_old.setVal(1.0);
        beta0_new.setVal(1.0);

        t_new = t_old + 1.;

        // NodalProj is going to operate on unew, where we've temporarily stored the initial data. Let's instead copy this initial data to umid.
        for (int lev = 0; lev <= finest_level; ++lev)
            std::swap(umid[lev], unew[lev]);

        // hgproject
        NodalProj(initial_projection_comp, rhcc_for_nodalproj);

        // swap pointers to restore initial data to uold and new data to unew
        for (int lev = 0; lev <= finest_level; ++lev) {
            std::swap(uold[lev], umid[lev]);
            std::swap(umid[lev], unew[lev]);
        }

        // write unew
        WritePlotFile(3, t_new, dt, dummy, dummy, dummy, dummy, unew, uold,
                      uold);

        // I think now can compare to uold and see if it's the same?
        {
            int lev = finest_level;
            Real norm = 0.;
            umid[lev].setVal(0.);

            MultiFab::Copy(umid[lev], uold[lev], 0, 0, AMREX_SPACEDIM, 0);
            MultiFab::Subtract(umid[lev], unew[lev], 0, 0, AMREX_SPACEDIM, 0);

            for (auto comp = 0; comp < AMREX_SPACEDIM; ++comp)
                norm += umid[lev].norm2(comp) / uold[lev].norm2(comp);

            Print() << "\nRelative error = " << norm << std::endl;

            norm = 0.;
            for (auto comp = 0; comp < AMREX_SPACEDIM; ++comp) {
                norm = umid[lev].norm0(comp);

                Print() << "\nAbsolute error, dim " << comp << " = " << norm
                        << std::endl;
            }
        }

    } else {
        // mac projection -- here pi is cell-centered and u is MAC

        Vector<MultiFab> macpi(finest_level + 1);
        Vector<MultiFab> macrhs(finest_level + 1);

        beta0_old.setVal(1.0);
        beta0_new.setVal(1.0);

        for (int lev = 0; lev <= finest_level; ++lev) {
            // cell-centered MultiFabs
            macpi[lev].define(grids[lev], dmap[lev], 1, 1);
            macpi[lev].setVal(0.);
            macrhs[lev].define(grids[lev], dmap[lev], 1, 1);
            macrhs[lev].setVal(0.);

            sold[lev].setVal(0.);
            snew[lev].setVal(0.);
            sold[lev].setVal(1., Rho, 1, 1);
            snew[lev].setVal(1., Rho, 1, 1);
        }

        // macproject
        auto is_predictor = 0;
        MacProj(umac_new, macpi, macrhs, beta0_old, is_predictor);

        // I think now can compare to umac_old and see if it's the same?
        {
            int lev = finest_level;
            Real norm = 0.;
            for (auto comp = 0; comp < AMREX_SPACEDIM; ++comp) {
                umac_mid[lev][comp].setVal(0.);
                MultiFab::Copy(umac_mid[lev][comp], umac_old[lev][comp], 0, 0,
                               1, 0);
                MultiFab::Subtract(umac_mid[lev][comp], umac_new[lev][comp], 0,
                                   0, 1, 0);

                norm +=
                    umac_mid[lev][comp].norm2() / umac_old[lev][comp].norm2();
            }
            Print() << "\nRelative error = " << norm << std::endl;
        }

        // write umac_new
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab& umac_mf = umac_new[lev][0];
            MultiFab& vmac_mf = umac_new[lev][1];
            MultiFab& wmac_mf = umac_new[lev][2];
            MultiFab& utemp_mf = utemp[lev];

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(umac_mf, true); mfi.isValid(); ++mfi) {
                const Box& tilebox = mfi.tilebox();

                auto utemp_arr = utemp_mf.array(mfi);
                auto umac_arr = umac_mf.array(mfi);
                auto vmac_arr = vmac_mf.array(mfi);
#if AMREX_SPACEDIM == 3
                auto wmac_arr = wmac_mf.array(mfi);
#endif

                amrex::ParallelFor(tilebox,
                [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
                {

                    utemp_arr(i,j,k,0) = 0.5_rt * (umac_arr(i,j,k) + umac_arr(i+1,j,k));
                    utemp_arr(i,j,k,1) = 0.5_rt * (vmac_arr(i,j,k) + vmac_arr(i,j+1,k));
#if (AMREX_SPACEDIM==3)
                    utemp_arr(i,j,k,2) = 0.5_rt * (wmac_arr(i,j,k) + wmac_arr(i,j,k+1));
#endif

                });

            }
        }

        WritePlotFile(3, t_new, dt, dummy, dummy, dummy, dummy, utemp, uold,
                      uold);
    }
}


#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

// compute eta_rho at edge- and cell-centers
void Maestro::MakeEtarho(const Vector<MultiFab>& etarho_flux) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeEtarho()", MakeEtarho);

#ifdef AMREX_USE_CUDA
    bool launched;
    if (deterministic_nodal_solve) {
        launched = !Gpu::notInLaunchRegion();
        // turn off GPU
        if (launched) Gpu::setLaunchRegion(false);
    }
#endif

    // Local variables
    BaseState<Real> etarhosum_s(base_geom.max_radial_level + 1,
                                base_geom.nr_fine + 1);
    etarhosum_s.setVal(0.0);
    auto etarhosum = etarhosum_s.array();

    // this stores how many cells there are laterally at each level
    BaseState<int> ncell_s(base_geom.max_radial_level + 1);
    auto ncell = ncell_s.array();

    for (int lev = 0; lev <= finest_level; ++lev) {
        // Get the index space of the domain
        const Box& domainBox = geom[lev].Domain();

        // compute number of cells at any given height for each level
        if (AMREX_SPACEDIM == 2) {
            ncell(lev) = domainBox.bigEnd(0) + 1;
        } else if (AMREX_SPACEDIM == 3) {
            ncell(lev) = (domainBox.bigEnd(0) + 1) * (domainBox.bigEnd(1) + 1);
        }

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction)
#endif
        for (MFIter mfi(sold[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid tile region
            const Box& tilebox = mfi.tilebox();

            const Array4<const Real> etarhoflux_arr =
                etarho_flux[lev].array(mfi);

#if (AMREX_SPACEDIM == 2)
            int zlo = tilebox.loVect3d()[2];
            ParallelFor(tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                if (k == zlo) {
                    amrex::HostDevice::Atomic::Add(&(etarhosum(lev, j)),
                                                   etarhoflux_arr(i, j, k));
                }
            });

            // we only add the contribution at the top edge if we are at the top of the domain
            // this prevents double counting
            auto top_edge = false;
            for (auto i = 1; i <= base_geom.numdisjointchunks(lev); ++i) {
                if (tilebox.hiVect3d()[1] == base_geom.r_end_coord(lev, i)) {
                    top_edge = true;
                }
            }

            if (top_edge) {
                const int k = 0;
                const auto ybx = mfi.nodaltilebox(1);
                const int j = ybx.hiVect3d()[1];
                const int lo = ybx.loVect3d()[0];
                const int hi = ybx.hiVect3d()[0];

                ParallelFor(hi - lo + 1, [=] AMREX_GPU_DEVICE(int n) {
                    int i = n + lo;
                    amrex::HostDevice::Atomic::Add(&(etarhosum(lev, j)),
                                                   etarhoflux_arr(i, j, k));
                });
                Gpu::synchronize();
            }
#else
            ParallelFor(tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                amrex::HostDevice::Atomic::Add(&(etarhosum(lev, k)),
                                               etarhoflux_arr(i, j, k));
            });

            // we only add the contribution at the top edge if we are at the top of the domain
            // this prevents double counting
            auto top_edge = false;

            for (auto i = 1; i <= base_geom.numdisjointchunks(lev); ++i) {
                if (tilebox.hiVect3d()[2] == base_geom.r_end_coord(lev, i)) {
                    top_edge = true;
                }
            }

            if (top_edge) {
                const auto zbx = mfi.nodaltilebox(2);
                int zhi = zbx.hiVect3d()[2];
                ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    if (k == zhi) {
                        amrex::HostDevice::Atomic::Add(&(etarhosum(lev, k)),
                                                       etarhoflux_arr(i, j, k));
                    }
                });
            }
#endif
        }
    }

    ParallelDescriptor::ReduceRealSum(
        etarhosum.dataPtr(),
        (base_geom.nr_fine + 1) * (base_geom.max_radial_level + 1));

    etarho_ec.setVal(0.0);
    etarho_cc.setVal(0.0);

    auto etarho_ec_arr = etarho_ec.array();
    auto etarho_cc_arr = etarho_cc.array();
    const auto etarhosum_arr = etarhosum_s.const_array();

    for (auto n = 0; n <= base_geom.finest_radial_level; ++n) {
        for (auto i = 1; i <= base_geom.numdisjointchunks(n); ++i) {
            const int lo = base_geom.r_start_coord(n, i);
            const int hi = base_geom.r_end_coord(n, i) + 1;
            const auto ncell_lev = ncell(n);
            ParallelFor(hi - lo + 1, [=] AMREX_GPU_DEVICE(int j) {
                int r = j + lo;
                etarho_ec_arr(n, r) = etarhosum_arr(n, r) / Real(ncell_lev);
            });
            Gpu::synchronize();
        }
    }

    // These calls shouldn't be needed since the planar algorithm doesn't use
    // these outside of this function, but this is just to be safe in case
    // things change in the future.
    RestrictBase(etarho_ec, false);
    FillGhostBase(etarho_ec, false);

    // make the cell-centered etarho_cc by averaging etarho to centers
    for (auto n = 0; n <= base_geom.finest_radial_level; ++n) {
        for (auto i = 1; i <= base_geom.numdisjointchunks(n); ++i) {
            const int lo = base_geom.r_start_coord(n, i);
            const int hi = base_geom.r_end_coord(n, i);
            ParallelFor(hi - lo + 1, [=] AMREX_GPU_DEVICE(int j) {
                int r = j + lo;
                etarho_cc_arr(n, r) =
                    0.5 * (etarho_ec_arr(n, r) + etarho_ec_arr(n, r + 1));
            });
            Gpu::synchronize();
        }
    }

    // These calls shouldn't be needed since the planar algorithm only uses
    // etarho_cc to make_psi, and then we fill ghost cells in make_psi, but
    // this is just to be safe in case things change in the future
    RestrictBase(etarho_cc, true);
    FillGhostBase(etarho_cc, true);

#ifdef AMREX_USE_CUDA
    if (deterministic_nodal_solve) {
        // turn GPU back on
        if (launched) Gpu::setLaunchRegion(true);
    }
#endif
}

#if AMREX_SPACEDIM == 3
void Maestro::MakeEtarhoSphr(
    const Vector<MultiFab>& scal_old,
    const Vector<MultiFab>& scal_new,
    const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& umac,
    const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& w0mac) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeEtarhoSphr()", MakeEtarhoSphr);

    const int max_lev = base_geom.max_radial_level + 1;
    const int nrf = base_geom.nr_fine + 1;

    Vector<MultiFab> eta_cart(finest_level + 1);
    Vector<MultiFab> rho0_nph_cart(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        eta_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        rho0_nph_cart[lev].define(grids[lev], dmap[lev], 1, 0);
    }

    BaseState<Real> rho0_nph(max_lev, base_geom.nr_fine);
    rho0_nph.copy(0.5 * (rho0_old + rho0_new));

    Put1dArrayOnCart(rho0_nph, rho0_nph_cart, false, false, bcs_f, 0);

    for (int lev = 0; lev <= finest_level; ++lev) {
        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(scal_old[lev], TilingIfNotGPU()); mfi.isValid();
             ++mfi) {
            // Get the index space of the valid region
            const Box& tilebox = mfi.tilebox();

            const Array4<const Real> rho_old = scal_old[lev].array(mfi);
            const Array4<const Real> rho_new = scal_new[lev].array(mfi);
            const Array4<const Real> rho0_nph_cart_arr =
                rho0_nph_cart[lev].array(mfi);
            const Array4<const Real> umac_arr = umac[lev][0].array(mfi);
            const Array4<const Real> vmac = umac[lev][1].array(mfi);
            const Array4<const Real> wmac = umac[lev][2].array(mfi);
            const Array4<const Real> w0macx = w0mac[lev][0].array(mfi);
            const Array4<const Real> w0macy = w0mac[lev][1].array(mfi);
            const Array4<const Real> w0macz = w0mac[lev][2].array(mfi);
            const Array4<const Real> normal_arr = normal[lev].array(mfi);
            const Array4<Real> eta_cart_arr = eta_cart[lev].array(mfi);

            ParallelFor(tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Real U_dot_er = 0.5 *
                                    (umac_arr(i, j, k) + umac_arr(i + 1, j, k) +
                                     w0macx(i, j, k) + w0macx(i + 1, j, k)) *
                                    normal_arr(i, j, k, 0) +
                                0.5 *
                                    (vmac(i, j, k) + vmac(i, j + 1, k) +
                                     w0macy(i, j, k) + w0macy(i, j + 1, k)) *
                                    normal_arr(i, j, k, 1) +
                                0.5 *
                                    (wmac(i, j, k) + wmac(i, j, k + 1) +
                                     w0macz(i, j, k) + w0macz(i, j, k + 1)) *
                                    normal_arr(i, j, k, 2);

                // construct time-centered [ rho' (U dot e_r) ]
                eta_cart_arr(i, j, k) =
                    (0.5 * (rho_old(i, j, k) + rho_new(i, j, k)) -
                     rho0_nph_cart_arr(i, j, k)) *
                    U_dot_er;
            });
        }  // end MFIter loop
    }      // end loop over levels


    // average fine data onto coarser cells & fill ghost cells
    AverageDown(eta_cart, 0, 1);
    FillPatch(t_old, eta_cart, eta_cart, eta_cart, 0, 0, 1, 0, bcs_f);

    // compute etarho_cc as the average of eta_cart = [ rho' (U dot e_r) ]
    Average(eta_cart, etarho_cc, 0);

    const auto& r_cc_loc = base_geom.r_cc_loc;
    const auto& r_edge_loc = base_geom.r_edge_loc;
    auto etarho_ec_arr = etarho_ec.array();
    const auto etarho_cc_arr = etarho_cc.const_array();

    // put eta on base state edges
    // note that in spherical the base state has no refinement
    // the 0th value of etarho = 0, since U dot . e_r must be
    // zero at the center (since e_r is not defined there)
    ParallelFor(nrf, [=] AMREX_GPU_DEVICE(int r) {
        if (r == 0) {
            etarho_ec_arr(0, r) = 0.0;
        } else if (r == nrf - 1) {
            // probably should do some better extrapolation here eventually
            etarho_ec_arr(0, r) = etarho_cc_arr(0, r - 1);
        } else {
            if (spherical) {
                Real dr1 = r_cc_loc(0, r) - r_edge_loc(0, r);
                Real dr2 = r_edge_loc(0, r) - r_cc_loc(0, r - 1);
                etarho_ec_arr(0, r) = (dr2 * etarho_cc_arr(0, r) +
                                       dr1 * etarho_cc_arr(0, r - 1)) /
                                      (dr1 + dr2);
            } else {
                etarho_ec_arr(0, r) =
                    0.5 * (etarho_cc_arr(0, r) + etarho_cc_arr(0, r - 1));
            }
        }
    });
    Gpu::synchronize();
}
#endif

void Maestro::MakeEtarhoPlanar(
    const Vector<MultiFab>& scal_old, const Vector<MultiFab>& scal_new,
    const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& umac) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeEtarhoPlanar()", MakeEtarhoPlanar);

    const int max_lev = base_geom.max_radial_level + 1;

    Vector<MultiFab> eta_cart(finest_level + 1);
    Vector<MultiFab> rho0_nph_cart(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        eta_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        rho0_nph_cart[lev].define(grids[lev], dmap[lev], 1, 0);
    }

    BaseState<Real> rho0_nph(max_lev, base_geom.nr_fine);
    rho0_nph.copy(0.5 * (rho0_old + rho0_new));

    Put1dArrayOnCart(rho0_nph, rho0_nph_cart, false, false, bcs_f, 0);

    for (int lev = 0; lev <= finest_level; ++lev) {
        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(scal_old[lev], TilingIfNotGPU()); mfi.isValid();
             ++mfi) {
            // Get the index space of the valid region
            const Box& tilebox = mfi.tilebox();

            const Array4<const Real> rho_old = scal_old[lev].array(mfi);
            const Array4<const Real> rho_new = scal_new[lev].array(mfi);
            const Array4<const Real> rho0_nph_cart_arr =
                rho0_nph_cart[lev].array(mfi);
#if AMREX_SPACEDIM == 2
            const Array4<const Real> vmac = umac[lev][1].array(mfi);
#else
            const Array4<const Real> vmac = umac[lev][2].array(mfi);
#endif
            const Array4<Real> eta_cart_arr = eta_cart[lev].array(mfi);

            ParallelFor(tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

#if AMREX_SPACEDIM == 2
                Real U_dot_er = 0.5 * (vmac(i, j, k) + vmac(i, j + 1, k));
#else
                Real U_dot_er = 0.5 *(vmac(i, j, k) + vmac(i, j, k + 1));
#endif
                // construct time-centered [ rho' (U dot e_r) ]
                eta_cart_arr(i, j, k) =
                    (0.5 * (rho_old(i, j, k) + rho_new(i, j, k)) -
                     rho0_nph_cart_arr(i, j, k)) *
                    U_dot_er;
            });
        }  // end MFIter loop
    }      // end loop over levels

    // average fine data onto coarser cells & fill ghost cells
    AverageDown(eta_cart, 0, 1);
    FillPatch(t_old, eta_cart, eta_cart, eta_cart, 0, 0, 1, 0, bcs_f);

    // compute etarho_cc as the average of eta_cart = [ rho' (U dot e_r) ]
    Average(eta_cart, etarho_cc, 0);
}

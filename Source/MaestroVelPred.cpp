#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void Maestro::VelPred(
    Vector<MultiFab>& utilde, const Vector<MultiFab>& ufull,
    const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& utrans,
    Vector<std::array<MultiFab, AMREX_SPACEDIM> >& umac,
    const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& w0mac,
    const Vector<MultiFab>& force) {

    amrex::ignore_unused(w0mac);

    // timer for profiling
    BL_PROFILE_VAR("Maestro::VelPred()", VelPred);

    for (int lev = 0; lev <= finest_level; ++lev) {
        // Get the index space and grid spacing of the domain
        const Box& domainBox = geom[lev].Domain();
        const auto dx = geom[lev].CellSizeArray();

        // get references to the MultiFabs at level lev
        MultiFab& utilde_mf = utilde[lev];
        const MultiFab& ufull_mf = ufull[lev];
        MultiFab& umac_mf = umac[lev][0];
        const MultiFab& utrans_mf = utrans[lev][0];
        const MultiFab& vtrans_mf = utrans[lev][1];
        MultiFab& vmac_mf = umac[lev][1];

        MultiFab Ipu, Imu, Ipv, Imv;
        Ipu.define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        Imu.define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        Ipv.define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        Imv.define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);

        MultiFab Ipfx, Imfx, Ipfy, Imfy;
        Ipfx.define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        Imfx.define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        Ipfy.define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        Imfy.define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);

        MultiFab ulx, urx, uimhx, uly, ury, uimhy;
        ulx.define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        urx.define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        uimhx.define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        uly.define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        ury.define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        uimhy.define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
#if (AMREX_SPACEDIM == 3)
        const MultiFab& wtrans_mf = utrans[lev][2];
        MultiFab& wmac_mf = umac[lev][2];
        const MultiFab& w0macx_mf = w0mac[lev][0];
        const MultiFab& w0macy_mf = w0mac[lev][1];
        const MultiFab& w0macz_mf = w0mac[lev][2];

        MultiFab Ipw, Imw;
        Ipw.define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        Imw.define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);

        MultiFab Ipfz, Imfz;
        Ipfz.define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        Imfz.define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);

        MultiFab ulz, urz, uimhz;
        ulz.define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        urz.define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        uimhz.define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);

        MultiFab uimhyz, uimhzy, vimhxz, vimhzx, wimhxy, wimhyx;
        uimhyz.define(grids[lev], dmap[lev], 1, 1);
        uimhzy.define(grids[lev], dmap[lev], 1, 1);
        vimhxz.define(grids[lev], dmap[lev], 1, 1);
        vimhzx.define(grids[lev], dmap[lev], 1, 1);
        wimhxy.define(grids[lev], dmap[lev], 1, 1);
        wimhyx.define(grids[lev], dmap[lev], 1, 1);
#endif
        const MultiFab& force_mf = force[lev];
        const MultiFab& w0_mf = w0_cart[lev];

#if (AMREX_SPACEDIM == 2)

        MultiFab u_mf, v_mf;

        if (ppm_type == 0) {
            u_mf.define(grids[lev], dmap[lev], 1, utilde[lev].nGrow());
            v_mf.define(grids[lev], dmap[lev], 1, utilde[lev].nGrow());

            MultiFab::Copy(u_mf, utilde[lev], 0, 0, 1, utilde[lev].nGrow());
            MultiFab::Copy(v_mf, utilde[lev], 1, 0, 1, utilde[lev].nGrow());

        } else if (ppm_type == 1 || ppm_type == 2) {
            u_mf.define(grids[lev], dmap[lev], 1, ufull[lev].nGrow());
            v_mf.define(grids[lev], dmap[lev], 1, ufull[lev].nGrow());

            MultiFab::Copy(u_mf, ufull[lev], 0, 0, 1, ufull[lev].nGrow());
            MultiFab::Copy(v_mf, ufull[lev], 1, 0, 1, ufull[lev].nGrow());
        }

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(utilde_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region
            const Box& obx = amrex::grow(mfi.tilebox(), 1);

            if (ppm_type == 0) {
                // we're going to reuse Ip here as slopex as it has the
                // correct number of ghost zones
                Slopex(obx, utilde_mf.array(mfi), Ipu.array(mfi), domainBox,
                       bcs_u, AMREX_SPACEDIM, 0);
            } else {
                PPM(obx, utilde_mf.array(mfi), u_mf.array(mfi), v_mf.array(mfi),
                    Ipu.array(mfi), Imu.array(mfi), domainBox, bcs_u, dx, false,
                    0, 0);

                if (ppm_trace_forces == 1) {
                    PPM(obx, force_mf.array(mfi), u_mf.array(mfi),
                        v_mf.array(mfi), Ipfx.array(mfi), Imfx.array(mfi),
                        domainBox, bcs_u, dx, false, 0, 0);
                }
            }

            if (ppm_type == 0) {
                // we're going to reuse Im here as slopey as it has the
                // correct number of ghost zones
                Slopey(obx, utilde_mf.array(mfi), Imv.array(mfi), domainBox,
                       bcs_u, AMREX_SPACEDIM, 0);
            } else {
                PPM(obx, utilde_mf.array(mfi), u_mf.array(mfi), v_mf.array(mfi),
                    Ipv.array(mfi), Imv.array(mfi), domainBox, bcs_u, dx, false,
                    1, 1);

                if (ppm_trace_forces == 1) {
                    PPM(obx, force_mf.array(mfi), u_mf.array(mfi),
                        v_mf.array(mfi), Ipv.array(mfi), Imv.array(mfi),
                        domainBox, bcs_u, dx, false, 1, 1);
                }
            }

            Gpu::synchronize();

            VelPredInterface(mfi, utilde_mf.array(mfi), ufull_mf.array(mfi),
                             utrans_mf.array(mfi), vtrans_mf.array(mfi),
                             Imu.array(mfi), Ipu.array(mfi), Imv.array(mfi),
                             Ipv.array(mfi), ulx.array(mfi), urx.array(mfi),
                             uimhx.array(mfi), uly.array(mfi), ury.array(mfi),
                             uimhy.array(mfi), domainBox, dx);

            Gpu::synchronize();

            VelPredVelocities(mfi, utilde_mf.array(mfi), utrans_mf.array(mfi),
                              vtrans_mf.array(mfi), umac_mf.array(mfi),
                              vmac_mf.array(mfi), Imfx.array(mfi),
                              Ipfx.array(mfi), Imfy.array(mfi), Ipfy.array(mfi),
                              ulx.array(mfi), urx.array(mfi), uimhx.array(mfi),
                              uly.array(mfi), ury.array(mfi), uimhy.array(mfi),
                              force_mf.array(mfi), w0_mf.array(mfi), domainBox,
                              dx);
        }  // end MFIter loop

#elif (AMREX_SPACEDIM == 3)

        MultiFab u_mf, v_mf, w_mf;

        if (ppm_type == 0) {
            u_mf.define(grids[lev], dmap[lev], 1, utilde[lev].nGrow());
            v_mf.define(grids[lev], dmap[lev], 1, utilde[lev].nGrow());
            w_mf.define(grids[lev], dmap[lev], 1, utilde[lev].nGrow());

            MultiFab::Copy(u_mf, utilde[lev], 0, 0, 1, utilde[lev].nGrow());
            MultiFab::Copy(v_mf, utilde[lev], 1, 0, 1, utilde[lev].nGrow());
            MultiFab::Copy(w_mf, utilde[lev], 2, 0, 1, utilde[lev].nGrow());

        } else if (ppm_type == 1 || ppm_type == 2) {
            u_mf.define(grids[lev], dmap[lev], 1, ufull[lev].nGrow());
            v_mf.define(grids[lev], dmap[lev], 1, ufull[lev].nGrow());
            w_mf.define(grids[lev], dmap[lev], 1, ufull[lev].nGrow());

            MultiFab::Copy(u_mf, ufull[lev], 0, 0, 1, ufull[lev].nGrow());
            MultiFab::Copy(v_mf, ufull[lev], 1, 0, 1, ufull[lev].nGrow());
            MultiFab::Copy(w_mf, ufull[lev], 2, 0, 1, ufull[lev].nGrow());
        }
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(utilde_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region
            const Box& obx = amrex::grow(mfi.tilebox(), 1);

            // x-direction
            if (ppm_type == 0) {
                // we're going to reuse Ipu here as slopex as it has the
                // correct number of ghost zones
                Slopex(obx, utilde_mf.array(mfi), Ipu.array(mfi), domainBox,
                       bcs_u, AMREX_SPACEDIM, 0);

            } else {
                PPM(obx, utilde_mf.array(mfi), u_mf.array(mfi), v_mf.array(mfi),
                    w_mf.array(mfi), Ipu.array(mfi), Imu.array(mfi), domainBox,
                    bcs_u, dx, false, 0, 0);

                if (ppm_trace_forces == 1) {
                    PPM(obx, force_mf.array(mfi), u_mf.array(mfi),
                        v_mf.array(mfi), w_mf.array(mfi), Ipfx.array(mfi),
                        Imfx.array(mfi), domainBox, bcs_u, dx, false, 0, 0);
                }
            }

            // y-direction
            if (ppm_type == 0) {
                // we're going to reuse Imv here as slopey as it has the
                // correct number of ghost zones
                Slopey(obx, utilde_mf.array(mfi), Imv.array(mfi), domainBox,
                       bcs_u, AMREX_SPACEDIM, 0);

            } else {
                PPM(obx, utilde_mf.array(mfi), u_mf.array(mfi), v_mf.array(mfi),
                    w_mf.array(mfi), Ipv.array(mfi), Imv.array(mfi), domainBox,
                    bcs_u, dx, false, 1, 1);

                if (ppm_trace_forces == 1) {
                    PPM(obx, force_mf.array(mfi), u_mf.array(mfi),
                        v_mf.array(mfi), w_mf.array(mfi), Ipfy.array(mfi),
                        Imfy.array(mfi), domainBox, bcs_u, dx, false, 1, 1);
                }
            }

            // z-direction
            if (ppm_type == 0) {
                // we're going to reuse Imw here as slopey as it has the
                // correct number of ghost zones

                Slopez(obx, utilde_mf.array(mfi), Imw.array(mfi), domainBox,
                       bcs_u, AMREX_SPACEDIM, 0);

            } else {
                PPM(obx, utilde_mf.array(mfi), u_mf.array(mfi), v_mf.array(mfi),
                    w_mf.array(mfi), Ipw.array(mfi), Imw.array(mfi), domainBox,
                    bcs_u, dx, false, 2, 2);

                if (ppm_trace_forces == 1) {
                    PPM(obx, force_mf.array(mfi), u_mf.array(mfi),
                        v_mf.array(mfi), w_mf.array(mfi), Ipfz.array(mfi),
                        Imfz.array(mfi), domainBox, bcs_u, dx, false, 2, 2);
                }
            }

            Gpu::synchronize();

            VelPredInterface(mfi, utilde_mf.array(mfi), ufull_mf.array(mfi),
                             utrans_mf.array(mfi), vtrans_mf.array(mfi),
                             wtrans_mf.array(mfi), Imu.array(mfi),
                             Ipu.array(mfi), Imv.array(mfi), Ipv.array(mfi),
                             Imw.array(mfi), Ipw.array(mfi), ulx.array(mfi),
                             urx.array(mfi), uimhx.array(mfi), uly.array(mfi),
                             ury.array(mfi), uimhy.array(mfi), ulz.array(mfi),
                             urz.array(mfi), uimhz.array(mfi), domainBox, dx);

            Gpu::synchronize();

            VelPredTransverse(
                mfi, utilde_mf.array(mfi), utrans_mf.array(mfi),
                vtrans_mf.array(mfi), wtrans_mf.array(mfi), ulx.array(mfi),
                urx.array(mfi), uimhx.array(mfi), uly.array(mfi),
                ury.array(mfi), uimhy.array(mfi), ulz.array(mfi),
                urz.array(mfi), uimhz.array(mfi), uimhyz.array(mfi),
                uimhzy.array(mfi), vimhxz.array(mfi), vimhzx.array(mfi),
                wimhxy.array(mfi), wimhyx.array(mfi), domainBox, dx);

            Gpu::synchronize();

            VelPredVelocities(
                mfi, utilde_mf.array(mfi), utrans_mf.array(mfi),
                vtrans_mf.array(mfi), wtrans_mf.array(mfi), umac_mf.array(mfi),
                vmac_mf.array(mfi), wmac_mf.array(mfi), w0macx_mf.array(mfi),
                w0macy_mf.array(mfi), w0macz_mf.array(mfi), Imfx.array(mfi),
                Ipfx.array(mfi), Imfy.array(mfi), Ipfy.array(mfi),
                Imfz.array(mfi), Ipfz.array(mfi), ulx.array(mfi),
                urx.array(mfi), uly.array(mfi), ury.array(mfi), ulz.array(mfi),
                urz.array(mfi), uimhyz.array(mfi), uimhzy.array(mfi),
                vimhxz.array(mfi), vimhzx.array(mfi), wimhxy.array(mfi),
                wimhyx.array(mfi), force_mf.array(mfi), w0_mf.array(mfi),
                domainBox, dx);
        }  // end MFIter loop

#endif  // AMREX_SPACEDIM
    }   // end loop over levels

    // edge_restriction
    AverageDownFaces(umac);
}

#if (AMREX_SPACEDIM == 2)

void Maestro::VelPredInterface(
    const MFIter& mfi, Array4<const Real> const utilde,
    Array4<const Real> const ufull, Array4<const Real> const utrans,
    Array4<const Real> const vtrans, Array4<const Real> const Imu,
    Array4<const Real> const Ipu, Array4<const Real> const Imv,
    Array4<const Real> const Ipv, Array4<Real> const ulx,
    Array4<Real> const urx, Array4<Real> const uimhx, Array4<Real> const uly,
    Array4<Real> const ury, Array4<Real> const uimhy, const Box& domainBox,
    const amrex::GpuArray<Real, AMREX_SPACEDIM> dx) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::VelPredInterface()", VelPredInterface);

    // NOTE: for ppm_type == 0, slopex == Ipu, slopey == Imv

    ////////////////////////////////////
    // Create u_{\i-\half\e_x}^x, etc.
    ////////////////////////////////////

    const Box& obx = amrex::grow(mfi.tilebox(), 1);
    const Box& mxbx = amrex::growLo(obx, 0, -1);
    const Box& mybx = amrex::growLo(obx, 1, -1);

    const Real dt2 = 0.5 * dt;

    const Real hx = dx[0];
    const Real hy = dx[1];

    const auto rel_eps_local = rel_eps;

    // x-direction
    const auto domlo = domainBox.loVect3d();
    const auto domhi = domainBox.hiVect3d();

    int bclo = phys_bc[0];
    int bchi = phys_bc[AMREX_SPACEDIM];

    ParallelFor(mxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        if (ppm_type == 0) {
            Real maxu = amrex::max(0.0, ufull(i - 1, j, k, 0));
            Real minu = amrex::min(0.0, ufull(i, j, k, 0));
            // extrapolate both components of velocity to left face
            ulx(i, j, k, 0) = utilde(i - 1, j, k, 0) +
                              (0.5 - (dt2 / hx) * maxu) * Ipu(i - 1, j, k, 0);
            ulx(i, j, k, 1) = utilde(i - 1, j, k, 1) +
                              (0.5 - (dt2 / hx) * maxu) * Ipu(i - 1, j, k, 1);
            // extrapolate both components of velocity to right face
            urx(i, j, k, 0) = utilde(i, j, k, 0) -
                              (0.5 + (dt2 / hx) * minu) * Ipu(i, j, k, 0);
            urx(i, j, k, 1) = utilde(i, j, k, 1) -
                              (0.5 + (dt2 / hx) * minu) * Ipu(i, j, k, 1);
        } else if (ppm_type == 1 || ppm_type == 2) {
            // extrapolate both components of velocity to left face
            ulx(i, j, k, 0) = Ipu(i - 1, j, k, 0);
            ulx(i, j, k, 1) = Ipv(i - 1, j, k, 0);
            // extrapolate both components of velocity to right face
            urx(i, j, k, 0) = Imu(i, j, k, 0);
            urx(i, j, k, 1) = Imv(i, j, k, 0);
        }

        // impose lo side bc's
        if (i == domlo[0]) {
            switch (bclo) {
                case Inflow:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        ulx(i, j, k, n) = utilde(i - 1, j, k, n);
                        urx(i, j, k, n) = utilde(i - 1, j, k, n);
                    }
                    break;
                case SlipWall:
                case Symmetry:
                    ulx(i, j, k, 0) = 0.0;
                    urx(i, j, k, 0) = 0.0;
                    ulx(i, j, k, 1) = urx(i, j, k, 1);
                    break;
                case NoSlipWall:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        ulx(i, j, k, n) = 0.0;
                        urx(i, j, k, n) = 0.0;
                    }
                    break;
                case Outflow:
                    urx(i, j, k, 0) = amrex::min(urx(i, j, k, 0), 0.0);
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        urx(i, j, k, n) = ulx(i, j, k, n);
                    }
                    break;
                case Interior:
                    break;
            }

            // impose hi side bc's
        } else if (i == domhi[0] + 1) {
            switch (bchi) {
                case Inflow:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        ulx(i, j, k, n) = utilde(i, j, k, n);
                        urx(i, j, k, n) = utilde(i, j, k, n);
                    }
                    break;
                case SlipWall:
                case Symmetry:
                    ulx(i, j, k, 0) = 0.0;
                    urx(i, j, k, 0) = 0.0;
                    urx(i, j, k, 1) = ulx(i, j, k, 1);
                    break;
                case NoSlipWall:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        ulx(i, j, k, n) = 0.0;
                        urx(i, j, k, n) = 0.0;
                    }
                    break;
                case Outflow:
                    ulx(i, j, k, 0) = amrex::max(ulx(i, j, k, 0), 0.0);
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        urx(i, j, k, n) = ulx(i, j, k, n);
                    }
                    break;
                case Interior:
                    break;
            }
        }

        // No need to compute uimh(:,:,0) since it's equal to utrans-w0
        // upwind using full velocity to get transverse component of uimhx
        // Note: utrans already contains w0
        uimhx(i, j, k, 1) =
            utrans(i, j, k) > 0.0 ? ulx(i, j, k, 1) : urx(i, j, k, 1);
        uimhx(i, j, k, 1) = amrex::Math::abs(utrans(i, j, k)) < rel_eps_local
                                ? 0.5 * (ulx(i, j, k, 1) + urx(i, j, k, 1))
                                : uimhx(i, j, k, 1);
    });

    // y-direction
    bclo = phys_bc[1];
    bchi = phys_bc[AMREX_SPACEDIM + 1];

    ParallelFor(mybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        if (ppm_type == 0) {
            Real maxu = amrex::max(0.0, ufull(i, j - 1, k, 1));
            Real minu = amrex::min(0.0, ufull(i, j, k, 1));
            // extrapolate both components of velocity to left face
            uly(i, j, k, 0) = utilde(i, j - 1, k, 0) +
                              (0.5 - (dt2 / hy) * maxu) * Imv(i, j - 1, k, 0);
            uly(i, j, k, 1) = utilde(i, j - 1, k, 1) +
                              (0.5 - (dt2 / hy) * maxu) * Imv(i, j - 1, k, 1);
            // extrapolate both components of velocity to right face
            ury(i, j, k, 0) = utilde(i, j, k, 0) -
                              (0.5 + (dt2 / hy) * minu) * Imv(i, j, k, 0);
            ury(i, j, k, 1) = utilde(i, j, k, 1) -
                              (0.5 + (dt2 / hy) * minu) * Imv(i, j, k, 1);
        } else if (ppm_type == 1 || ppm_type == 2) {
            // extrapolate both components of velocity to left face
            uly(i, j, k, 0) = Ipu(i, j - 1, k, 1);
            uly(i, j, k, 1) = Ipv(i, j - 1, k, 1);
            // extrapolate both components of velocity to right face
            ury(i, j, k, 0) = Imu(i, j, k, 1);
            ury(i, j, k, 1) = Imv(i, j, k, 1);
        }

        // impose lo side bc's
        if (j == domlo[1]) {
            switch (bclo) {
                case Inflow:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        uly(i, j, k, n) = utilde(i, j - 1, k, n);
                        ury(i, j, k, n) = utilde(i, j - 1, k, n);
                    }
                    break;
                case SlipWall:
                case Symmetry:
                    uly(i, j, k, 0) = ury(i, j, k, 0);
                    uly(i, j, k, 1) = 0.0;
                    ury(i, j, k, 1) = 0.0;
                    break;
                case NoSlipWall:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        uly(i, j, k, n) = 0.0;
                        ury(i, j, k, n) = 0.0;
                    }
                    break;
                case Outflow:
                    ury(i, j, k, 1) = amrex::min(ury(i, j, k, 1), 0.0);
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        uly(i, j, k, n) = ury(i, j, k, n);
                    }
                    break;
                case Interior:
                    break;
            }

            // impose hi side bc's
        } else if (j == domhi[1] + 1) {
            switch (bchi) {
                case Inflow:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        uly(i, j, k, n) = utilde(i, j, k, n);
                        ury(i, j, k, n) = utilde(i, j, k, n);
                    }
                    break;
                case SlipWall:
                case Symmetry:
                    ury(i, j, k, 0) = uly(i, j, k, 0);
                    uly(i, j, k, 1) = 0.0;
                    ury(i, j, k, 1) = 0.0;
                    break;
                case NoSlipWall:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        uly(i, j, k, n) = 0.0;
                        ury(i, j, k, n) = 0.0;
                    }
                    break;
                case Outflow:
                    uly(i, j, k, 1) = amrex::max(uly(i, j, k, 1), 0.0);
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        ury(i, j, k, n) = uly(i, j, k, n);
                    }
                    break;
                case Interior:
                    break;
            }
        }
        // No need to compute uimh(:,:,1) since it's equal to utrans-w0
        // upwind using full velocity to get transverse component of uimhy
        // Note: utrans already contains w0
        uimhy(i, j, k, 0) =
            vtrans(i, j, k) > 0.0 ? uly(i, j, k, 0) : ury(i, j, k, 0);
        uimhy(i, j, k, 0) = amrex::Math::abs(vtrans(i, j, k)) < rel_eps_local
                                ? 0.5 * (uly(i, j, k, 0) + ury(i, j, k, 0))
                                : uimhy(i, j, k, 0);
    });
}

void Maestro::VelPredVelocities(
    const MFIter& mfi, Array4<const Real> const utilde,
    Array4<const Real> const utrans, Array4<const Real> const vtrans,
    Array4<Real> const umac, Array4<Real> const vmac,
    Array4<const Real> const Imfx, Array4<const Real> const Ipfx,
    Array4<const Real> const Imfy, Array4<const Real> const Ipfy,
    Array4<const Real> const ulx, Array4<const Real> const urx,
    Array4<const Real> const uimhx, Array4<const Real> const uly,
    Array4<const Real> const ury, Array4<const Real> const uimhy,
    Array4<const Real> const force, Array4<const Real> const w0_cart_in,
    const Box& domainBox, const amrex::GpuArray<Real, AMREX_SPACEDIM> dx) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::VelPredVelocities()", VelPredVelocities);

    //******************************************************************
    // Create umac and vmac
    //******************************************************************

    const Box& xbx = mfi.nodaltilebox(0);
    const Box& ybx = mfi.nodaltilebox(1);

    const Real dt2 = 0.5 * dt;
    const Real dt4 = 0.25 * dt;

    const Real hx = dx[0];
    const Real hy = dx[1];

    const auto rel_eps_local = rel_eps;

    // x-direction
    const auto domlo = domainBox.loVect3d();
    const auto domhi = domainBox.hiVect3d();

    int bclo = phys_bc[0];
    int bchi = phys_bc[AMREX_SPACEDIM];

    ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        // use the traced force if ppm_trace_forces = 1
        Real fl = ppm_trace_forces == 0 ? force(i - 1, j, k, 0)
                                        : Ipfx(i - 1, j, k, 0);
        Real fr = ppm_trace_forces == 0 ? force(i, j, k, 0) : Imfx(i, j, k, 0);

        // extrapolate to edges
        Real umacl = ulx(i, j, k, 0) -
                     (dt4 / hy) *
                         (vtrans(i - 1, j + 1, k) + vtrans(i - 1, j, k)) *
                         (uimhy(i - 1, j + 1, k, 0) - uimhy(i - 1, j, k, 0)) +
                     dt2 * fl;
        Real umacr = urx(i, j, k, 0) -
                     (dt4 / hy) * (vtrans(i, j + 1, k) + vtrans(i, j, k)) *
                         (uimhy(i, j + 1, k, 0) - uimhy(i, j, k, 0)) +
                     dt2 * fr;

        // solve Riemann problem using full velocity
        umac(i, j, k) = 0.5 * (umacl + umacr) > 0.0 ? umacl : umacr;
        umac(i, j, k) =
            (umacl <= 0.0 && umacr >= 0.0) ||
                    (amrex::Math::abs(umacl + umacr) < rel_eps_local)
                ? 0.0
                : umac(i, j, k);

        // impose lo side bc's
        if (i == domlo[0]) {
            switch (bclo) {
                case Inflow:
                    umac(i, j, k) = utilde(i - 1, j, k, 0);
                    break;
                case SlipWall:
                case NoSlipWall:
                case Symmetry:
                    umac(i, j, k) = 0.0;
                    break;
                case Outflow:
                    umac(i, j, k) = amrex::min(umacr, 0.0);
                    break;
                case Interior:
                    break;
            }

            // impose hi side bc's
        } else if (i == domhi[0] + 1) {
            switch (bchi) {
                case Inflow:
                    umac(i, j, k) = utilde(i, j, k, 0);
                    break;
                case SlipWall:
                case NoSlipWall:
                case Symmetry:
                    umac(i, j, k) = 0.0;
                    break;
                case Outflow:
                    umac(i, j, k) = amrex::max(umacl, 0.0);
                    break;
                case Interior:
                    break;
            }
        }
    });

    // y-direction
    bclo = phys_bc[1];
    bchi = phys_bc[AMREX_SPACEDIM + 1];

    ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        // use the traced force if ppm_trace_forces = 1
        Real fl = ppm_trace_forces == 0 ? force(i, j - 1, k, 1)
                                        : Ipfy(i, j - 1, k, 1);
        Real fr = ppm_trace_forces == 0 ? force(i, j, k, 1) : Imfy(i, j, k, 1);

        // extrapolate to edges
        Real vmacl = uly(i, j, k, 1) -
                     (dt4 / hx) *
                         (utrans(i + 1, j - 1, k) + utrans(i, j - 1, k)) *
                         (uimhx(i + 1, j - 1, k, 1) - uimhx(i, j - 1, k, 1)) +
                     dt2 * fl;
        Real vmacr = ury(i, j, k, 1) -
                     (dt4 / hx) * (utrans(i + 1, j, k) + utrans(i, j, k)) *
                         (uimhx(i + 1, j, k, 1) - uimhx(i, j, k, 1)) +
                     dt2 * fr;

        // solve Riemann problem using full velocity
        bool test =
            (vmacl + w0_cart_in(i, j, k, AMREX_SPACEDIM - 1) <= 0.0 &&
             vmacr + w0_cart_in(i, j, k, AMREX_SPACEDIM - 1) >= 0.0) ||
            (amrex::Math::abs(vmacl + vmacr +
                              2 * w0_cart_in(i, j, k, AMREX_SPACEDIM - 1)) <
             rel_eps_local);
        vmac(i, j, k) =
            0.5 * (vmacl + vmacr) + w0_cart_in(i, j, k, AMREX_SPACEDIM - 1) >
                    0.0
                ? vmacl
                : vmacr;
        vmac(i, j, k) = test ? 0.0 : vmac(i, j, k);

        // impose lo side bc's
        if (j == domlo[1]) {
            switch (bclo) {
                case Inflow:
                    vmac(i, j, k) = utilde(i, j - 1, k, 1);
                    break;
                case SlipWall:
                case NoSlipWall:
                case Symmetry:
                    vmac(i, j, k) = 0.0;
                    break;
                case Outflow:
                    vmac(i, j, k) = amrex::min(vmacr, 0.0);
                    break;
                case Interior:
                    break;
            }

            // impose hi side bc's
        } else if (j == domhi[1] + 1) {
            switch (bchi) {
                case Inflow:
                    vmac(i, j, k) = utilde(i, j, k, 1);
                    break;
                case SlipWall:
                case NoSlipWall:
                case Symmetry:
                    vmac(i, j, k) = 0.0;
                    break;
                case Outflow:
                    vmac(i, j, k) = amrex::max(vmacl, 0.0);
                    break;
                case Interior:
                    break;
            }
        }
    });
}

#else

void Maestro::VelPredInterface(
    const MFIter& mfi, Array4<const Real> const utilde,
    Array4<const Real> const ufull, Array4<const Real> const utrans,
    Array4<const Real> const vtrans, Array4<const Real> const wtrans,
    Array4<const Real> const Imu, Array4<const Real> const Ipu,
    Array4<const Real> const Imv, Array4<const Real> const Ipv,
    Array4<const Real> const Imw, Array4<const Real> const Ipw,
    Array4<Real> const ulx, Array4<Real> const urx, Array4<Real> const uimhx,
    Array4<Real> const uly, Array4<Real> const ury, Array4<Real> const uimhy,
    Array4<Real> const ulz, Array4<Real> const urz, Array4<Real> const uimhz,
    const Box& domainBox, const amrex::GpuArray<Real, AMREX_SPACEDIM> dx) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::VelPredInterface()", VelPredInterface);

    ////////////////////////////////////
    // Create u_{\i-\half\e_x}^x, etc.
    ////////////////////////////////////

    // NOTE: for ppm_type == 0, slopex == Ipu, slopey == Imv

    // normal predictor states
    // Allocated from lo:hi+1 in the normal direction
    // lo-1:hi+1 in the transverse directions

    const Box& obx = amrex::grow(mfi.tilebox(), 1);
    const Box& mxbx = amrex::growLo(obx, 0, -1);
    const Box& mybx = amrex::growLo(obx, 1, -1);
    const Box& mzbx = amrex::growLo(obx, 2, -1);

    const Real dt2 = 0.5 * dt;

    const Real hx = dx[0];
    const Real hy = dx[1];
    const Real hz = dx[2];

    const auto rel_eps_local = rel_eps;

    // x-direction
    const auto domlo = domainBox.loVect3d();
    const auto domhi = domainBox.hiVect3d();

    int bclo = phys_bc[0];
    int bchi = phys_bc[AMREX_SPACEDIM];

    ParallelFor(mxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        if (ppm_type == 0) {
            Real maxu = 0.5 - dt2 * amrex::max(0.0, ufull(i - 1, j, k, 0)) / hx;
            Real minu = 0.5 + dt2 * amrex::min(0.0, ufull(i, j, k, 0)) / hx;

            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                // extrapolate all components of velocity to left face
                ulx(i, j, k, n) =
                    utilde(i - 1, j, k, n) + maxu * Ipu(i - 1, j, k, n);
                // extrapolate all components of velocity to right face
                urx(i, j, k, n) = utilde(i, j, k, n) - minu * Ipu(i, j, k, n);
            }
        } else if (ppm_type == 1 || ppm_type == 2) {
            // extrapolate all components of velocity to left face
            ulx(i, j, k, 0) = Ipu(i - 1, j, k, 0);
            ulx(i, j, k, 1) = Ipv(i - 1, j, k, 0);
            ulx(i, j, k, 2) = Ipw(i - 1, j, k, 0);

            // extrapolate all components of velocity to right face
            urx(i, j, k, 0) = Imu(i, j, k, 0);
            urx(i, j, k, 1) = Imv(i, j, k, 0);
            urx(i, j, k, 2) = Imw(i, j, k, 0);
        }

        // impose lo side bc's
        if (i == domlo[0]) {
            switch (bclo) {
                case Inflow:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        ulx(i, j, k, n) = utilde(i - 1, j, k, n);
                        urx(i, j, k, n) = utilde(i - 1, j, k, n);
                    }
                    break;
                case SlipWall:
                case Symmetry:
                    ulx(i, j, k, 0) = 0.0;
                    urx(i, j, k, 0) = 0.0;
                    ulx(i, j, k, 1) = urx(i, j, k, 1);
                    ulx(i, j, k, 2) = urx(i, j, k, 2);
                    break;
                case NoSlipWall:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        ulx(i, j, k, n) = 0.0;
                        urx(i, j, k, n) = 0.0;
                    }
                    break;
                case Outflow:
                    urx(i, j, k, 0) = amrex::min(urx(i, j, k, 0), 0.0);
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        ulx(i, j, k, n) = urx(i, j, k, n);
                    }
                    break;
                case Interior:
                    break;
            }

            // impose hi side bc's
        } else if (i == domhi[0] + 1) {
            switch (bchi) {
                case Inflow:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        ulx(i, j, k, n) = utilde(i, j, k, n);
                        urx(i, j, k, n) = utilde(i, j, k, n);
                    }
                    break;
                case SlipWall:
                case Symmetry:
                    ulx(i, j, k, 0) = 0.0;
                    urx(i, j, k, 0) = 0.0;
                    urx(i, j, k, 1) = ulx(i, j, k, 1);
                    urx(i, j, k, 2) = ulx(i, j, k, 2);
                    break;
                case NoSlipWall:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        ulx(i, j, k, n) = 0.0;
                        urx(i, j, k, n) = 0.0;
                    }
                    break;
                case Outflow:
                    ulx(i, j, k, 0) = amrex::max(ulx(i, j, k, 0), 0.0);
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        urx(i, j, k, n) = ulx(i, j, k, n);
                    }
                    break;
                case Interior:
                    break;
            }
        }

        // No need to compute uimhx(:,:,:,0) since it's equal to utrans-w0
        // upwind using full velocity to get transverse components of uimhx
        // Note: utrans already contains w0
        uimhx(i, j, k, 1) =
            utrans(i, j, k) > 0.0 ? ulx(i, j, k, 1) : urx(i, j, k, 1);
        uimhx(i, j, k, 1) = amrex::Math::abs(utrans(i, j, k)) < rel_eps_local
                                ? 0.5 * (ulx(i, j, k, 1) + urx(i, j, k, 1))
                                : uimhx(i, j, k, 1);

        uimhx(i, j, k, 2) =
            utrans(i, j, k) > 0.0 ? ulx(i, j, k, 2) : urx(i, j, k, 2);
        uimhx(i, j, k, 2) = amrex::Math::abs(utrans(i, j, k)) < rel_eps_local
                                ? 0.5 * (ulx(i, j, k, 2) + urx(i, j, k, 2))
                                : uimhx(i, j, k, 2);
    });

    // y-direction
    bclo = phys_bc[1];
    bchi = phys_bc[AMREX_SPACEDIM + 1];

    ParallelFor(mybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        if (ppm_type == 0) {
            Real maxu =
                (0.5 - dt2 * amrex::max(0.0, ufull(i, j - 1, k, 1)) / hy);
            Real minu = (0.5 + dt2 * amrex::min(0.0, ufull(i, j, k, 1)) / hy);

            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                // extrapolate all components of velocity to left face
                uly(i, j, k, n) =
                    utilde(i, j - 1, k, n) + maxu * Imv(i, j - 1, k, n);
                // extrapolate all components of velocity to right face
                ury(i, j, k, n) = utilde(i, j, k, n) - minu * Imv(i, j, k, n);
            }

        } else if (ppm_type == 1 || ppm_type == 2) {
            // extrapolate all components of velocity to left face
            uly(i, j, k, 0) = Ipu(i, j - 1, k, 1);
            uly(i, j, k, 1) = Ipv(i, j - 1, k, 1);
            uly(i, j, k, 2) = Ipw(i, j - 1, k, 1);

            // extrapolate all components of velocity to right face
            ury(i, j, k, 0) = Imu(i, j, k, 1);
            ury(i, j, k, 1) = Imv(i, j, k, 1);
            ury(i, j, k, 2) = Imw(i, j, k, 1);
        }

        // impose lo side bc's
        if (j == domlo[1]) {
            switch (bclo) {
                case Inflow:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        uly(i, j, k, n) = utilde(i, j - 1, k, n);
                        ury(i, j, k, n) = utilde(i, j - 1, k, n);
                    }
                    break;
                case SlipWall:
                case Symmetry:
                    uly(i, j, k, 0) = ury(i, j, k, 0);
                    uly(i, j, k, 1) = 0.0;
                    ury(i, j, k, 1) = 0.0;
                    uly(i, j, k, 2) = ury(i, j, k, 2);
                    break;
                case NoSlipWall:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        uly(i, j, k, n) = 0.0;
                        ury(i, j, k, n) = 0.0;
                    }
                    break;
                case Outflow:
                    ury(i, j, k, 1) = amrex::min(ury(i, j, k, 1), 0.0);
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        uly(i, j, k, n) = ury(i, j, k, n);
                    }
                    break;
                case Interior:
                    break;
            }

            // impose hi side bc's
        } else if (j == domhi[1] + 1) {
            switch (bchi) {
                case Inflow:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        uly(i, j, k, n) = utilde(i, j, k, n);
                        ury(i, j, k, n) = utilde(i, j, k, n);
                    }
                    break;
                case SlipWall:
                case Symmetry:
                    ury(i, j, k, 0) = uly(i, j, k, 0);
                    uly(i, j, k, 1) = 0.0;
                    ury(i, j, k, 1) = 0.0;
                    ury(i, j, k, 2) = uly(i, j, k, 2);
                    break;
                case NoSlipWall:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        uly(i, j, k, n) = 0.0;
                        ury(i, j, k, n) = 0.0;
                    }
                    break;
                case Outflow:
                    uly(i, j, k, 1) = amrex::max(uly(i, j, k, 1), 0.0);
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        ury(i, j, k, n) = uly(i, j, k, n);
                    }
                    break;
                case Interior:
                    break;
            }
        }

        // No need to compute uimhy(:,:,:,1) since it's equal to vtrans-w0
        // upwind using full velocity to get transverse components of uimhy
        // Note: vtrans already contains w0
        uimhy(i, j, k, 0) =
            vtrans(i, j, k) > 0.0 ? uly(i, j, k, 0) : ury(i, j, k, 0);
        uimhy(i, j, k, 0) = amrex::Math::abs(vtrans(i, j, k)) < rel_eps_local
                                ? 0.5 * (uly(i, j, k, 0) + ury(i, j, k, 0))
                                : uimhy(i, j, k, 0);

        uimhy(i, j, k, 2) =
            vtrans(i, j, k) > 0.0 ? uly(i, j, k, 2) : ury(i, j, k, 2);
        uimhy(i, j, k, 2) = amrex::Math::abs(vtrans(i, j, k)) < rel_eps_local
                                ? 0.5 * (uly(i, j, k, 2) + ury(i, j, k, 2))
                                : uimhy(i, j, k, 2);
    });

    // z-direction
    bclo = phys_bc[2];
    bchi = phys_bc[AMREX_SPACEDIM + 2];

    ParallelFor(mzbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        if (ppm_type == 0) {
            Real maxu = 0.5 - dt2 * amrex::max(0.0, ufull(i, j, k - 1, 2)) / hz;
            Real minu = 0.5 + dt2 * amrex::min(0.0, ufull(i, j, k, 2)) / hz;

            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                // extrapolate all components of velocity to left face
                ulz(i, j, k, n) =
                    utilde(i, j, k - 1, n) + maxu * Imw(i, j, k - 1, n);
                // extrapolate all components of velocity to right face
                urz(i, j, k, n) = utilde(i, j, k, n) - minu * Imw(i, j, k, n);
            }

        } else if (ppm_type == 1 || ppm_type == 2) {
            // extrapolate all components of velocity to left face
            ulz(i, j, k, 0) = Ipu(i, j, k - 1, 2);
            ulz(i, j, k, 1) = Ipv(i, j, k - 1, 2);
            ulz(i, j, k, 2) = Ipw(i, j, k - 1, 2);

            // extrapolate all components of velocity to right face
            urz(i, j, k, 0) = Imu(i, j, k, 2);
            urz(i, j, k, 1) = Imv(i, j, k, 2);
            urz(i, j, k, 2) = Imw(i, j, k, 2);
        }

        // impose lo side bc's
        if (k == domlo[2]) {
            switch (bclo) {
                case Inflow:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        ulz(i, j, k, n) = utilde(i, j, k - 1, n);
                        urz(i, j, k, n) = utilde(i, j, k - 1, n);
                    }
                    break;
                case SlipWall:
                case Symmetry:
                    ulz(i, j, k, 0) = urz(i, j, k, 0);
                    ulz(i, j, k, 1) = urz(i, j, k, 1);
                    ulz(i, j, k, 2) = 0.0;
                    urz(i, j, k, 2) = 0.0;
                    break;
                case NoSlipWall:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        ulz(i, j, k, n) = 0.0;
                        urz(i, j, k, n) = 0.0;
                    }
                    break;
                case Outflow:
                    urz(i, j, k, 2) = amrex::min(urz(i, j, k, 2), 0.0);
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        ulz(i, j, k, n) = urz(i, j, k, n);
                    }
                    break;
                case Interior:
                    break;
            }

            // impose hi side bc's
        } else if (k == domhi[2] + 1) {
            switch (bchi) {
                case Inflow:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        ulz(i, j, k, n) = utilde(i, j, k, n);
                        urz(i, j, k, n) = utilde(i, j, k, n);
                    }
                    break;
                case SlipWall:
                case Symmetry:
                    urz(i, j, k, 0) = ulz(i, j, k, 0);
                    urz(i, j, k, 1) = ulz(i, j, k, 1);
                    ulz(i, j, k, 2) = 0.0;
                    urz(i, j, k, 2) = 0.0;
                    break;
                case NoSlipWall:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        ulz(i, j, k, n) = 0.0;
                        urz(i, j, k, n) = 0.0;
                    }
                    break;
                case Outflow:
                    ulz(i, j, k, 2) = amrex::max(ulz(i, j, k, 2), 0.0);
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        urz(i, j, k, n) = ulz(i, j, k, n);
                    }
                    break;
                case Interior:
                    break;
            }
        }

        // No need to compute uimhz(:,:,:,2) since it's equal to wtrans-w0
        // upwind using full velocity to get transverse components of uimhz
        // Note: wtrans already contains w0
        uimhz(i, j, k, 0) =
            wtrans(i, j, k) > 0.0 ? ulz(i, j, k, 0) : urz(i, j, k, 0);
        uimhz(i, j, k, 0) = amrex::Math::abs(wtrans(i, j, k)) < rel_eps_local
                                ? 0.5 * (ulz(i, j, k, 0) + urz(i, j, k, 0))
                                : uimhz(i, j, k, 0);

        uimhz(i, j, k, 1) =
            wtrans(i, j, k) > 0.0 ? ulz(i, j, k, 1) : urz(i, j, k, 1);
        uimhz(i, j, k, 1) = amrex::Math::abs(wtrans(i, j, k)) < rel_eps_local
                                ? 0.5 * (ulz(i, j, k, 1) + urz(i, j, k, 1))
                                : uimhz(i, j, k, 1);
    });
}

void Maestro::VelPredTransverse(
    const MFIter& mfi, Array4<const Real> const utilde,
    Array4<const Real> const utrans, Array4<const Real> const vtrans,
    Array4<const Real> const wtrans, Array4<const Real> const ulx,
    Array4<const Real> const urx, Array4<const Real> const uimhx,
    Array4<const Real> const uly, Array4<const Real> const ury,
    Array4<const Real> const uimhy, Array4<const Real> const ulz,
    Array4<const Real> const urz, Array4<const Real> const uimhz,
    Array4<Real> const uimhyz, Array4<Real> const uimhzy,
    Array4<Real> const vimhxz, Array4<Real> const vimhzx,
    Array4<Real> const wimhxy, Array4<Real> const wimhyx, const Box& domainBox,
    const amrex::GpuArray<Real, AMREX_SPACEDIM> dx) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::VelPredTransverse()", VelPredTransverse);

    //////////////////////////////////////
    // Create u_{\i-\half\e_y}^{y|z}, etc.
    //////////////////////////////////////

    const Real dt6 = dt / 6.0;

    const Real hx = dx[0];
    const Real hy = dx[1];
    const Real hz = dx[2];
    const auto domlo = domainBox.loVect3d();
    const auto domhi = domainBox.hiVect3d();
    const auto rel_eps_local = rel_eps;

    GpuArray<int, AMREX_SPACEDIM * 2> physbc;
    for (int n = 0; n < AMREX_SPACEDIM * 2; ++n) {
        physbc[n] = phys_bc[n];
    }

    // uimhyz, 1, 2
    Box imhbox = amrex::grow(mfi.tilebox(), 0, 1);
    imhbox = amrex::growHi(imhbox, 1, 1);

    ParallelFor(imhbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        // extrapolate to faces
        Real ulyz = uly(i, j, k, 0) -
                    (dt6 / hz) *
                        (wtrans(i, j - 1, k + 1) + wtrans(i, j - 1, k)) *
                        (uimhz(i, j - 1, k + 1, 0) - uimhz(i, j - 1, k, 0));
        Real uryz = ury(i, j, k, 0) -
                    (dt6 / hz) * (wtrans(i, j, k + 1) + wtrans(i, j, k)) *
                        (uimhz(i, j, k + 1, 0) - uimhz(i, j, k, 0));

        // impose lo side bc's
        if (j == domlo[1]) {
            switch (physbc[1]) {
                case Inflow:
                    ulyz = utilde(i, j - 1, k, 0);
                    uryz = utilde(i, j - 1, k, 0);
                    break;
                case SlipWall:
                case Symmetry:
                case Outflow:
                    ulyz = uryz;
                    break;
                case NoSlipWall:
                    ulyz = 0.0;
                    uryz = 0.0;
                    break;
                case Interior:
                    break;
            }

            // impose hi side bc's
        } else if (j == domhi[1] + 1) {
            switch (physbc[AMREX_SPACEDIM + 1]) {
                case Inflow:
                    ulyz = utilde(i, j, k, 0);
                    uryz = utilde(i, j, k, 0);
                    break;
                case SlipWall:
                case Symmetry:
                case Outflow:
                    uryz = ulyz;
                    break;
                case NoSlipWall:
                    ulyz = 0.0;
                    uryz = 0.0;
                    break;
                case Interior:
                    break;
            }
        }

        // upwind using full velocity
        uimhyz(i, j, k) = vtrans(i, j, k) > 0.0 ? ulyz : uryz;
        uimhyz(i, j, k) = amrex::Math::abs(vtrans(i, j, k)) < rel_eps_local
                              ? 0.5 * (ulyz + uryz)
                              : uimhyz(i, j, k);
    });

    // uimhzy, 1, 3
    imhbox = amrex::grow(mfi.tilebox(), 0, 1);
    imhbox = amrex::growHi(imhbox, 2, 1);

    ParallelFor(imhbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        // extrapolate to faces
        Real ulzy = ulz(i, j, k, 0) -
                    (dt6 / hy) *
                        (vtrans(i, j + 1, k - 1) + vtrans(i, j, k - 1)) *
                        (uimhy(i, j + 1, k - 1, 0) - uimhy(i, j, k - 1, 0));
        Real urzy = urz(i, j, k, 0) -
                    (dt6 / hy) * (vtrans(i, j + 1, k) + vtrans(i, j, k)) *
                        (uimhy(i, j + 1, k, 0) - uimhy(i, j, k, 0));

        // impose lo side bc's
        if (k == domlo[2]) {
            switch (physbc[2]) {
                case Inflow:
                    ulzy = utilde(i, j, k - 1, 0);
                    urzy = utilde(i, j, k - 1, 0);
                    break;
                case SlipWall:
                case Symmetry:
                case Outflow:
                    ulzy = urzy;
                    break;
                case NoSlipWall:
                    ulzy = 0.0;
                    urzy = 0.0;
                    break;
                case Interior:
                    break;
            }

            // impose hi side bc's
        } else if (k == domhi[2] + 1) {
            switch (physbc[AMREX_SPACEDIM + 2]) {
                case Inflow:
                    ulzy = utilde(i, j, k, 0);
                    urzy = utilde(i, j, k, 0);
                    break;
                case SlipWall:
                case Symmetry:
                case Outflow:
                    urzy = ulzy;
                    break;
                case NoSlipWall:
                    ulzy = 0.0;
                    urzy = 0.0;
                    break;
                case Interior:
                    break;
            }
        }

        // upwind using full velocity
        uimhzy(i, j, k) = wtrans(i, j, k) > 0.0 ? ulzy : urzy;
        uimhzy(i, j, k) = amrex::Math::abs(wtrans(i, j, k)) < rel_eps_local
                              ? 0.5 * (ulzy + urzy)
                              : uimhzy(i, j, k);
    });

    // vimhxz, 2, 1
    imhbox = amrex::grow(mfi.tilebox(), 1, 1);
    imhbox = amrex::growHi(imhbox, 0, 1);

    ParallelFor(imhbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        // extrapolate to faces
        Real vlxz = ulx(i, j, k, 1) -
                    (dt6 / hz) *
                        (wtrans(i - 1, j, k + 1) + wtrans(i - 1, j, k)) *
                        (uimhz(i - 1, j, k + 1, 1) - uimhz(i - 1, j, k, 1));
        Real vrxz = urx(i, j, k, 1) -
                    (dt6 / hz) * (wtrans(i, j, k + 1) + wtrans(i, j, k)) *
                        (uimhz(i, j, k + 1, 1) - uimhz(i, j, k, 1));

        // impose lo side bc's
        if (i == domlo[0]) {
            switch (physbc[0]) {
                case Inflow:
                    vlxz = utilde(i - 1, j, k, 1);
                    vrxz = utilde(i - 1, j, k, 1);
                    break;
                case SlipWall:
                case Symmetry:
                case Outflow:
                    vlxz = vrxz;
                    break;
                case NoSlipWall:
                    vlxz = 0.0;
                    vrxz = 0.0;
                    break;
                case Interior:
                    break;
            }

            // impose hi side bc's
        } else if (i == domhi[0] + 1) {
            switch (physbc[AMREX_SPACEDIM]) {
                case Inflow:
                    vlxz = utilde(i, j, k, 1);
                    vrxz = utilde(i, j, k, 1);
                    break;
                case SlipWall:
                case Symmetry:
                case Outflow:
                    vrxz = vlxz;
                    break;
                case NoSlipWall:
                    vlxz = 0.0;
                    vrxz = 0.0;
                    break;
                case Interior:
                    break;
            }
        }

        // upwind using full velocity
        vimhxz(i, j, k) = utrans(i, j, k) > 0.0 ? vlxz : vrxz;
        vimhxz(i, j, k) = amrex::Math::abs(utrans(i, j, k)) < rel_eps_local
                              ? 0.5 * (vlxz + vrxz)
                              : vimhxz(i, j, k);
    });

    // vimhzx, 2, 3
    imhbox = amrex::grow(mfi.tilebox(), 1, 1);
    imhbox = amrex::growHi(imhbox, 2, 1);

    ParallelFor(imhbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        // extrapolate to faces
        Real vlzx = ulz(i, j, k, 1) -
                    (dt6 / hx) *
                        (utrans(i + 1, j, k - 1) + utrans(i, j, k - 1)) *
                        (uimhx(i + 1, j, k - 1, 1) - uimhx(i, j, k - 1, 1));
        Real vrzx = urz(i, j, k, 1) -
                    (dt6 / hx) * (utrans(i + 1, j, k) + utrans(i, j, k)) *
                        (uimhx(i + 1, j, k, 1) - uimhx(i, j, k, 1));

        // impose lo side bc's
        if (k == domlo[2]) {
            switch (physbc[2]) {
                case Inflow:
                    vlzx = utilde(i, j, k - 1, 1);
                    vrzx = utilde(i, j, k - 1, 1);
                    break;
                case SlipWall:
                case Symmetry:
                case Outflow:
                    vlzx = vrzx;
                    break;
                case NoSlipWall:
                    vlzx = 0.0;
                    vrzx = 0.0;
                    break;
                case Interior:
                    break;
            }

            // impose hi side bc's
        } else if (k == domhi[2] + 1) {
            switch (physbc[AMREX_SPACEDIM + 2]) {
                case Inflow:
                    vlzx = utilde(i, j, k, 1);
                    vrzx = utilde(i, j, k, 1);
                    break;
                case SlipWall:
                case Symmetry:
                case Outflow:
                    vrzx = vlzx;
                    break;
                case NoSlipWall:
                    vlzx = 0.0;
                    vrzx = 0.0;
                    break;
                case Interior:
                    break;
            }
        }

        // upwind using full velocity
        vimhzx(i, j, k) = wtrans(i, j, k) > 0.0 ? vlzx : vrzx;
        vimhzx(i, j, k) = amrex::Math::abs(wtrans(i, j, k)) < rel_eps_local
                              ? 0.5 * (vlzx + vrzx)
                              : vimhzx(i, j, k);
    });

    // wimhxy, 3, 1
    imhbox = amrex::grow(mfi.tilebox(), 2, 1);
    imhbox = amrex::growHi(imhbox, 0, 1);

    ParallelFor(imhbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        // extrapolate to faces
        Real wlxy = ulx(i, j, k, 2) -
                    (dt6 / hy) *
                        (vtrans(i - 1, j + 1, k) + vtrans(i - 1, j, k)) *
                        (uimhy(i - 1, j + 1, k, 2) - uimhy(i - 1, j, k, 2));
        Real wrxy = urx(i, j, k, 2) -
                    (dt6 / hy) * (vtrans(i, j + 1, k) + vtrans(i, j, k)) *
                        (uimhy(i, j + 1, k, 2) - uimhy(i, j, k, 2));

        // impose lo side bc's
        if (i == domlo[0]) {
            switch (physbc[0]) {
                case Inflow:
                    wlxy = utilde(i - 1, j, k, 2);
                    wrxy = utilde(i - 1, j, k, 2);
                    break;
                case SlipWall:
                case Symmetry:
                case Outflow:
                    wlxy = wrxy;
                    break;
                case NoSlipWall:
                    wlxy = 0.0;
                    wrxy = 0.0;
                    break;
                case Interior:
                    break;
            }

            // impose hi side bc's
        } else if (i == domhi[0] + 1) {
            switch (physbc[AMREX_SPACEDIM]) {
                case Inflow:
                    wlxy = utilde(i, j, k, 2);
                    wrxy = utilde(i, j, k, 2);
                    break;
                case SlipWall:
                case Symmetry:
                case Outflow:
                    wrxy = wlxy;
                    break;
                case NoSlipWall:
                    wlxy = 0.0;
                    wrxy = 0.0;
                    break;
                case Interior:
                    break;
            }
        }

        // upwind using full velocity
        wimhxy(i, j, k) = utrans(i, j, k) > 0.0 ? wlxy : wrxy;
        wimhxy(i, j, k) = amrex::Math::abs(utrans(i, j, k)) < rel_eps_local
                              ? 0.5 * (wlxy + wrxy)
                              : wimhxy(i, j, k);
    });

    // wimhyx, 3, 2
    imhbox = amrex::grow(mfi.tilebox(), 2, 1);
    imhbox = amrex::growHi(imhbox, 1, 1);

    ParallelFor(imhbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        // extrapolate to faces
        Real wlyx = uly(i, j, k, 2) -
                    (dt6 / hx) *
                        (utrans(i + 1, j - 1, k) + utrans(i, j - 1, k)) *
                        (uimhx(i + 1, j - 1, k, 2) - uimhx(i, j - 1, k, 2));
        Real wryx = ury(i, j, k, 2) -
                    (dt6 / hx) * (utrans(i + 1, j, k) + utrans(i, j, k)) *
                        (uimhx(i + 1, j, k, 2) - uimhx(i, j, k, 2));

        // impose lo side bc's
        if (j == domlo[1]) {
            switch (physbc[1]) {
                case Inflow:
                    wlyx = utilde(i, j - 1, k, 2);
                    wryx = utilde(i, j - 1, k, 2);
                    break;
                case SlipWall:
                case Symmetry:
                case Outflow:
                    wlyx = wryx;
                    break;
                case NoSlipWall:
                    wlyx = 0.0;
                    wryx = 0.0;
                    break;
                case Interior:
                    break;
            }

            // impose hi side bc's
        } else if (j == domhi[1] + 1) {
            switch (physbc[AMREX_SPACEDIM + 1]) {
                case Inflow:
                    wlyx = utilde(i, j, k, 2);
                    wryx = utilde(i, j, k, 2);
                    break;
                case SlipWall:
                case Symmetry:
                case Outflow:
                    wryx = wlyx;
                    break;
                case NoSlipWall:
                    wlyx = 0.0;
                    wryx = 0.0;
                    break;
                case Interior:
                    break;
            }
        }

        // upwind using full velocity
        wimhyx(i, j, k) = vtrans(i, j, k) > 0.0 ? wlyx : wryx;
        wimhyx(i, j, k) = amrex::Math::abs(vtrans(i, j, k)) < rel_eps_local
                              ? 0.5 * (wlyx + wryx)
                              : wimhyx(i, j, k);
    });
}

void Maestro::VelPredVelocities(
    const MFIter& mfi, Array4<const Real> const utilde,
    Array4<const Real> const utrans, Array4<const Real> const vtrans,
    Array4<const Real> const wtrans, Array4<Real> const umac,
    Array4<Real> const vmac, Array4<Real> const wmac,
    Array4<const Real> const w0macx, Array4<const Real> const w0macy,
    Array4<const Real> const w0macz, Array4<const Real> const Imfx,
    Array4<const Real> const Ipfx, Array4<const Real> const Imfy,
    Array4<const Real> const Ipfy, Array4<const Real> const Imfz,
    Array4<const Real> const Ipfz, Array4<const Real> const ulx,
    Array4<const Real> const urx, Array4<const Real> const uly,
    Array4<const Real> const ury, Array4<const Real> const ulz,
    Array4<const Real> const urz, Array4<const Real> const uimhyz,
    Array4<const Real> const uimhzy, Array4<const Real> const vimhxz,
    Array4<const Real> const vimhzx, Array4<const Real> const wimhxy,
    Array4<const Real> const wimhyx, Array4<const Real> const force,
    Array4<const Real> const w0_cart_in, const Box& domainBox,
    const amrex::GpuArray<Real, AMREX_SPACEDIM> dx) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::VelPredVelocities()", VelPredVelocities);

    //******************************************************************
    // Create umac and vmac
    //******************************************************************

    const Box& xbx = mfi.nodaltilebox(0);
    const Box& ybx = mfi.nodaltilebox(1);
    const Box& zbx = mfi.nodaltilebox(2);

    const Real dt2 = 0.5 * dt;
    const Real dt4 = 0.25 * dt;

    const Real hx = dx[0];
    const Real hy = dx[1];
    const Real hz = dx[2];

    const auto rel_eps_local = rel_eps;

    const auto domlo = domainBox.loVect3d();
    const auto domhi = domainBox.hiVect3d();

    GpuArray<int, AMREX_SPACEDIM * 2> physbc;
    for (int n = 0; n < AMREX_SPACEDIM * 2; ++n) {
        physbc[n] = phys_bc[n];
    }

    // x-direction
    ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        // use the traced force if ppm_trace_forces = 1
        Real fl = ppm_trace_forces == 0 ? force(i - 1, j, k, 0)
                                        : Ipfx(i - 1, j, k, 0);
        Real fr = ppm_trace_forces == 0 ? force(i, j, k, 0) : Imfx(i, j, k, 0);

        // extrapolate to edges
        Real umacl =
            ulx(i, j, k, 0) -
            (dt4 / hy) * (vtrans(i - 1, j + 1, k) + vtrans(i - 1, j, k)) *
                (uimhyz(i - 1, j + 1, k) - uimhyz(i - 1, j, k)) -
            (dt4 / hz) * (wtrans(i - 1, j, k + 1) + wtrans(i - 1, j, k)) *
                (uimhzy(i - 1, j, k + 1) - uimhzy(i - 1, j, k)) +
            dt2 * fl;
        Real umacr = urx(i, j, k, 0) -
                     (dt4 / hy) * (vtrans(i, j + 1, k) + vtrans(i, j, k)) *
                         (uimhyz(i, j + 1, k) - uimhyz(i, j, k)) -
                     (dt4 / hz) * (wtrans(i, j, k + 1) + wtrans(i, j, k)) *
                         (uimhzy(i, j, k + 1) - uimhzy(i, j, k)) +
                     dt2 * fr;

        if (spherical) {
            // solve Riemann problem using full velocity
            bool test =
                (umacl + w0macx(i, j, k) <= 0.0 &&
                 umacr + w0macx(i, j, k) >= 0.0) ||
                (amrex::Math::abs(umacl + umacr + 2.0 * w0macx(i, j, k)) <
                 rel_eps_local);
            umac(i, j, k) =
                0.5 * (umacl + umacr) + w0macx(i, j, k) > 0.0 ? umacl : umacr;
            umac(i, j, k) = test ? 0.0 : umac(i, j, k);
        } else {
            // solve Riemann problem using full velocity
            bool test = (umacl <= 0.0 && umacr >= 0.0) ||
                        (amrex::Math::abs(umacl + umacr) < rel_eps_local);
            umac(i, j, k) = 0.5 * (umacl + umacr) > 0.0 ? umacl : umacr;
            umac(i, j, k) = test ? 0.0 : umac(i, j, k);
        }

        // impose lo side bc's
        if (i == domlo[0]) {
            switch (physbc[0]) {
                case Inflow:
                    umac(i, j, k) = utilde(i - 1, j, k, 0);
                    break;
                case SlipWall:
                case Symmetry:
                case NoSlipWall:
                    umac(i, j, k) = 0.0;
                    break;
                case Outflow:
                    umac(i, j, k) = amrex::min(umacr, 0.0);
                    break;
                case Interior:
                    break;
            }

            // impose hi side bc's
        } else if (i == domhi[0] + 1) {
            switch (physbc[AMREX_SPACEDIM]) {
                case Inflow:
                    umac(i, j, k) = utilde(i, j, k, 0);
                    break;
                case SlipWall:
                case Symmetry:
                case NoSlipWall:
                    umac(i, j, k) = 0.0;
                    break;
                case Outflow:
                    umac(i, j, k) = amrex::max(umacl, 0.0);
                    break;
                case Interior:
                    break;
            }
        }
    });

    // y-direction
    ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        // use the traced force if ppm_trace_forces = 1
        Real fl = ppm_trace_forces == 0 ? force(i, j - 1, k, 1)
                                        : Ipfy(i, j - 1, k, 1);
        Real fr = ppm_trace_forces == 0 ? force(i, j, k, 1) : Imfy(i, j, k, 1);

        // extrapolate to edges
        Real vmacl =
            uly(i, j, k, 1) -
            (dt4 / hx) * (utrans(i + 1, j - 1, k) + utrans(i, j - 1, k)) *
                (vimhxz(i + 1, j - 1, k) - vimhxz(i, j - 1, k)) -
            (dt4 / hz) * (wtrans(i, j - 1, k + 1) + wtrans(i, j - 1, k)) *
                (vimhzx(i, j - 1, k + 1) - vimhzx(i, j - 1, k)) +
            dt2 * fl;
        Real vmacr = ury(i, j, k, 1) -
                     (dt4 / hx) * (utrans(i + 1, j, k) + utrans(i, j, k)) *
                         (vimhxz(i + 1, j, k) - vimhxz(i, j, k)) -
                     (dt4 / hz) * (wtrans(i, j, k + 1) + wtrans(i, j, k)) *
                         (vimhzx(i, j, k + 1) - vimhzx(i, j, k)) +
                     dt2 * fr;

        if (spherical) {
            // solve Riemann problem using full velocity
            bool test =
                (vmacl + w0macy(i, j, k) <= 0.0 &&
                 vmacr + w0macy(i, j, k) >= 0.0) ||
                (amrex::Math::abs(vmacl + vmacr + 2.0 * w0macy(i, j, k)) <
                 rel_eps_local);
            vmac(i, j, k) =
                0.5 * (vmacl + vmacr) + w0macy(i, j, k) > 0.0 ? vmacl : vmacr;
            vmac(i, j, k) = test ? 0.0 : vmac(i, j, k);
        } else {
            // solve Riemann problem using full velocity
            bool test = (vmacl <= 0.0 && vmacr >= 0.0) ||
                        (amrex::Math::abs(vmacl + vmacr) < rel_eps_local);
            vmac(i, j, k) = 0.5 * (vmacl + vmacr) > 0.0 ? vmacl : vmacr;
            vmac(i, j, k) = test ? 0.0 : vmac(i, j, k);
        }

        // impose lo side bc's
        if (j == domlo[1]) {
            switch (physbc[1]) {
                case Inflow:
                    vmac(i, j, k) = utilde(i, j - 1, k, 1);
                    break;
                case SlipWall:
                case Symmetry:
                case NoSlipWall:
                    vmac(i, j, k) = 0.0;
                    break;
                case Outflow:
                    vmac(i, j, k) = amrex::min(vmacr, 0.0);
                    break;
                case Interior:
                    break;
            }

            // impose hi side bc's
        } else if (j == domhi[1] + 1) {
            switch (physbc[AMREX_SPACEDIM + 1]) {
                case Inflow:
                    vmac(i, j, k) = utilde(i, j, k, 1);
                    break;
                case SlipWall:
                case Symmetry:
                case NoSlipWall:
                    vmac(i, j, k) = 0.0;
                    break;
                case Outflow:
                    vmac(i, j, k) = amrex::max(vmacl, 0.0);
                    break;
                case Interior:
                    break;
            }
        }
    });

    // z-direction
    ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        // use the traced force if ppm_trace_forces = 1
        Real fl = ppm_trace_forces == 0 ? force(i, j, k - 1, 2)
                                        : Ipfz(i, j, k - 1, 2);
        Real fr = ppm_trace_forces == 0 ? force(i, j, k, 2) : Imfz(i, j, k, 2);

        // extrapolate to edges
        Real wmacl =
            ulz(i, j, k, 2) -
            (dt4 / hx) * (utrans(i + 1, j, k - 1) + utrans(i, j, k - 1)) *
                (wimhxy(i + 1, j, k - 1) - wimhxy(i, j, k - 1)) -
            (dt4 / hy) * (vtrans(i, j + 1, k - 1) + vtrans(i, j, k - 1)) *
                (wimhyx(i, j + 1, k - 1) - wimhyx(i, j, k - 1)) +
            dt2 * fl;
        Real wmacr = urz(i, j, k, 2) -
                     (dt4 / hx) * (utrans(i + 1, j, k) + utrans(i, j, k)) *
                         (wimhxy(i + 1, j, k) - wimhxy(i, j, k)) -
                     (dt4 / hy) * (vtrans(i, j + 1, k) + vtrans(i, j, k)) *
                         (wimhyx(i, j + 1, k) - wimhyx(i, j, k)) +
                     dt2 * fr;

        if (spherical) {
            // solve Riemann problem using full velocity
            bool test =
                (wmacl + w0macz(i, j, k) <= 0.0 &&
                 wmacr + w0macz(i, j, k) >= 0.0) ||
                (amrex::Math::abs(wmacl + wmacr + 2.0 * w0macz(i, j, k)) <
                 rel_eps_local);
            wmac(i, j, k) =
                0.5 * (wmacl + wmacr) + w0macz(i, j, k) > 0.0 ? wmacl : wmacr;
            wmac(i, j, k) = test ? 0.0 : wmac(i, j, k);
        } else {
            // solve Riemann problem using full velocity
            bool test =
                (wmacl + w0_cart_in(i, j, k, AMREX_SPACEDIM - 1) <= 0.0 &&
                 wmacr + w0_cart_in(i, j, k, AMREX_SPACEDIM - 1) >= 0.0) ||
                (amrex::Math::abs(wmacl + wmacr +
                                  2.0 *
                                      w0_cart_in(i, j, k, AMREX_SPACEDIM - 1)) <
                 rel_eps_local);
            wmac(i, j, k) =
                0.5 * (wmacl + wmacr) +
                            w0_cart_in(i, j, k, AMREX_SPACEDIM - 1) >
                        0.0
                    ? wmacl
                    : wmacr;
            wmac(i, j, k) = test ? 0.0 : wmac(i, j, k);
        }

        // impose hi side bc's
        if (k == domlo[2]) {
            switch (physbc[2]) {
                case Inflow:
                    wmac(i, j, k) = utilde(i, j, k - 1, 2);
                    break;
                case SlipWall:
                case Symmetry:
                case NoSlipWall:
                    wmac(i, j, k) = 0.0;
                    break;
                case Outflow:
                    wmac(i, j, k) = amrex::min(wmacr, 0.0);
                    break;
                case Interior:
                    break;
            }

            // impose lo side bc's
        } else if (k == domhi[2] + 1) {
            switch (physbc[AMREX_SPACEDIM + 2]) {
                case Inflow:
                    wmac(i, j, k) = utilde(i, j, k, 2);
                    break;
                case SlipWall:
                case Symmetry:
                case NoSlipWall:
                    wmac(i, j, k) = 0.0;
                    break;
                case Outflow:
                    wmac(i, j, k) = amrex::max(wmacl, 0.0);
                    break;
                case Interior:
                    break;
            }
        }
    });
}

#endif

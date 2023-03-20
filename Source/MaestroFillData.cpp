
#include <Maestro.H>
#include <Maestro_F.H>
#include <PhysBCFunctMaestro.H>

using namespace amrex;

// call FillPatch for all levels
void Maestro::FillPatch(Real time, Vector<MultiFab>& mf,
                        Vector<MultiFab>& mf_old, Vector<MultiFab>& mf_new,
                        int srccomp, int destcomp, int ncomp, int startbccomp,
                        const Vector<BCRec>& bcs_in, int variable_type) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::FillPatch()", FillPatch);

    for (int lev = 0; lev <= finest_level; ++lev) {
        FillPatch(lev, time, mf[lev], mf_old, mf_new, srccomp, destcomp, ncomp,
                  startbccomp, bcs_in, variable_type);
    }
}

// compute a new multifab by coping in phi from valid region and filling ghost cells
// works for single level and 2-level cases
// (fill fine grid ghost by interpolating from coarse)
// srccomp of the source component
// destcomp is the destination component AND the bc component
void Maestro::FillPatch(int lev, Real time, MultiFab& mf,
                        Vector<MultiFab>& mf_old, Vector<MultiFab>& mf_new,
                        int srccomp, int destcomp, int ncomp, int startbccomp,
                        const Vector<BCRec>& bcs_in, int variable_type) {
    Vector<BCRec> bcs{bcs_in.begin() + startbccomp,
                      bcs_in.begin() + startbccomp + ncomp};

    if (lev == 0) {
        Vector<MultiFab*> smf;
        Vector<Real> stime;
        GetData(0, time, smf, stime, mf_old, mf_new);

        PhysBCFunctMaestro physbc;

        if (variable_type == 1) {  // velocity
            physbc.define(geom[lev], bcs, BndryFuncArrayMaestro(VelFill));
        } else {  // scalar
            physbc.define(geom[lev], bcs, BndryFuncArrayMaestro(ScalarFill));
        }

        FillPatchSingleLevel(mf, time, smf, stime, srccomp, destcomp, ncomp,
                             geom[lev], physbc, 0);
    } else {
        Vector<MultiFab*> cmf, fmf;
        Vector<Real> ctime, ftime;
        GetData(lev - 1, time, cmf, ctime, mf_old, mf_new);
        GetData(lev, time, fmf, ftime, mf_old, mf_new);

        PhysBCFunctMaestro cphysbc;
        PhysBCFunctMaestro fphysbc;

        if (variable_type == 1) {  // velocity
            cphysbc.define(geom[lev - 1], bcs, BndryFuncArrayMaestro(VelFill));
            fphysbc.define(geom[lev], bcs, BndryFuncArrayMaestro(VelFill));
        } else {  // scalar
            cphysbc.define(geom[lev - 1], bcs,
                           BndryFuncArrayMaestro(ScalarFill));
            fphysbc.define(geom[lev], bcs, BndryFuncArrayMaestro(ScalarFill));
        }

        Interpolater* mapper = &cell_cons_interp;
        FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime, srccomp, destcomp,
                           ncomp, geom[lev - 1], geom[lev], cphysbc, 0, fphysbc,
                           0, refRatio(lev - 1), mapper, bcs, 0);
    }
}

// fill an entire multifab by interpolating from the coarser level
// this comes into play when a new level of refinement appears
// srccomp of the source component
// destcomp is the destination component AND the bc component
void Maestro::FillCoarsePatch(int lev, Real time, MultiFab& mf,
                              Vector<MultiFab>& mf_old,
                              Vector<MultiFab>& mf_new, int srccomp,
                              int destcomp, int ncomp, const Vector<BCRec>& bcs,
                              int variable_type) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::FillCoarsePatch()", FillCoarsePatch);

    AMREX_ASSERT(lev > 0);

    Vector<MultiFab*> cmf;
    Vector<Real> ctime;
    GetData(lev - 1, time, cmf, ctime, mf_old, mf_new);

    if (cmf.size() != 1) {
        Abort("FillCoarsePatch: how did this happen?");
    }

    PhysBCFunctMaestro cphysbc;
    PhysBCFunctMaestro fphysbc;

    if (variable_type == 1) {  // velocity
        cphysbc.define(geom[lev - 1], bcs, BndryFuncArrayMaestro(VelFill));
        fphysbc.define(geom[lev], bcs, BndryFuncArrayMaestro(VelFill));
    } else {  // scalar
        cphysbc.define(geom[lev - 1], bcs, BndryFuncArrayMaestro(ScalarFill));
        fphysbc.define(geom[lev], bcs, BndryFuncArrayMaestro(ScalarFill));
    }

    Interpolater* mapper = &cell_cons_interp;

    InterpFromCoarseLevel(mf, time, *cmf[0], srccomp, destcomp, ncomp,
                          geom[lev - 1], geom[lev], cphysbc, 0, fphysbc, 0,
                          refRatio(lev - 1), mapper, bcs, 0);
}

// utility to copy in data from mf_old and/or mf_new into mf
// if time=t_old we copy mf_old into mf
// if time=t_new we copy mf_new into mf
// otherwise copy copy in both mf_old and mf_new into mf and the fillpatch
// routines know to interpolate in time.  However in MAESTRO since we don't
// subcycle I'm not sure if we need this capability?
void Maestro::GetData(int lev, Real time, Vector<MultiFab*>& mf,
                      Vector<Real>& mftime, Vector<MultiFab>& mf_old,
                      Vector<MultiFab>& mf_new) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::GetData()", GetData);

    mf.clear();
    mftime.clear();

    const Real teps = (t_new - t_old) * 1.e-3;

    if (time >= t_new - teps && time <= t_new + teps) {
        mf.push_back(&mf_new[lev]);
        mftime.push_back(t_new);
    } else if (time >= t_old - teps && time <= t_old + teps) {
        mf.push_back(&mf_old[lev]);
        mftime.push_back(t_old);
    } else {
        mf.push_back(&mf_old[lev]);
        mf.push_back(&mf_new[lev]);
        mftime.push_back(t_old);
        mftime.push_back(t_new);
    }
}

// set covered coarse cells to be the average of overlying fine cells
void Maestro::AverageDown(Vector<MultiFab>& mf, int comp, int ncomp) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::AverageDown()", AverageDown);

    for (int lev = finest_level - 1; lev >= 0; --lev) {
        average_down(mf[lev + 1], mf[lev], geom[lev + 1], geom[lev], comp,
                     ncomp, refRatio(lev));
    }
}

// set covered faces to be the average of overlying fine faces
void Maestro::AverageDownFaces(
    Vector<std::array<MultiFab, AMREX_SPACEDIM> >& edge) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::AverageDownFaces()", AverageDownFaces);

    for (int lev = finest_level - 1; lev >= 0; --lev) {
        Vector<const MultiFab*> edge_f(AMREX_SPACEDIM);
        Vector<MultiFab*> edge_c(AMREX_SPACEDIM);

        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            edge_f[dir] = &(edge[lev + 1][dir]);
            edge_c[dir] = &(edge[lev][dir]);
        }

        average_down_faces(edge_f, edge_c, refRatio(lev));
    }
}

// fill in ONE ghost cell for all components of a face-centered (MAC) velocity
// field behind physical boundaries.  Does not modify the velocities on the boundary
void Maestro::FillUmacGhost(
    Vector<std::array<MultiFab, AMREX_SPACEDIM> >& umac_in, int level) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::FillUmacGhost()", FillUmacGhost);

    int start_lev;
    int end_lev;
    if (level == -1) {
        start_lev = 0;
        end_lev = finest_level;
    } else {
        start_lev = level;
        end_lev = level;
    }

    for (int lev = start_lev; lev <= end_lev; ++lev) {
        // Get the index space of the domain
        const Box& domainBox = geom[lev].Domain();
        const auto domlo = domainBox.loVect3d();
        const auto domhi = domainBox.hiVect3d();

        const int* AMREX_RESTRICT physbc_p = phys_bc.dataPtr();

        // get references to the MultiFabs at level lev
        MultiFab& sold_mf =
            sold[lev];  // need a cell-centered MF for the MFIter

        // DO NOT TILE THIS SUBROUTINE
        // this just filling ghost cells so the fortran logic has to be reworked
        // to properly capture the corner terms
        for (MFIter mfi(sold_mf, false); mfi.isValid(); ++mfi) {
            // Get the index space of the valid (cell-centered) region
            const auto xbx = mfi.grownnodaltilebox(0, 1);
            const auto ybx = mfi.grownnodaltilebox(1, 1);
#if (AMREX_SPACEDIM == 3)
            const auto zbx = mfi.grownnodaltilebox(2, 1);
#endif

            const Array4<Real> umac = umac_in[lev][0].array(mfi);
            const Array4<Real> vmac = umac_in[lev][1].array(mfi);
#if (AMREX_SPACEDIM == 3)
            const Array4<Real> wmac = umac_in[lev][2].array(mfi);
#endif

            ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                // lo x-faces
                if (i == domlo[0] - 1) {
                    switch (physbc_p[0]) {
                        case Inflow:
                            umac(i, j, k) = umac(i + 1, j, k);
                            vmac(i, j, k) = 0.0;
#if (AMREX_SPACEDIM == 3)
                            wmac(i, j, k) = 0.0;
#endif
                            break;
                        case SlipWall:
                        case NoSlipWall:
                            umac(i, j, k) = 0.0;
                            vmac(i, j, k) = 0.0;
#if (AMREX_SPACEDIM == 3)
                            wmac(i, j, k) = 0.0;
#endif
                            break;
                        case Outflow:
                            umac(i, j, k) = umac(i + 1, j, k);
                            vmac(i, j, k) = vmac(i + 1, j, k);
#if (AMREX_SPACEDIM == 3)
                            wmac(i, j, k) = wmac(i + 1, j, k);
#endif
                            break;
                        case Symmetry:
                            umac(i, j, k) = -umac(i + 2, j, k);
                            vmac(i, j, k) = vmac(i + 1, j, k);
#if (AMREX_SPACEDIM == 3)
                            wmac(i, j, k) = wmac(i + 1, j, k);
#endif
                            break;
                        case Interior:
                            // do nothing
                            break;
                    }
                }

                // hi x-faces
                if (i == domhi[0] + 2) {
                    switch (physbc_p[AMREX_SPACEDIM]) {
                        case Inflow:
                            umac(i, j, k) = umac(i - 1, j, k);
                            vmac(i - 1, j, k) = 0.0;
#if (AMREX_SPACEDIM == 3)
                            wmac(i - 1, j, k) = 0.0;
#endif
                            break;
                        case SlipWall:
                        case NoSlipWall:
                            umac(i, j, k) = 0.0;
                            vmac(i - 1, j, k) = 0.0;
#if (AMREX_SPACEDIM == 3)
                            wmac(i - 1, j, k) = 0.0;
#endif
                            break;
                        case Outflow:
                            umac(i, j, k) = umac(i - 1, j, k);
                            vmac(i - 1, j, k) = vmac(i - 2, j, k);
#if (AMREX_SPACEDIM == 3)
                            wmac(i - 1, j, k) = wmac(i - 2, j, k);
#endif
                            break;
                        case Symmetry:
                            umac(i, j, k) = -umac(i - 2, j, k);
                            vmac(i - 1, j, k) = vmac(i - 2, j, k);
#if (AMREX_SPACEDIM == 3)
                            wmac(i - 1, j, k) = wmac(i - 2, j, k);
#endif
                            break;
                        case Interior:
                            // do nothing
                            break;
                    }
                }
            });

            Gpu::synchronize();

            ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                // lo y-faces
                if (j == domlo[1] - 1) {
                    switch (physbc_p[1]) {
                        case Inflow:
                            umac(i, j, k) = 0.0;
                            vmac(i, j, k) = vmac(i, j + 1, k);
#if (AMREX_SPACEDIM == 3)
                            wmac(i, j, k) = 0.0;
#endif
                            break;
                        case SlipWall:
                        case NoSlipWall:
                            umac(i, j, k) = 0.0;
                            vmac(i, j, k) = 0.0;
#if (AMREX_SPACEDIM == 3)
                            wmac(i, j, k) = 0.0;
#endif
                            break;
                        case Outflow:
                            umac(i, j, k) = umac(i, j + 1, k);
                            vmac(i, j, k) = vmac(i, j + 1, k);
#if (AMREX_SPACEDIM == 3)
                            wmac(i, j, k) = wmac(i, j + 1, k);
#endif
                            break;
                        case Symmetry:
                            umac(i, j, k) = umac(i, j + 1, k);
                            vmac(i, j, k) = -vmac(i, j + 2, k);
#if (AMREX_SPACEDIM == 3)
                            wmac(i, j, k) = wmac(i, j + 1, k);
#endif
                            break;
                        case Interior:
                            // do nothing
                            break;
                    }
                }

                // hi y-faces
                if (j == domhi[1] + 2) {
                    switch (physbc_p[AMREX_SPACEDIM + 1]) {
                        case Inflow:
                            umac(i, j - 1, k) = 0.0;
                            vmac(i, j, k) = vmac(i, j - 1, k);
#if (AMREX_SPACEDIM == 3)
                            wmac(i, j - 1, k) = 0.0;
#endif
                            break;
                        case SlipWall:
                        case NoSlipWall:
                            umac(i, j - 1, k) = 0.0;
                            vmac(i, j, k) = 0.0;
#if (AMREX_SPACEDIM == 3)
                            wmac(i, j - 1, k) = 0.0;
#endif
                            break;
                        case Outflow:
                            umac(i, j - 1, k) = umac(i, j - 2, k);
                            vmac(i, j, k) = vmac(i, j - 1, k);
#if (AMREX_SPACEDIM == 3)
                            wmac(i, j - 1, k) = wmac(i, j - 2, k);
#endif
                            break;
                        case Symmetry:
                            umac(i, j - 1, k) = umac(i, j - 2, k);
                            vmac(i, j, k) = -vmac(i, j - 2, k);
#if (AMREX_SPACEDIM == 3)
                            wmac(i, j - 1, k) = wmac(i, j - 2, k);
#endif
                            break;
                        case Interior:
                            // do nothing
                            break;
                    }
                }
            });

#if (AMREX_SPACEDIM == 3)

            Gpu::synchronize();

            ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                // lo z-faces
                if (k == domlo[2] - 1) {
                    switch (physbc_p[2]) {
                        case Inflow:
                            umac(i, j, k) = 0.0;
                            vmac(i, j, k) = 0.0;
                            wmac(i, j, k) = wmac(i, j, k + 1);
                            break;
                        case SlipWall:
                        case NoSlipWall:
                            umac(i, j, k) = 0.0;
                            vmac(i, j, k) = 0.0;
                            wmac(i, j, k) = 0.0;
                            break;
                        case Outflow:
                            umac(i, j, k) = umac(i, j, k + 1);
                            vmac(i, j, k) = vmac(i, j, k + 1);
                            wmac(i, j, k) = wmac(i, j, k + 1);
                            break;
                        case Symmetry:
                            umac(i, j, k) = umac(i, j, k + 1);
                            vmac(i, j, k) = vmac(i, j, k + 1);
                            wmac(i, j, k) = -wmac(i, j, k + 2);
                            break;
                        case Interior:
                            // do nothing
                            break;
                    }
                }

                // hi z-faces
                if (k == domhi[2] + 2) {
                    switch (physbc_p[2 + AMREX_SPACEDIM]) {
                        case Inflow:
                            umac(i, j, k - 1) = 0.0;
                            vmac(i, j, k - 1) = 0.0;
                            wmac(i, j, k) = wmac(i, j, k - 1);
                            break;
                        case SlipWall:
                        case NoSlipWall:
                            umac(i, j, k - 1) = 0.0;
                            vmac(i, j, k - 1) = 0.0;
                            wmac(i, j, k) = 0.0;
                            break;
                        case Outflow:
                            umac(i, j, k - 1) = umac(i, j, k - 2);
                            vmac(i, j, k - 1) = vmac(i, j, k - 2);
                            wmac(i, j, k) = wmac(i, j, k - 1);
                            break;
                        case Symmetry:
                            umac(i, j, k - 1) = umac(i, j, k - 2);
                            vmac(i, j, k - 1) = vmac(i, j, k - 2);
                            wmac(i, j, k) = -wmac(i, j, k - 2);
                            break;
                        case Interior:
                            // do nothing
                            break;
                    }
                }
            });
#endif
        }
    }
}

// fill in all ghost cells for an edge-based MAC velocity field
void Maestro::FillPatchUedge(
    Vector<std::array<MultiFab, AMREX_SPACEDIM> >& uedge) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::FillPatchUedge()", FillPatchUedge);

    // in IAMR and original MAESTRO this routine was called "create_umac_grown"

    int nGrow = uedge[0][0].nGrow();

    for (int lev = 0; lev <= finest_level; ++lev) {
        // for refined levels we need to "fillpatch" the MAC velocity field
        if (lev > 0) {
            // create a BoxArray whose boxes include only the cell-centered ghost cells
            // associated with the fine grids
            BoxList f_bndry_bl = amrex::GetBndryCells(grids[lev], nGrow);
            BoxArray f_bndry_ba(std::move(f_bndry_bl));
            f_bndry_ba.maxSize(32);

            // create a coarsened version of the fine ghost cell BoxArray
            BoxArray c_bndry_ba = f_bndry_ba;
            c_bndry_ba.coarsen(refRatio(lev - 1));

            // recreate the fine ghost cell BoxArray so it overlaps perfectly
            // with the coarse ghost cell BoxArray
            // (if there was an odd number of fine ghost cells previously this
            // makes the coarse and fine BoxArrays cover the same physical locations)
            f_bndry_ba = c_bndry_ba;
            f_bndry_ba.refine(refRatio(lev - 1));

            for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
                // crse_src & fine_src must have same parallel distribution.
                // We'll use the KnapSack distribution for the fine_src_ba.
                // Since fine_src_ba should contain more points, this'll lead
                // to a better distribution.
                BoxArray crse_src_ba(c_bndry_ba);
                BoxArray fine_src_ba(f_bndry_ba);

                // grow BoxArrays nodally
                crse_src_ba.surroundingNodes(dir);
                fine_src_ba.surroundingNodes(dir);

                // number of boxes and weights used for KnapSack distribution
                const int N = fine_src_ba.size();
                std::vector<long> wgts(N);

#ifdef _OPENMP
#pragma omp parallel for
#endif
                // set weights equal to number of points in the box
                for (int i = 0; i < N; ++i) {
                    wgts[i] = fine_src_ba[i].numPts();
                }

                // This DM won't be put into the cache.
                DistributionMapping dm;
                dm.KnapSackProcessorMap(wgts, ParallelDescriptor::NProcs());

                // allocate coarse and fine umac in the boundary region
                MultiFab crse_src;
                MultiFab fine_src;

                crse_src.define(crse_src_ba, dm, 1, 0);
                fine_src.define(fine_src_ba, dm, 1, 0);

                crse_src.setVal(1.e200);
                fine_src.setVal(1.e200);

                // We want to fill crse_src from coarse uedge
                crse_src.ParallelCopy(uedge[lev - 1][dir], 0, 0, 1, 1, 0);

#ifdef _OPENMP
#pragma omp parallel
#endif
                for (MFIter mfi(crse_src); mfi.isValid(); ++mfi) {
                    const int nComp = 1;
                    const Box& box = crse_src[mfi].box();
                    IntVect ratio = refRatio(lev - 1);

                    const auto lo = box.loVect3d();
                    const auto hi = box.hiVect3d();

                    const Array4<Real> fine = fine_src.array(mfi);
                    const Array4<Real> crse = crse_src.array(mfi);

                    // For edge-based data, fill fine values with piecewise-constant interp of coarse data.
                    // Operate only on faces that overlap--ie, only fill the fine faces that make up each
                    // coarse face, leave the in-between faces alone.
#if (AMREX_SPACEDIM == 2)
                    int k = lo[2];
                    if (dir == 0) {
                        for (auto n = 0; n < nComp; ++n) {
                            for (auto j = lo[1]; j <= hi[1]; ++j) {
                                int jj = ratio[1] * j;
                                for (auto i = lo[0]; i <= hi[0]; ++i) {
                                    int ii = ratio[0] * i;
                                    for (auto L = 0; L < ratio[1]; ++L) {
                                        fine(ii, jj + L, k, n) =
                                            crse(i, j, k, n);
                                    }
                                }
                            }
                        }
                    } else {
                        for (auto n = 0; n < nComp; ++n) {
                            for (auto j = lo[1]; j <= hi[1]; ++j) {
                                int jj = ratio[1] * j;
                                for (auto i = lo[0]; i <= hi[0]; ++i) {
                                    int ii = ratio[0] * i;
                                    for (auto L = 0; L < ratio[0]; ++L) {
                                        fine(ii + L, jj, k, n) =
                                            crse(i, j, k, n);
                                    }
                                }
                            }
                        }
                    }
#elif (AMREX_SPACEDIM == 3)
                    if (dir == 0) {
                        for (auto n = 0; n < nComp; ++n) {
                            for (auto k = lo[2]; k <= hi[2]; ++k) {
                                int kk = ratio[2] * k;
                                for (auto j = lo[1]; j <= hi[1]; ++j) {
                                    int jj = ratio[1] * j;
                                    for (auto i = lo[0]; i <= hi[0]; ++i) {
                                        int ii = ratio[0] * i;
                                        for (auto P = 0; P < ratio[2]; ++P) {
                                            for (auto L = 0; L < ratio[1];
                                                 ++L) {
                                                fine(ii, jj + L, kk + P, n) =
                                                    crse(i, j, k, n);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    } else if (dir == 1) {
                        for (auto n = 0; n < nComp; ++n) {
                            for (auto k = lo[2]; k <= hi[2]; ++k) {
                                int kk = ratio[2] * k;
                                for (auto j = lo[1]; j <= hi[1]; ++j) {
                                    int jj = ratio[1] * j;
                                    for (auto i = lo[0]; i <= hi[0]; ++i) {
                                        int ii = ratio[0] * i;
                                        for (auto P = 0; P < ratio[2]; ++P) {
                                            for (auto L = 0; L < ratio[0];
                                                 ++L) {
                                                fine(ii + L, jj, kk + P, n) =
                                                    crse(i, j, k, n);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    } else {
                        for (auto n = 0; n < nComp; ++n) {
                            for (auto k = lo[2]; k <= hi[2]; ++k) {
                                int kk = ratio[2] * k;
                                for (auto j = lo[1]; j <= hi[1]; ++j) {
                                    int jj = ratio[1] * j;
                                    for (auto i = lo[0]; i <= hi[0]; ++i) {
                                        int ii = ratio[0] * i;
                                        for (auto P = 0; P < ratio[1]; ++P) {
                                            for (auto L = 0; L < ratio[0];
                                                 ++L) {
                                                fine(ii + L, jj + P, kk, n) =
                                                    crse(i, j, k, n);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
#endif
                }
                crse_src.clear();

                // Replace pc-interpd fine data with preferred u_mac data at
                // this level u_mac valid only on surrounding faces of valid
                // region - this op will not fill grow region.
                fine_src.ParallelCopy(uedge[lev][dir]);

                // Interpolate unfilled grow cells using best data from
                // surrounding faces of valid region, and pc-interpd data
                // on fine edges overlaying coarse edges.
#ifdef _OPENMP
#pragma omp parallel
#endif
                for (MFIter mfi(fine_src); mfi.isValid(); ++mfi) {
                    const int nComp = 1;
                    const Box& fbox = fine_src[mfi].box();
                    IntVect ratio = refRatio(lev - 1);

                    const auto flo = fbox.loVect3d();
                    const auto fhi = fbox.hiVect3d();

                    const Array4<Real> fine = fine_src.array(mfi);

                    // Do linear in dir, pc transverse to dir, leave alone the fine values
                    // lining up with coarse edges--assume these have been set to hold the
                    // values you want to interpolate to the rest.

#if (AMREX_SPACEDIM == 2)
                    int k = flo[2];
                    if (dir == 0) {
                        for (auto n = 0; n < nComp; ++n) {
                            for (auto j = flo[1]; j <= fhi[1]; j += ratio[1]) {
                                for (auto i = flo[0]; i <= fhi[0] - ratio[dir];
                                     i += ratio[0]) {
                                    Real df = fine(i + ratio[dir], j, k, n) -
                                              fine(i, j, k, n);
                                    for (auto M = 1; M < ratio[dir]; ++M) {
                                        Real val =
                                            fine(i, j, k, n) +
                                            df * Real(M) / Real(ratio[dir]);
                                        for (auto P = amrex::max(j, flo[1]);
                                             P <= amrex::min(j + ratio[1] - 1,
                                                             fhi[1]);
                                             ++P) {
                                            fine(i + M, P, k, n) = val;
                                        }
                                    }
                                }
                            }
                        }
                    } else {
                        for (auto n = 0; n < nComp; ++n) {
                            for (auto j = flo[1]; j <= fhi[1] - ratio[dir];
                                 j += ratio[1]) {
                                for (auto i = flo[0]; i <= fhi[0];
                                     i += ratio[0]) {
                                    Real df = fine(i, j + ratio[dir], k, n) -
                                              fine(i, j, k, n);
                                    for (auto M = 1; M < ratio[dir]; ++M) {
                                        Real val =
                                            fine(i, j, k, n) +
                                            df * Real(M) / Real(ratio[dir]);
                                        for (auto P = amrex::max(i, flo[0]);
                                             P <= amrex::min(i + ratio[0] - 1,
                                                             fhi[0]);
                                             ++P) {
                                            fine(P, j + M, k, n) = val;
                                        }
                                    }
                                }
                            }
                        }
                    }
#elif (AMREX_SPACEDIM == 3)
                    if (dir == 0) {
                        for (auto n = 0; n < nComp; ++n) {
                            for (auto k = flo[2]; k <= fhi[2]; k += ratio[2]) {
                                for (auto j = flo[1]; j <= fhi[1];
                                     j += ratio[1]) {
                                    for (auto i = flo[0];
                                         i <= fhi[0] - ratio[dir];
                                         i += ratio[0]) {
                                        Real df =
                                            fine(i + ratio[dir], j, k, n) -
                                            fine(i, j, k, n);
                                        for (auto M = 1; M < ratio[dir]; ++M) {
                                            Real val =
                                                fine(i, j, k, n) +
                                                df * Real(M) / Real(ratio[dir]);
                                            for (auto P = amrex::max(j, flo[1]);
                                                 P <=
                                                 amrex::min(j + ratio[1] - 1,
                                                            fhi[1]);
                                                 ++P) {
                                                for (auto L =
                                                         amrex::max(k, flo[2]);
                                                     L <= amrex::min(
                                                              k + ratio[2] - 1,
                                                              fhi[2]);
                                                     ++L) {
                                                    fine(i + M, P, L, n) = val;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    } else if (dir == 1) {
                        for (auto n = 0; n < nComp; ++n) {
                            for (auto k = flo[2]; k <= fhi[2]; k += ratio[2]) {
                                for (auto j = flo[1]; j <= fhi[1] - ratio[dir];
                                     j += ratio[1]) {
                                    for (auto i = flo[0]; i <= fhi[0];
                                         i += ratio[0]) {
                                        Real df =
                                            fine(i, j + ratio[dir], k, n) -
                                            fine(i, j, k, n);
                                        for (auto M = 1; M < ratio[dir]; ++M) {
                                            Real val =
                                                fine(i, j, k, n) +
                                                df * Real(M) / Real(ratio[dir]);
                                            for (auto P = amrex::max(i, flo[0]);
                                                 P <=
                                                 amrex::min(i + ratio[0] - 1,
                                                            fhi[0]);
                                                 ++P) {
                                                for (auto L =
                                                         amrex::max(k, flo[2]);
                                                     L <= amrex::min(
                                                              k + ratio[2] - 1,
                                                              fhi[2]);
                                                     ++L) {
                                                    fine(P, j + M, L, n) = val;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    } else {
                        for (auto n = 0; n < nComp; ++n) {
                            for (auto k = flo[2]; k <= fhi[2] - ratio[dir];
                                 k += ratio[2]) {
                                for (auto j = flo[1]; j <= fhi[1];
                                     j += ratio[1]) {
                                    for (auto i = flo[0]; i <= fhi[0];
                                         i += ratio[0]) {
                                        Real df =
                                            fine(i, j, k + ratio[dir], n) -
                                            fine(i, j, k, n);
                                        for (auto M = 1; M < ratio[dir]; ++M) {
                                            Real val =
                                                fine(i, j, k, n) +
                                                df * Real(M) / Real(ratio[dir]);
                                            for (auto P = amrex::max(i, flo[0]);
                                                 P <=
                                                 amrex::min(i + ratio[0] - 1,
                                                            fhi[0]);
                                                 ++P) {
                                                for (auto L =
                                                         amrex::max(j, flo[1]);
                                                     L <= amrex::min(
                                                              j + ratio[1] - 1,
                                                              fhi[1]);
                                                     ++L) {
                                                    fine(P, L, k + M, n) = val;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
#endif
                }

                // make a copy of the original fine uedge but with no ghost cells
                MultiFab uedge_f_save(uedge[lev][dir].boxArray(),
                                      uedge[lev][dir].DistributionMap(), 1, 0);
                uedge_f_save.ParallelCopy(uedge[lev][dir]);

                // copy in the grown data into fine uedge
                uedge[lev][dir].ParallelCopy(fine_src, 0, 0, 1, 0, nGrow);

                // copy the original valid region back into uedge
                // to we don't change the values on the C-F interface
                uedge[lev][dir].ParallelCopy(uedge_f_save);
            }

        }  // end if

        // fill periodic ghost cells
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            uedge[lev][d].FillBoundary(geom[lev].periodicity());
        }

        // fill ghost cells behind physical boundaries
        FillUmacGhost(uedge, lev);

    }  // end loop over levels
}

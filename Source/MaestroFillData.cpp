
#include <Maestro.H>

using namespace amrex;

// call FillPatch for all levels
void
Maestro::FillPatch (Real time,
                    Vector<MultiFab>& mf,
                    Vector<MultiFab>& mf_old,
                    Vector<MultiFab>& mf_new,
                    int srccomp, int destcomp, int ncomp, int startbccomp,
                    const Vector<BCRec>& bcs_in)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::FillPatch()",FillPatch);

    for (int lev=0; lev<=finest_level; ++lev) {
        FillPatch(lev, time, mf[lev], mf_old, mf_new, srccomp, destcomp, ncomp,
                  startbccomp, bcs_in);
    }
}

// compute a new multifab by coping in phi from valid region and filling ghost cells
// works for single level and 2-level cases
// (fill fine grid ghost by interpolating from coarse)
// srccomp of the source component
// destcomp is the destination component AND the bc component
void
Maestro::FillPatch (int lev, Real time, MultiFab& mf,
                    Vector<MultiFab>& mf_old,
                    Vector<MultiFab>& mf_new,
                    int srccomp, int destcomp, int ncomp, int startbccomp,
                    const Vector<BCRec>& bcs_in)
{

    Vector<BCRec> bcs{bcs_in.begin()+startbccomp,bcs_in.begin()+startbccomp+ncomp};

    if (lev == 0)
    {
        Vector<MultiFab*> smf;
        Vector<Real> stime;
        GetData(0, time, smf, stime, mf_old, mf_new);

        PhysBCFunctMaestro physbc(geom[lev],bcs,BndryFuncArray(phifill));
        FillPatchSingleLevel(mf, time, smf, stime, srccomp, destcomp, ncomp,
                             geom[lev], physbc, 0);
    }
    else
    {
        Vector<MultiFab*> cmf, fmf;
        Vector<Real> ctime, ftime;
        GetData(lev-1, time, cmf, ctime, mf_old, mf_new);
        GetData(lev, time, fmf, ftime, mf_old, mf_new);

        PhysBCFunctMaestro cphysbc(geom[lev-1],bcs,BndryFuncArray(phifill));
        PhysBCFunctMaestro fphysbc(geom[lev  ],bcs,BndryFuncArray(phifill));

        Interpolater* mapper = &cell_cons_interp;

        FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                           srccomp, destcomp, ncomp, geom[lev-1], geom[lev],
                           cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                           mapper, bcs, 0);
    }
}

// fill an entire multifab by interpolating from the coarser level
// this comes into play when a new level of refinement appears
// srccomp of the source component
// destcomp is the destination component AND the bc component
void
Maestro::FillCoarsePatch (int lev, Real time, MultiFab& mf,
                          Vector<MultiFab>& mf_old,
                          Vector<MultiFab>& mf_new,
                          int srccomp, int destcomp, int ncomp,
                          const Vector<BCRec>& bcs, bool is_vel)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::FillCoarsePatch()",FillCoarsePatch);

    AMREX_ASSERT(lev > 0);

    Vector<MultiFab*> cmf;
    Vector<Real> ctime;
    GetData(lev-1, time, cmf, ctime, mf_old, mf_new);

    if (cmf.size() != 1) {
        Abort("FillCoarsePatch: how did this happen?");
    }

    Interpolater* mapper = &cell_cons_interp;

    if (is_vel){
        PhysBCFunctMaestro cphysbc(geom[lev-1],bcs,BndryFuncArray(velfill));
        PhysBCFunctMaestro fphysbc(geom[lev  ],bcs,BndryFuncArray(velfill));
        
        InterpFromCoarseLevel(mf, time, *cmf[0], srccomp, destcomp, ncomp, geom[lev-1], geom[lev],
                              cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                              mapper, bcs, 0);
    } else {
        PhysBCFunctMaestro cphysbc(geom[lev-1],bcs,BndryFuncArray(phifill));
        PhysBCFunctMaestro fphysbc(geom[lev  ],bcs,BndryFuncArray(phifill));

        InterpFromCoarseLevel(mf, time, *cmf[0], srccomp, destcomp, ncomp, geom[lev-1], geom[lev],
                              cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                              mapper, bcs, 0);
    }
}

// utility to copy in data from mf_old and/or mf_new into mf
// if time=t_old we copy mf_old into mf
// if time=t_new we copy mf_new into mf
// otherwise copy copy in both mf_old and mf_new into mf and the fillpatch
// routines know to interpolate in time.  However in MAESTRO since we don't
// subcycle I'm not sure if we need this capability?
void
Maestro::GetData (int lev, Real time,
                  Vector<MultiFab*>& mf,
                  Vector<Real>& mftime,
                  Vector<MultiFab>& mf_old,
                  Vector<MultiFab>& mf_new)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::GetData()",GetData);

    mf.clear();
    mftime.clear();

    const Real teps = (t_new - t_old) * 1.e-3;

    if (time >= t_new - teps && time <= t_new + teps)
    {
        mf.push_back(&mf_new[lev]);
        mftime.push_back(t_new);
    }
    else if (time >= t_old - teps && time <= t_old + teps)
    {
        mf.push_back(&mf_old[lev]);
        mftime.push_back(t_old);
    }
    else
    {
        mf.push_back(&mf_old[lev]);
        mf.push_back(&mf_new[lev]);
        mftime.push_back(t_old);
        mftime.push_back(t_new);
    }
}

// set covered coarse cells to be the average of overlying fine cells
void
Maestro::AverageDown (Vector<MultiFab>& mf,
                      int comp,
                      int ncomp)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::AverageDown()",AverageDown);

    for (int lev = finest_level-1; lev >= 0; --lev) {
        average_down(mf[lev+1], mf[lev],
                     geom[lev+1], geom[lev],
                     comp, ncomp, refRatio(lev));
    }
}

// set covered faces to be the average of overlying fine faces
void
Maestro::AverageDownFaces (Vector<std::array< MultiFab, AMREX_SPACEDIM > >& edge)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::AverageDownFaces()",AverageDownFaces);

    for (int lev = finest_level-1; lev >= 0; --lev) {

        Vector<const MultiFab*> edge_f(AMREX_SPACEDIM);
        Vector<      MultiFab*> edge_c(AMREX_SPACEDIM);

        for (int dir=0; dir<AMREX_SPACEDIM; ++dir) {
            edge_f[dir] = &(edge[lev+1][dir]);
            edge_c[dir] = &(edge[lev  ][dir]);
        }

        average_down_faces(edge_f, edge_c, refRatio(lev));
    }

}


// fill in ONE ghost cell for all components of a face-centered (MAC) velocity
// field behind physical boundaries.  Does not modify the velocities on the boundary
void
Maestro::FillUmacGhost (Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
                        int level)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::FillUmacGhost()",FillUmacGhost);

    int start_lev;
    int end_lev;
    if (level == -1) {
        start_lev = 0;
        end_lev = finest_level;
    }
    else {
        start_lev = level;
        end_lev = level;
    }

    for (int lev=start_lev; lev<=end_lev; ++lev) {

        // Get the index space of the domain
        const Box& domainBox = geom[lev].Domain();

        // get references to the MultiFabs at level lev
        MultiFab& sold_mf = sold[lev];         // need a cell-centered MF for the MFIter
        MultiFab& umacx_mf = umac[lev][0];
        MultiFab& umacy_mf = umac[lev][1];
#if (AMREX_SPACEDIM == 3)
        MultiFab& umacz_mf = umac[lev][2];
#endif

        // DO NOT TILE THIS SUBROUTINE
        // this just filling ghost cells so the fortran logic has to be reworked
        // to properly capture the corner terms
        for ( MFIter mfi(sold_mf, false); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid (cell-centered) region
            const Box& tilebox = mfi.tilebox();

            fill_umac_ghost(ARLIM_3D(domainBox.loVect()), ARLIM_3D(domainBox.hiVect()),
                            ARLIM_3D(tilebox.loVect()), ARLIM_3D(tilebox.hiVect()),
                            BL_TO_FORTRAN_3D(umacx_mf[mfi]),
                            BL_TO_FORTRAN_3D(umacy_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
                            BL_TO_FORTRAN_3D(umacz_mf[mfi]),
#endif
                            phys_bc.dataPtr());

        }
    }
}

// fill in all ghost cells for an edge-based MAC velocity field
void
Maestro::FillPatchUedge (Vector<std::array< MultiFab, AMREX_SPACEDIM > >& uedge)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::FillPatchUedge()",FillPatchUedge);

    // in IAMR and original MAESTRO this routine was called "create_umac_grown"

    int nGrow = uedge[0][0].nGrow();

    for (int lev=0; lev<= finest_level; ++lev) {

        // for refined levels we need to "fillpatch" the MAC velocity field
        if (lev > 0) {

            // create a BoxArray whose boxes include only the cell-centered ghost cells
            // associated with the fine grids
            BoxList f_bndry_bl = amrex::GetBndryCells(grids[lev],nGrow);
            BoxArray f_bndry_ba(std::move(f_bndry_bl));
            f_bndry_ba.maxSize(32);

            // create a coarsened version of the fine ghost cell BoxArray
            BoxArray c_bndry_ba = f_bndry_ba;
            c_bndry_ba.coarsen(refRatio(lev-1));

            // recreate the fine ghost cell BoxArray so it overlaps perfectly
            // with the coarse ghost cell BoxArray
            // (if there was an odd number of fine ghost cells previously this
            // makes the coarse and fine BoxArrays cover the same physical locations)
            f_bndry_ba = c_bndry_ba;
            f_bndry_ba.refine(refRatio(lev-1));

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
                dm.KnapSackProcessorMap(wgts,ParallelDescriptor::NProcs());

                // allocate coarse and fine umac in the boundary region
                MultiFab crse_src;
                MultiFab fine_src;

                crse_src.define(crse_src_ba, dm, 1, 0);
                fine_src.define(fine_src_ba, dm, 1, 0);

                crse_src.setVal(1.e200);
                fine_src.setVal(1.e200);

                // We want to fill crse_src from coarse uedge
                crse_src.copy(uedge[lev-1][dir],0,0,1,1,0);

#ifdef _OPENMP
#pragma omp parallel
#endif
                for (MFIter mfi(crse_src); mfi.isValid(); ++mfi) {
                    const int nComp = 1;
                    const Box& box   = crse_src[mfi].box();
                    IntVect rr = refRatio(lev-1);
                    const int* rat = rr.getVect();
                    // For edge-based data, fill fine values with piecewise-constant interp of coarse data.
                    // Operate only on faces that overlap--ie, only fill the fine faces that make up each
                    // coarse face, leave the in-between faces alone.
                    PC_EDGE_INTERP(box.loVect(), box.hiVect(), &nComp, rat, &dir,
                                   BL_TO_FORTRAN_FAB(crse_src[mfi]),
                                   BL_TO_FORTRAN_FAB(fine_src[mfi]));
                }
                crse_src.clear();
                //
                // Replace pc-interpd fine data with preferred u_mac data at
                // this level u_mac valid only on surrounding faces of valid
                // region - this op will not fill grow region.
                //
                fine_src.copy(uedge[lev][dir]);
                //
                // Interpolate unfilled grow cells using best data from
                // surrounding faces of valid region, and pc-interpd data
                // on fine edges overlaying coarse edges.
                //
#ifdef _OPENMP
#pragma omp parallel
#endif
                for (MFIter mfi(fine_src); mfi.isValid(); ++mfi) {
                    const int nComp = 1;
                    const Box& fbox  = fine_src[mfi].box();
                    IntVect rr = refRatio(lev-1);
                    const int* rat = rr.getVect();
                    // Do linear in dir, pc transverse to dir, leave alone the fine values
                    // lining up with coarse edges--assume these have been set to hold the
                    // values you want to interpolate to the rest.
                    EDGE_INTERP(fbox.loVect(), fbox.hiVect(), &nComp, rat, &dir,
                                BL_TO_FORTRAN_FAB(fine_src[mfi]));
                }

                // make a copy of the original fine uedge but with no ghost cells
                MultiFab uedge_f_save(uedge[lev][dir].boxArray(),uedge[lev][dir].DistributionMap(), 1,0);
                uedge_f_save.copy(uedge[lev][dir]);

                // copy in the grown data into fine uedge
                uedge[lev][dir].copy(fine_src,0,0,1,0,nGrow);

                // copy the original valid region back into uedge
                // to we don't change the values on the C-F interface
                uedge[lev][dir].copy(uedge_f_save);
            }

        }         // end if

        // fill periodic ghost cells
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            uedge[lev][d].FillBoundary(geom[lev].periodicity());
        }

        // fill ghost cells behind physical boundaries
        FillUmacGhost(uedge,lev);

    }     // end loop over levels

}

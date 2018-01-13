
#include <Maestro.H>

using namespace amrex;

// compute a new multifab by coping in phi from valid region and filling ghost cells
// works for single level and 2-level cases
// (fill fine grid ghost by interpolating from coarse)
// scomp of the source component
// dcomp is the destination component AND the bc component
void
Maestro::FillPatch (int lev, Real time, MultiFab& mf, 
                    Vector<MultiFab>& mf_old,
                    Vector<MultiFab>& mf_new,
                    int scomp, int dcomp, int ncomp,
                    Vector<BCRec> bcs)
{
    if (lev == 0)
    {
        Vector<MultiFab*> smf;
        Vector<Real> stime;
        GetData(0, time, smf, stime, mf_old, mf_new);

        PhysBCFunctMaestro physbc(geom[lev],bcs,BndryFunctBase(phifill));
        FillPatchSingleLevel(mf, time, smf, stime, scomp, dcomp, ncomp,
                             geom[lev], physbc);
    }
    else
    {
        Vector<MultiFab*> cmf, fmf;
        Vector<Real> ctime, ftime;
        GetData(lev-1, time, cmf, ctime, mf_old, mf_new);
        GetData(lev  , time, fmf, ftime, mf_old, mf_new);

        PhysBCFunctMaestro cphysbc(geom[lev-1],bcs,BndryFunctBase(phifill));
        PhysBCFunctMaestro fphysbc(geom[lev  ],bcs,BndryFunctBase(phifill));

        Interpolater* mapper = &cell_cons_interp;

        FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                           scomp, dcomp, ncomp, geom[lev-1], geom[lev],
                           cphysbc, fphysbc, refRatio(lev-1),
                           mapper, bcs);
    }
}

// fill an entire multifab by interpolating from the coarser level
// this comes into play when a new level of refinement appears
// scomp of the source component
// dcomp is the destination component AND the bc component
void
Maestro::FillCoarsePatch (int lev, Real time, MultiFab& mf,
                          Vector<MultiFab>& mf_old,
                          Vector<MultiFab>& mf_new,
                          int scomp, int dcomp, int ncomp,
                          Vector<BCRec> bcs)
{
    AMREX_ASSERT(lev > 0);

    Vector<MultiFab*> cmf;
    Vector<Real> ctime;
    GetData(lev-1, time, cmf, ctime, mf_old, mf_new);
    
    if (cmf.size() != 1) {
        Abort("FillCoarsePatch: how did this happen?");
    }

    PhysBCFunctMaestro cphysbc(geom[lev-1],bcs,BndryFunctBase(phifill));
    PhysBCFunctMaestro fphysbc(geom[lev  ],bcs,BndryFunctBase(phifill));

    Interpolater* mapper = &cell_cons_interp;
    InterpFromCoarseLevel(mf, time, *cmf[0], scomp, dcomp, ncomp, geom[lev-1], geom[lev],
                          cphysbc, fphysbc, refRatio(lev-1),
                          mapper, bcs);
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
    for (int lev = finest_level-1; lev >= 0; --lev)
    {
        average_down(mf[lev+1], mf[lev],
                     geom[lev+1], geom[lev],
                     comp, ncomp, refRatio(lev));
    }
}

// more flexible version of AverageDown() that lets you average down across multiple levels
void
Maestro::AverageDownTo (int crse_lev,
                        Vector<MultiFab>& mf,
                        int comp,
                        int ncomp)
{
    average_down(mf[crse_lev+1], mf[crse_lev],
                 geom[crse_lev+1], geom[crse_lev],
                 comp, ncomp, refRatio(crse_lev));
}


// fill in ONE ghost cell for all components of a face-centered (MAC) velocity
// field behind physical boundaries.  Does not modify the velocities on the boundary
void
Maestro::FillUmacGhost (Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac)
{
    for (int lev=0; lev<=finest_level; ++lev) {
        MultiFab& sold_mf = sold[lev];  // need a cell-centered MF for the MFIter        
        MultiFab& umacx_mf = umac[lev][0];
        MultiFab& umacy_mf = umac[lev][1];
#if (AMREX_SPACEDIM == 3)
        MultiFab& umacz_mf = umac[lev][2];
#endif

        for ( MFIter mfi(sold_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid (cell-centered) region
            const Box& validBox = mfi.validbox();
            const Box& domainBox = geom[lev].Domain();

            fill_umac_ghost(ARLIM_3D(domainBox.loVect()), ARLIM_3D(domainBox.hiVect()),
                            ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                            BL_TO_FORTRAN_3D(umacx_mf[mfi]),
                            BL_TO_FORTRAN_3D(umacy_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
                            BL_TO_FORTRAN_3D(umacz_mf[mfi]),
#endif
                            lo_bc.dataPtr(),hi_bc.data());

        }
    }
}

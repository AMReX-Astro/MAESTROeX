
#include <Maestro.H>

using namespace amrex;

// compute a new multifab by coping in phi from valid region and filling ghost cells
// works for single level and 2-level cases
// (fill fine grid ghost by interpolating from coarse)
void
Maestro::FillPatch (int lev, Real time, MultiFab& mf, 
                    int icomp, int ncomp, Vector<BCRec> bcs)
{
    if (lev == 0)
    {
        Vector<MultiFab*> smf;
        Vector<Real> stime;
        GetData(0, time, smf, stime);

        MaestroPhysBCFunct physbc(geom[lev],bcs,BndryFunctBase(phifill));
        FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
                             geom[lev], physbc);
    }
    else
    {
        Vector<MultiFab*> cmf, fmf;
        Vector<Real> ctime, ftime;
        GetData(lev-1, time, cmf, ctime);
        GetData(lev  , time, fmf, ftime);

        MaestroPhysBCFunct cphysbc(geom[lev-1],bcs,BndryFunctBase(phifill));
        MaestroPhysBCFunct fphysbc(geom[lev  ],bcs,BndryFunctBase(phifill));

        Interpolater* mapper = &cell_cons_interp;

        FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                           0, icomp, ncomp, geom[lev-1], geom[lev],
                           cphysbc, fphysbc, refRatio(lev-1),
                           mapper, bcs);
    }
}

// fill an entire multifab by interpolating from the coarser level
// this comes into play when a new level of refinement appears
void
Maestro::FillCoarsePatch (int lev, Real time, MultiFab& mf, 
                          int icomp, int ncomp, Vector<BCRec> bcs)
{
    BL_ASSERT(lev > 0);

    Vector<MultiFab*> cmf;
    Vector<Real> ctime;
    GetData(lev-1, time, cmf, ctime);
    
    if (cmf.size() != 1) {
        Abort("FillCoarsePatch: how did this happen?");
    }

    MaestroPhysBCFunct cphysbc(geom[lev-1],bcs,BndryFunctBase(phifill));
    MaestroPhysBCFunct fphysbc(geom[lev  ],bcs,BndryFunctBase(phifill));

    Interpolater* mapper = &cell_cons_interp;
    InterpFromCoarseLevel(mf, time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],
                          cphysbc, fphysbc, refRatio(lev-1),
                          mapper, bcs);
}

// utility to copy in data from sold and/or snew into another multifab
void
Maestro::GetData (int lev, Real time, Vector<MultiFab*>& data, Vector<Real>& datatime)
{
    data.clear();
    datatime.clear();

    const Real teps = (t_new - t_old) * 1.e-3;

    if (time > t_new - teps && time < t_new + teps)
    {
        data.push_back(snew[lev].get());
        datatime.push_back(t_new);
    }
    else if (time > t_old - teps && time < t_old + teps)
    {
        data.push_back(sold[lev].get());
        datatime.push_back(t_old);
    }
    else
    {
        data.push_back(sold[lev].get());
        data.push_back(snew[lev].get());
        datatime.push_back(t_old);
        datatime.push_back(t_new);
    }
}

// set covered coarse cells to be the average of overlying fine cells
void
Maestro::AverageDown (Vector<std::unique_ptr<MultiFab> >& mf)
{
    for (int lev = finest_level-1; lev >= 0; --lev)
    {
        average_down(*mf[lev+1], *mf[lev],
                     geom[lev+1], geom[lev],
                     0, mf[lev]->nComp(), refRatio(lev));
    }
}

// more flexible version of AverageDown() that lets you average down across multiple levels
void
Maestro::AverageDownTo (int crse_lev, Vector<std::unique_ptr<MultiFab> >& mf)
{
    average_down(*mf[crse_lev+1], *mf[crse_lev],
                 geom[crse_lev+1], geom[crse_lev],
                 0, mf[crse_lev]->nComp(), refRatio(crse_lev));
}

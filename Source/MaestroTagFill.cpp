
#include <Maestro.H>

using namespace amrex;

// tag all cells for refinement
// overrides the pure virtual function in AmrCore
void
Maestro::ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow)
{
    static bool first = true;
    static Vector<Real> phierr;

    // only do this during the first call to ErrorEst
    if (first)
    {
        first = false;
        // read in an array of "phierr", which is the tagging threshold
        // in this example, we tag values of "phi" which are greater than phierr
        // for that particular level
        // in subroutine state_error, you could use more elaborate tagging, such
        // as more advanced logical expressions, or gradients, etc.
        ParmParse pp("maestro");
        int n = pp.countval("phierr");
        if (n > 0) {
            pp.getarr("phierr", phierr, 0, n);
        }
    }

    if (lev >= phierr.size()) return;

    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

    const Real* dx      = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();

    const MultiFab& state = *snew[lev];

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Vector<int>  itags;
	
        for (MFIter mfi(state,true); mfi.isValid(); ++mfi)
        {
            const Box& tilebox  = mfi.tilebox();

            TagBox&     tagfab  = tags[mfi];
	    
            // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
            // So we are going to get a temporary integer array.
            // set itags initially to 'untagged' everywhere
            // we define itags over the tilebox region
            tagfab.get_itags(itags, tilebox);
	    
            // data pointer and index space
            int*        tptr    = itags.dataPtr();
            const int*  tlo     = tilebox.loVect();
            const int*  thi     = tilebox.hiVect();

            // tag cells for refinement
            state_error(tptr,  ARLIM_3D(tlo), ARLIM_3D(thi),
                        BL_TO_FORTRAN_3D(state[mfi]),
                        &tagval, &clearval, 
                        ARLIM_3D(tilebox.loVect()), ARLIM_3D(tilebox.hiVect()), 
                        ZFILL(dx), ZFILL(prob_lo), &time, &phierr[lev]);
            //
            // Now update the tags in the TagBox in the tilebox region
            // to be equal to itags
            //
            tagfab.tags_and_untags(itags, tilebox);
        }
    }
}

// Make a new level using provided BoxArray and DistributionMapping and 
// fill with interpolated coarse level data.
// overrides the pure virtual function in AmrCore
void
Maestro::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
				 const DistributionMapping& dm)
{
    const int nghost = snew[lev-1]->nGrow();
    
    snew[lev].reset(new MultiFab(ba, dm, NSCAL, nghost));
    sold[lev].reset(new MultiFab(ba, dm, NSCAL, nghost));

    t_new = time;
    t_old = time - 1.e200;

    if (lev > 0 && do_reflux) {
        flux_reg_s[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, NSCAL));
        flux_reg_u[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, AMREX_SPACEDIM));
    }

    FillCoarsePatch(lev, time, *snew[lev], 0, NSCAL, bcs_s);
}

// Remake an existing level using provided BoxArray and DistributionMapping and 
// fill with existing fine and coarse data.
// overrides the pure virtual function in AmrCore
void
Maestro::RemakeLevel (int lev, Real time, const BoxArray& ba,
		      const DistributionMapping& dm)
{
    const int nghost = snew[lev]->nGrow();

#if __cplusplus >= 201402L
    auto new_state = std::make_unique<MultiFab>(ba, dm, NSCAL, nghost);
    auto old_state = std::make_unique<MultiFab>(ba, dm, NSCAL, nghost);
#else
    std::unique_ptr<MultiFab> new_state(new MultiFab(ba, dm, NSCAL, nghost));
    std::unique_ptr<MultiFab> old_state(new MultiFab(ba, dm, NSCAL, nghost));
#endif

    FillPatch(lev, time, *new_state, 0, NSCAL, bcs_s);

    std::swap(new_state, snew[lev]);
    std::swap(old_state, sold[lev]);

    t_new = time;
    t_old = time - 1.e200;

    if (lev > 0 && do_reflux) {
        flux_reg_s[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, NSCAL));
        flux_reg_u[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, AMREX_SPACEDIM));
    }    
}



// compute a new multifab by coping in phi from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
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


#include <Maestro.H>

// tag all cells for refinement
// overrides the pure virtual function in AmrCore
void
Maestro::ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow)
{
    static bool first = true;
    static Array<Real> phierr;

    // only do this during the first call to ErrorEst
    if (first)
    {
        first = false;
        // read in an array of "phierr", which is the tagging threshold
        // in this example, we tag values of "phi" which are greater than phierr
        // for that particular level
        // in subroutine state_error, you could use more elaborate tagging, such
        // as more advanced logical expressions, or gradients, etc.
        ParmParse pp("adv");
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

    const MultiFab& state = *phi_new[lev];

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Array<int>  itags;
	
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

// Make a new level from scratch using provided BoxArray and DistributionMapping.
// Only used during initialization.
// overrides the pure virtual function in AmrCore
void Maestro::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
				       const DistributionMapping& dm)
{
    const int ncomp = 2;
    const int nghost = 0;

    phi_new[lev].reset(new MultiFab(ba, dm, ncomp, nghost));
    phi_old[lev].reset(new MultiFab(ba, dm, ncomp, nghost));

    t_new = time;
    t_old = time - 1.e200;

    if (lev > 0 && do_reflux) {
        flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp));
    }

    const Real* dx = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();
    Real cur_time = t_new;

    MultiFab& state = *phi_new[lev];

    for (MFIter mfi(state); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        const int* lo  = box.loVect();
        const int* hi  = box.hiVect();

        initdata(lev, cur_time, ARLIM_3D(lo), ARLIM_3D(hi),
                 BL_TO_FORTRAN_3D(state[mfi]), ZFILL(dx),
                 ZFILL(prob_lo));
    }
}

// Make a new level using provided BoxArray and DistributionMapping and 
// fill with interpolated coarse level data.
// overrides the pure virtual function in AmrCore
void
Maestro::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
				 const DistributionMapping& dm)
{
    const int ncomp = phi_new[lev-1]->nComp();
    const int nghost = phi_new[lev-1]->nGrow();
    
    phi_new[lev].reset(new MultiFab(ba, dm, ncomp, nghost));
    phi_old[lev].reset(new MultiFab(ba, dm, ncomp, nghost));

    t_new = time;
    t_old = time - 1.e200;

    if (lev > 0 && do_reflux) {
        flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp));
    }

    FillCoarsePatch(lev, time, *phi_new[lev], 0, ncomp);
}

// Remake an existing level using provided BoxArray and DistributionMapping and 
// fill with existing fine and coarse data.
// overrides the pure virtual function in AmrCore
void
Maestro::RemakeLevel (int lev, Real time, const BoxArray& ba,
		      const DistributionMapping& dm)
{
    const int ncomp = phi_new[lev]->nComp();
    const int nghost = phi_new[lev]->nGrow();

#if __cplusplus >= 201402L
    auto new_state = std::make_unique<MultiFab>(ba, dm, ncomp, nghost);
    auto old_state = std::make_unique<MultiFab>(ba, dm, ncomp, nghost);
#else
    std::unique_ptr<MultiFab> new_state(new MultiFab(ba, dm, ncomp, nghost));
    std::unique_ptr<MultiFab> old_state(new MultiFab(ba, dm, ncomp, nghost));
#endif

    FillPatch(lev, time, *new_state, 0, ncomp);

    std::swap(new_state, phi_new[lev]);
    std::swap(old_state, phi_old[lev]);

    t_new = time;
    t_old = time - 1.e200;

    if (lev > 0 && do_reflux) {
        flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp));
    }    
}



// compute a new multifab by coping in phi from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
void
Maestro::FillPatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{
    if (lev == 0)
    {
        Array<MultiFab*> smf;
        Array<Real> stime;
        GetData(0, time, smf, stime);

        MaestroPhysBCFunct physbc(geom[lev],bcs,BndryFunctBase(phifill));
        amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
                                    geom[lev], physbc);
    }
    else
    {
        Array<MultiFab*> cmf, fmf;
        Array<Real> ctime, ftime;
        GetData(lev-1, time, cmf, ctime);
        GetData(lev  , time, fmf, ftime);

        MaestroPhysBCFunct cphysbc(geom[lev-1],bcs,BndryFunctBase(phifill));
        MaestroPhysBCFunct fphysbc(geom[lev  ],bcs,BndryFunctBase(phifill));

        Interpolater* mapper = &cell_cons_interp;

        amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                  0, icomp, ncomp, geom[lev-1], geom[lev],
                                  cphysbc, fphysbc, refRatio(lev-1),
                                  mapper, bcs);
    }
}

// fill an entire multifab by interpolating from the coarser level
// this comes into play when a new level of refinement appears
void
Maestro::FillCoarsePatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{
    BL_ASSERT(lev > 0);

    Array<MultiFab*> cmf;
    Array<Real> ctime;
    GetData(lev-1, time, cmf, ctime);
    
    if (cmf.size() != 1) {
        amrex::Abort("FillCoarsePatch: how did this happen?");
    }

    MaestroPhysBCFunct cphysbc(geom[lev-1],bcs,BndryFunctBase(phifill));
    MaestroPhysBCFunct fphysbc(geom[lev  ],bcs,BndryFunctBase(phifill));

    Interpolater* mapper = &cell_cons_interp;
    amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],
                                 cphysbc, fphysbc, refRatio(lev-1),
                                 mapper, bcs);
}

// utility to copy in data from phi_old and/or phi_new into another multifab
void
Maestro::GetData (int lev, Real time, Array<MultiFab*>& data, Array<Real>& datatime)
{
    data.clear();
    datatime.clear();

    const Real teps = (t_new - t_old) * 1.e-3;

    if (time > t_new - teps && time < t_new + teps)
    {
        data.push_back(phi_new[lev].get());
        datatime.push_back(t_new);
    }
    else if (time > t_old - teps && time < t_old + teps)
    {
        data.push_back(phi_old[lev].get());
        datatime.push_back(t_old);
    }
    else
    {
        data.push_back(phi_old[lev].get());
        data.push_back(phi_new[lev].get());
        datatime.push_back(t_old);
        datatime.push_back(t_new);
    }
}

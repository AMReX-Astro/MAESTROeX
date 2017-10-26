
#include <Maestro.H>

using namespace amrex;

// check to see if we need to regrid, then regrid
void
Maestro::Regrid (int& istep)
{
    if (regrid_int > 0)  // We may need to regrid
    {
        if ( (istep-1) % regrid_int == 0)  // if we have hit regrid_int
        {
            // wallclock time
            const Real strt_total = ParallelDescriptor::second();
            
            // regrid could add newly refine levels (if finest_level < max_level)
            // so we save the previous finest level index
            regrid(0, t_new);
            
            // wallclock time
            Real end_total = ParallelDescriptor::second() - strt_total;
            
            // print wallclock time
            ParallelDescriptor::ReduceRealMax(end_total ,ParallelDescriptor::IOProcessorNumber());
            if (Verbose()) {
                Print() << "Time to regrid: " << end_total << '\n';
            }
        }
    }
}

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

// within a call to AmrCore::regrid, this function fills in data at a level
// that existed before, using pre-existing fine and interpolated coarse data
// overrides the pure virtual function in AmrCore
void
Maestro::RemakeLevel (int lev, Real time, const BoxArray& ba,
		      const DistributionMapping& dm)
{
    const int nghost_s = snew[lev]->nGrow();
    const int nghost_u = unew[lev]->nGrow();

#if __cplusplus >= 201402L
    auto snew_state = std::make_unique<MultiFab>(ba, dm,          NSCAL, nghost_s);
    auto sold_state = std::make_unique<MultiFab>(ba, dm,          NSCAL, nghost_s);
    auto unew_state = std::make_unique<MultiFab>(ba, dm, AMREX_SPACEDIM, nghost_u);
    auto uold_state = std::make_unique<MultiFab>(ba, dm, AMREX_SPACEDIM, nghost_u);
#else
    std::unique_ptr<MultiFab> snew_state(new MultiFab(ba, dm,          NSCAL, nghost_s));
    std::unique_ptr<MultiFab> sold_state(new MultiFab(ba, dm,          NSCAL, nghost_s));
    std::unique_ptr<MultiFab> unew_state(new MultiFab(ba, dm, AMREX_SPACEDIM, nghost_u));
    std::unique_ptr<MultiFab> uold_state(new MultiFab(ba, dm, AMREX_SPACEDIM, nghost_u));
#endif

    FillPatch(lev, time, *snew_state, sold, snew, 0, NSCAL, bcs_s);
    std::swap(snew_state, snew[lev]);
    std::swap(sold_state, sold[lev]);

    FillPatch(lev, time, *unew_state, uold, unew, 0, AMREX_SPACEDIM, bcs_u);
    std::swap(unew_state, unew[lev]);
    std::swap(uold_state, uold[lev]);

    t_new = time;
    t_old = time - 1.e200;

    if (lev > 0 && do_reflux) {
        flux_reg_s[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, NSCAL));
        flux_reg_u[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, AMREX_SPACEDIM));
    }    
}

// within a call to AmrCore::regrid, this function fills in data at a level
// that did NOT exist before, using interpolated coarse data
// overrides the pure virtual function in AmrCore
void
Maestro::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
				 const DistributionMapping& dm)
{
    const int nghost_s = snew[lev-1]->nGrow();
    const int nghost_u = unew[lev-1]->nGrow();
    
    snew[lev].reset(new MultiFab(ba, dm,          NSCAL, nghost_s));
    sold[lev].reset(new MultiFab(ba, dm,          NSCAL, nghost_s));
    unew[lev].reset(new MultiFab(ba, dm, AMREX_SPACEDIM, nghost_u));
    uold[lev].reset(new MultiFab(ba, dm, AMREX_SPACEDIM, nghost_u));

    t_new = time;
    t_old = time - 1.e200;

    if (lev > 0 && do_reflux) {
        flux_reg_s[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, NSCAL));
        flux_reg_u[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, AMREX_SPACEDIM));
    }

    FillCoarsePatch(lev, time, *snew[lev], sold, snew, 0,          NSCAL, bcs_s);
    FillCoarsePatch(lev, time, *unew[lev], uold, unew, 0, AMREX_SPACEDIM, bcs_u);
}

// within a call to AmrCore::regrid, this function deletes all data
// at a level of refinement that is no longer needed
// overrides the pure virtual function in AmrCore
void
Maestro::ClearLevel (int lev)
{
    sold[lev].reset(nullptr);
    snew[lev].reset(nullptr);

    uold[lev].reset(nullptr);
    unew[lev].reset(nullptr);

    S_cc_old[lev].reset(nullptr);
    S_cc_new[lev].reset(nullptr);

    gpi[lev].reset(nullptr);

    dSdt[lev].reset(nullptr);

    flux_reg_s[lev].reset(nullptr);
    flux_reg_u[lev].reset(nullptr);
}

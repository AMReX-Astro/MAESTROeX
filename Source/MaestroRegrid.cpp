
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

// within a call to AmrCore::regrid, this function fills in data at a level
// that existed before, using pre-existing fine and interpolated coarse data
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

    FillPatch(lev, time, *new_state, sold, snew, 0, NSCAL, bcs_s);

    std::swap(new_state, snew[lev]);
    std::swap(old_state, sold[lev]);

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
    const int nghost = snew[lev-1]->nGrow();
    
    snew[lev].reset(new MultiFab(ba, dm, NSCAL, nghost));
    sold[lev].reset(new MultiFab(ba, dm, NSCAL, nghost));

    t_new = time;
    t_old = time - 1.e200;

    if (lev > 0 && do_reflux) {
        flux_reg_s[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, NSCAL));
        flux_reg_u[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, AMREX_SPACEDIM));
    }

    FillCoarsePatch(lev, time, *snew[lev], sold, snew, 0, NSCAL, bcs_s);
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

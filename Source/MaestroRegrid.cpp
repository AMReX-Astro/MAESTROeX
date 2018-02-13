
#include <Maestro.H>

using namespace amrex;

// check to see if we need to regrid, then regrid
void
Maestro::Regrid ()
{
    // wallclock time
    const Real strt_total = ParallelDescriptor::second();
            
    // regrid could add newly refine levels (if finest_level < max_level)
    // so we save the previous finest level index
    regrid(0, t_new);
            
    if (spherical == 0) {
        finest_radial_level = finest_level;
        init_multilevel(&finest_level);
        // FIXME
        // we also need to redefine numdisjointchunks, r_start_coord, r_end_coord
        // and "regrid" the base state rho0, rhoh0, tempbar
        // call init_multilevel
        // look at MAESTRO/Source/varden.f90:750-1060
    }

    if (spherical == 1) {
        MakeNormal();
    }
    
    if (evolve_base_state) {
        // force rho0 to be the average of rho - FIXME
    }

    // compute cutoff coordinates
    compute_cutoff_coords(rho0_new.dataPtr());

    // make gravity
    make_grav_cell(grav_cell_new.dataPtr(),
                   rho0_new.dataPtr(),
                   r_cc_loc.dataPtr(),
                   r_edge_loc.dataPtr());

/* FIXME
           ! enforce HSE
           call enforce_HSE(rho0_old,p0_old,grav_cell)

           if (use_tfromp) then
              ! compute full state T = T(rho,p0,X)
              call makeTfromRhoP(sold,p0_old,mla,the_bc_tower%bc_tower_array,dx)
           else
              ! compute full state T = T(rho,h,X)
              call makeTfromRhoH(sold,p0_old,mla,the_bc_tower%bc_tower_array,dx)
           end if

           ! force tempbar to be the average of temp
           call average(mla,sold,tempbar,dx,temp_comp)

           ! gamma1bar needs to be recomputed
           call make_gamma1bar(mla,sold,gamma1bar,p0_old,dx)

           ! beta0_old needs to be recomputed
           call make_beta0(beta0_old,rho0_old,p0_old,gamma1bar,grav_cell)
*/

    // wallclock time
    Real end_total = ParallelDescriptor::second() - strt_total;
            
    // print wallclock time
    ParallelDescriptor::ReduceRealMax(end_total ,ParallelDescriptor::IOProcessorNumber());
    if (maestro_verbose > 0) {
        Print() << "Time to regrid: " << end_total << '\n';
    }
        
}

// tag all cells for refinement
// overrides the pure virtual function in AmrCore
void
Maestro::ErrorEst (int lev, TagBoxArray& tags, Real time, int ng)
{
    if (lev >= temperr.size()) return;

    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

    const Real* dx      = geom[lev].CellSize();

    const MultiFab& state = snew[lev];

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Vector<int>  itags;

	// loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
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
                        ZFILL(dx), &time, &temperr[lev]);
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
    const int ng_s = snew[lev].nGrow();
    const int ng_u = unew[lev].nGrow();
    const int ng_S = S_cc_new[lev].nGrow();
    const int ng_g = gpi[lev].nGrow();
    const int ng_d = dSdt[lev].nGrow();

    MultiFab snew_state    (ba, dm,          Nscal, ng_s);
    MultiFab sold_state    (ba, dm,          Nscal, ng_s);
    MultiFab unew_state    (ba, dm, AMREX_SPACEDIM, ng_u);
    MultiFab uold_state    (ba, dm, AMREX_SPACEDIM, ng_u);
    MultiFab S_cc_new_state(ba, dm,              1, ng_S);
    MultiFab S_cc_old_state(ba, dm,              1, ng_S);
    MultiFab gpi_state     (ba, dm, AMREX_SPACEDIM, ng_g);
    MultiFab dSdt_state    (ba, dm,              1, ng_d);

    FillPatch(lev, time, snew_state, sold, snew, 0, 0, Nscal, 0, bcs_s);
    std::swap(snew_state, snew[lev]);
    std::swap(sold_state, sold[lev]);

    FillPatch(lev, time, unew_state, uold, unew, 0, 0, AMREX_SPACEDIM, 0, bcs_u);
    std::swap(unew_state, unew[lev]);
    std::swap(uold_state, uold[lev]);

    FillPatch(lev, time, S_cc_new_state, S_cc_old, S_cc_new, 0, 0, 1, 0, bcs_f);
    std::swap(S_cc_new_state, S_cc_new[lev]);
    std::swap(S_cc_old_state, S_cc_old[lev]);

    FillPatch(lev, time, gpi_state, gpi, gpi, 0, 0, AMREX_SPACEDIM, 0, bcs_f);
    std::swap(gpi_state, gpi[lev]);

    FillPatch(lev, time, dSdt_state, dSdt, dSdt, 0, 0, 1, 0, bcs_f);
    std::swap(dSdt_state, dSdt[lev]);

    t_new = time;
    t_old = time - 1.e200;

    if (lev > 0 && do_reflux) {
        flux_reg_s[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, Nscal));
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
    snew[lev].define    (ba, dm,          Nscal, 0);
    sold[lev].define    (ba, dm,          Nscal, 0);
    unew[lev].define    (ba, dm, AMREX_SPACEDIM, 0);
    uold[lev].define    (ba, dm, AMREX_SPACEDIM, 0);
    S_cc_new[lev].define(ba, dm,              1, 0);
    S_cc_old[lev].define(ba, dm,              1, 0);
    gpi[lev].define     (ba, dm, AMREX_SPACEDIM, 0);
    dSdt[lev].define    (ba, dm,              1, 0);

    t_new = time;
    t_old = time - 1.e200;

    if (lev > 0 && do_reflux) {
        flux_reg_s[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, Nscal));
        flux_reg_u[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, AMREX_SPACEDIM));
    }

    FillCoarsePatch(lev, time,     snew[lev],     sold,     snew, 0, 0,          Nscal, bcs_s);
    FillCoarsePatch(lev, time,     unew[lev],     uold,     unew, 0, 0, AMREX_SPACEDIM, bcs_u);
    FillCoarsePatch(lev, time, S_cc_new[lev], S_cc_old, S_cc_new, 0, 0,              1, bcs_f);
    FillCoarsePatch(lev, time,      gpi[lev],      gpi,      gpi, 0, 0, AMREX_SPACEDIM, bcs_f);
    FillCoarsePatch(lev, time,     dSdt[lev],     dSdt,     dSdt, 0, 0,              1, bcs_f);
}

// within a call to AmrCore::regrid, this function deletes all data
// at a level of refinement that is no longer needed
// overrides the pure virtual function in AmrCore
void
Maestro::ClearLevel (int lev)
{
    sold[lev].clear();
    snew[lev].clear();

    uold[lev].clear();
    unew[lev].clear();

    S_cc_old[lev].clear();
    S_cc_new[lev].clear();

    gpi[lev].clear();

    dSdt[lev].clear();

    flux_reg_s[lev].reset(nullptr);
    flux_reg_u[lev].reset(nullptr);
}

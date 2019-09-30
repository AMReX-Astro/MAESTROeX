
#include <Maestro.H>

using namespace amrex;

// check to see if we need to regrid, then regrid
void
Maestro::Regrid ()
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Regrid()",Regrid);

    // wallclock time
    const Real strt_total = ParallelDescriptor::second();

    Vector<Real> rho0_temp ( (max_radial_level+1)*nr_fine );
    rho0_temp.shrink_to_fit();

    if (spherical == 0) {
        finest_radial_level = finest_level;

        // look at MAESTRO/Source/varden.f90:750-1060
	// Regrid psi, etarho_cc, etarho_ec, and w0.
	// We do not regrid these in spherical since the base state array only
	// contains 1 level of refinement.
	// We do not regrid these if evolve_base_state=F since they are
	// identically zero, set this way in initialize.
        if (evolve_base_state) {

	    // We must regrid psi, etarho_cc, etarho_ec, and w0
	    // before we call init_multilevel or else we lose
	    // track of where the old valid data was

            // FIXME: may need if statement for irregularly-spaced base states

	    regrid_base_state_cc(psi.dataPtr());
	    regrid_base_state_cc(etarho_cc.dataPtr());
	    regrid_base_state_edge(etarho_ec.dataPtr());
	    regrid_base_state_edge(w0.dataPtr());

        } else {
	    // evolve_base_state == F and spherical == 0

	    // Here we want to fill in the rho0 array so there is
	    // valid data in any new grid locations that are created
	    // during the regrid.
	    for (int i=0; i<rho0_old.size(); ++i) {
		rho0_temp[i] = rho0_old[i];
	    }

	    // We will copy rho0_temp back into the rho0 array after we regrid.
	    regrid_base_state_cc(rho0_temp.dataPtr());

        }

    	// regardless of evolve_base_state, if new grids were
    	// created, we need to initialize tempbar_init there, in
    	// case drive_initial_convection = T
    	regrid_base_state_cc(tempbar_init.dataPtr());
    } else {
        // Here we want to fill in the rho0 array so there is
        // valid data in any new grid locations that are created
        // during the regrid.
        for (int i=0; i<rho0_old.size(); ++i) {
            rho0_temp[i] = rho0_old[i]; 
        }
    }

    // regrid could add newly refine levels (if finest_level < max_level)
    // so we save the previous finest level index
    regrid(0, t_old);

    // Redefine numdisjointchunks, r_start_coord, r_end_coord
    if (spherical == 0) {
	TagArray();
    }
    init_multilevel(tag_array.dataPtr(),&finest_level);

    if (spherical == 1) {
        MakeNormal();
    	if (use_exact_base_state) {
    	    Abort("MaestroRegrid.cpp: need to fill cell_cc_to_r for spherical & exact_base_state");
    	}
    }

    if (evolve_base_state) {
        // force rho0 to be the average of rho
        Average(sold,rho0_old,Rho);
    } else {
    	for (int i=0; i<rho0_old.size(); ++i) {
    	    rho0_old[i] = rho0_temp[i];
    	}
    }

    // compute cutoff coordinates
    compute_cutoff_coords(rho0_old.dataPtr());

    // make gravity
    make_grav_cell(grav_cell_old.dataPtr(),
                   rho0_old.dataPtr(),
                   r_cc_loc.dataPtr(),
                   r_edge_loc.dataPtr());

    // enforce HSE
    enforce_HSE(rho0_old.dataPtr(),
                p0_old.dataPtr(),
                grav_cell_old.dataPtr(),
                r_cc_loc.dataPtr(),
                r_edge_loc.dataPtr());

    if (use_tfromp) {
        // compute full state T = T(rho,p0,X)
        TfromRhoP(sold,p0_old,0);
    } else {
        // compute full state T = T(rho,h,X)
        TfromRhoH(sold,p0_old);
    }

    // force tempbar to be the average of temp
    Average(sold,tempbar,Temp);

    // gamma1bar needs to be recomputed
    MakeGamma1bar(sold,gamma1bar_old,p0_old);

    // beta0_old needs to be recomputed
    if (use_exact_base_state) {
        make_beta0_irreg(beta0_old.dataPtr(), rho0_old.dataPtr(), p0_old.dataPtr(),
                         gamma1bar_old.dataPtr(), grav_cell_old.dataPtr(),
                         r_cc_loc.dataPtr(), r_edge_loc.dataPtr());
    } else {
        make_beta0(beta0_old.dataPtr(), rho0_old.dataPtr(), p0_old.dataPtr(),
                   gamma1bar_old.dataPtr(), grav_cell_old.dataPtr());
    }

    // wallclock time
    Real end_total = ParallelDescriptor::second() - strt_total;

    // print wallclock time
    ParallelDescriptor::ReduceRealMax(end_total,ParallelDescriptor::IOProcessorNumber());
    if (maestro_verbose > 0) {
        Print() << "Time to regrid: " << end_total << '\n';
    }

}

// re-compute tag_array since the actual grid structure changed due to buffering
// this is required in order to compute numdisjointchunks, r_start_coord, r_end_coord
void
Maestro::TagArray ()
{
    // this routine is not required for spherical
    if (spherical == 1) {
        return;
    }

    // timer for profiling
    BL_PROFILE_VAR("Maestro::TagArray()",TagArray);

    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

    for (int lev=1; lev<=max_radial_level; ++lev) {

	const MultiFab& state = sold[lev];

        {
            Vector<int>  itags;

            for (MFIter mfi(state); mfi.isValid(); ++mfi)
            {
                const Box& validBox = mfi.validbox();

                // re-compute tag_array since the actual grid structure changed due to buffering
                // this is required in order to compute numdisjointchunks, r_start_coord, r_end_coord
                retag_array(&tagval, &clearval,
                            ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                            &lev, tag_array.dataPtr());
            }
        }
    }
    ParallelDescriptor::ReduceIntMax(tag_array.dataPtr(),(max_radial_level+1)*nr_fine);
}

// tag all cells for refinement
// overrides the pure virtual function in AmrCore
void
Maestro::ErrorEst (int lev, TagBoxArray& tags, Real time, int ng)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::ErrorEst()",ErrorEst);

    // reset the tag_array (marks radii for planar tagging)
    std::fill(tag_array.begin(), tag_array.end(), 0);

    // convert temperature to perturbation values
    if (use_tpert_in_tagging) {
	PutInPertForm(lev,sold,tempbar,Temp,Temp,bcs_s,true);
    }

    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

    const Real* dx      = geom[lev].CellSize();

    const MultiFab& state = sold[lev];

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Vector<int>  itags;

        for (MFIter mfi(state, true); mfi.isValid(); ++mfi) {

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
	    // for planar problems, we keep track of when a cell at a particular
            // latitude is tagged using tag_array
            state_error(tptr,  ARLIM_3D(tlo), ARLIM_3D(thi),
                        BL_TO_FORTRAN_3D(state[mfi]),
                        &tagval, &clearval,
                        ARLIM_3D(tilebox.loVect()), ARLIM_3D(tilebox.hiVect()),
                        ZFILL(dx), &time,
			&lev, tag_array.dataPtr());

	    //
            // Now update the tags in the TagBox in the tilebox region
            // to be equal to itags
            //
            tagfab.tags_and_untags(itags, tilebox);
        }

        // for planar refinement, we need to gather tagged entries in arrays
        // from all processors and then re-tag tileboxes across each tagged
        // height
        if (spherical == 0) {
            ParallelDescriptor::ReduceIntMax(tag_array.dataPtr(),(max_radial_level+1)*nr_fine);

	    for (MFIter mfi(state, true); mfi.isValid(); ++mfi) {

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

                // tag all cells at a given height if any cells at that height were tagged
                tag_boxes(tptr, ARLIM_3D(tlo), ARLIM_3D(thi),
                          &tagval, &clearval,
                          ARLIM_3D(tilebox.loVect()), ARLIM_3D(tilebox.hiVect()),
                          ZFILL(dx), &time, &lev, tag_array.dataPtr());

                //
                // Now update the tags in the TagBox in the tilebox region
                // to be equal to itags
                //
                tagfab.tags_and_untags(itags, tilebox);
            }
	}
    }

    // convert back to full temperature states
    if (use_tpert_in_tagging) {
        PutInPertForm(lev,sold,tempbar,Temp,Temp,bcs_s,false);
    }
}

// within a call to AmrCore::regrid, this function fills in data at a level
// that existed before, using pre-existing fine and interpolated coarse data
// overrides the pure virtual function in AmrCore
void
Maestro::RemakeLevel (int lev, Real time, const BoxArray& ba,
                      const DistributionMapping& dm)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::RemakeLevel()",RemakeLevel);

    const int ng_s = snew[lev].nGrow();
    const int ng_u = unew[lev].nGrow();
    const int ng_S = S_cc_new[lev].nGrow();
    const int ng_g = gpi[lev].nGrow();
    const int ng_d = dSdt[lev].nGrow();
    const int ng_w = w0_cart[lev].nGrow();
    const int ng_r = rhcc_for_nodalproj[lev].nGrow();
    const int ng_p = pi[lev].nGrow();

    MultiFab snew_state              (ba, dm,          Nscal, ng_s);
    MultiFab sold_state              (ba, dm,          Nscal, ng_s);
    MultiFab unew_state              (ba, dm, AMREX_SPACEDIM, ng_u);
    MultiFab uold_state              (ba, dm, AMREX_SPACEDIM, ng_u);
    MultiFab S_cc_new_state          (ba, dm,              1, ng_S);
    MultiFab S_cc_old_state          (ba, dm,              1, ng_S);
    MultiFab gpi_state               (ba, dm, AMREX_SPACEDIM, ng_g);
    MultiFab dSdt_state              (ba, dm,              1, ng_d);
    MultiFab w0_cart_state           (ba, dm, AMREX_SPACEDIM, ng_w);
    MultiFab rhcc_for_nodalproj_state(ba, dm,              1, ng_r);
    MultiFab pi_state                (convert(ba,nodal_flag), dm, 1, ng_p);

    FillPatch(lev, time, sold_state, sold, sold, 0, 0, Nscal, 0, bcs_s);
    std::swap(sold_state, sold[lev]);
    std::swap(snew_state, snew[lev]);

    FillPatch(lev, time, uold_state, uold, uold, 0, 0, AMREX_SPACEDIM, 0, bcs_u);
    std::swap(uold_state, uold[lev]);
    std::swap(unew_state, unew[lev]);

    FillPatch(lev, time, S_cc_old_state, S_cc_old, S_cc_old, 0, 0, 1, 0, bcs_f);
    std::swap(S_cc_old_state, S_cc_old[lev]);
    std::swap(S_cc_new_state, S_cc_new[lev]);

    FillPatch(lev, time, gpi_state, gpi, gpi, 0, 0, AMREX_SPACEDIM, 0, bcs_f);
    std::swap(gpi_state, gpi[lev]);

    FillPatch(lev, time, dSdt_state, dSdt, dSdt, 0, 0, 1, 0, bcs_f);
    std::swap(dSdt_state, dSdt[lev]);

    std::swap(           w0_cart_state,            w0_cart[lev]);
    std::swap(rhcc_for_nodalproj_state, rhcc_for_nodalproj[lev]);
    std::swap(                pi_state,                 pi[lev]);

    if (spherical == 1) {
        const int ng_n = normal[lev].nGrow();
        const int ng_c = cell_cc_to_r[lev].nGrow();
        MultiFab normal_state(ba, dm, 3, ng_n);
        MultiFab cell_cc_to_r_state(ba, dm, 1, ng_c);
        std::swap(      normal_state,      normal[lev]);
        std::swap(cell_cc_to_r_state,cell_cc_to_r[lev]);
    }

    if (lev > 0 && reflux_type == 2) {
        flux_reg_s[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, Nscal));
    }
}

// within a call to AmrCore::regrid, this function fills in data at a level
// that did NOT exist before, using interpolated coarse data
// overrides the pure virtual function in AmrCore
void
Maestro::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
                                 const DistributionMapping& dm)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeNewLevelFromCoarse()",MakeNewLevelFromCoarse);

    sold[lev].define              (ba, dm,          Nscal, 0);
    snew[lev].define              (ba, dm,          Nscal, 0);
    uold[lev].define              (ba, dm, AMREX_SPACEDIM, 0);
    unew[lev].define              (ba, dm, AMREX_SPACEDIM, 0);
    S_cc_old[lev].define          (ba, dm,              1, 0);
    S_cc_new[lev].define          (ba, dm,              1, 0);
    gpi[lev].define               (ba, dm, AMREX_SPACEDIM, 0);
    dSdt[lev].define              (ba, dm,              1, 0);
    w0_cart[lev].define           (ba, dm, AMREX_SPACEDIM, 1);
    rhcc_for_nodalproj[lev].define(ba, dm,              1, 1);

    pi[lev].define(convert(ba,nodal_flag), dm, 1, 0);     // nodal

    if (spherical == 1) {
        normal      [lev].define(ba, dm, 3, 1);
        cell_cc_to_r[lev].define(ba, dm, 1, 0);
    }

    if (lev > 0 && reflux_type == 2) {
        flux_reg_s[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, Nscal));
    }

    FillCoarsePatch(lev, time,     sold[lev],     sold,     sold, 0, 0,          Nscal, bcs_s);
    FillCoarsePatch(lev, time,     uold[lev],     uold,     uold, 0, 0, AMREX_SPACEDIM, bcs_u, 1);
    FillCoarsePatch(lev, time, S_cc_old[lev], S_cc_old, S_cc_old, 0, 0,              1, bcs_f);
    FillCoarsePatch(lev, time,      gpi[lev],      gpi,      gpi, 0, 0, AMREX_SPACEDIM, bcs_f);
    FillCoarsePatch(lev, time,     dSdt[lev],     dSdt,     dSdt, 0, 0,              1, bcs_f);
}

// within a call to AmrCore::regrid, this function deletes all data
// at a level of refinement that is no longer needed
// overrides the pure virtual function in AmrCore
void
Maestro::ClearLevel (int lev)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::ClearLevel()",ClearLevel);

    sold[lev].clear();
    snew[lev].clear();
    uold[lev].clear();
    unew[lev].clear();
    S_cc_old[lev].clear();
    S_cc_new[lev].clear();
    gpi[lev].clear();
    dSdt[lev].clear();
    w0_cart[lev].clear();
    rhcc_for_nodalproj[lev].clear();
    pi[lev].clear();
    if (spherical == 1) {
        normal[lev].clear();
        cell_cc_to_r[lev].clear();
    }

    flux_reg_s[lev].reset(nullptr);
}

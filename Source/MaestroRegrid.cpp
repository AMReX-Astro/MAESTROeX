
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

// check to see if we need to regrid, then regrid
void Maestro::Regrid() {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Regrid()", Regrid);

    // wallclock time
    const Real strt_total = ParallelDescriptor::second();

    BaseState<Real> rho0_temp(base_geom.max_radial_level + 1,
                              base_geom.nr_fine);

    if (!spherical) {
        base_geom.finest_radial_level = finest_level;

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
            RegridBaseState(psi);
            RegridBaseState(etarho_cc);
            RegridBaseState(etarho_ec, true);
            RegridBaseState(w0, true);

        } else {
            // evolve_base_state == F and !spherical

            // Here we want to fill in the rho0 array so there is
            // valid data in any new grid locations that are created
            // during the regrid.
            rho0_temp.copy(rho0_old);

            // We will copy rho0_temp back into the rho0 array after we regrid.
            RegridBaseState(rho0_temp);
        }

        // regardless of evolve_base_state, if new grids were
        // created, we need to initialize tempbar_init there, in
        // case drive_initial_convection = T
        RegridBaseState(tempbar_init);
    } else {
        // Here we want to fill in the rho0 array so there is
        // valid data in any new grid locations that are created
        // during the regrid.
        rho0_temp.copy(rho0_old);
    }

    // regrid could add newly refine levels (if finest_level < max_level)
    // so we save the previous finest level index
    regrid(0, t_old);

    // Redefine numdisjointchunks, r_start_coord, r_end_coord
    if (!spherical) {
        TagArray();
    }
    BaseState<int> tag_array_b(tag_array, base_geom.max_radial_level + 1,
                               base_geom.nr_fine);
    base_geom.InitMultiLevel(finest_level, tag_array_b.array());

    if (spherical) {
        MakeNormal();
        if (use_exact_base_state) {
            Abort(
                "MaestroRegrid.cpp: need to fill cell_cc_to_r for spherical & "
                "exact_base_state");
        }
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
        w0_cart[lev].setVal(0.);
    }
    // put w0 on Cartesian cell-centers
    Put1dArrayOnCart(w0, w0_cart, true, true, bcs_u, 0, 1);

    if (evolve_base_state) {
        // force rho0 to be the average of rho
        Average(sold, rho0_old, Rho);
    } else {
        rho0_old.copy(rho0_temp);
    }

    // compute cutoff coordinates
    ComputeCutoffCoords(rho0_old);
    base_geom.ComputeCutoffCoords(rho0_old.array());

    // make gravity
    MakeGravCell(grav_cell_old, rho0_old);

    // enforce HSE
    EnforceHSE(rho0_old, p0_old, grav_cell_old);

    if (use_tfromp) {
        // compute full state T = T(rho,p0,X)
        TfromRhoP(sold, p0_old, false);
    } else {
        // compute full state T = T(rho,h,X)
        TfromRhoH(sold, p0_old);
    }

    // force tempbar to be the average of temp
    Average(sold, tempbar, Temp);

    // gamma1bar needs to be recomputed
    MakeGamma1bar(sold, gamma1bar_old, p0_old);

    // beta0_old needs to be recomputed
    MakeBeta0(beta0_old, rho0_old, p0_old, gamma1bar_old, grav_cell_old);

    // wallclock time
    Real end_total = ParallelDescriptor::second() - strt_total;

    // print wallclock time
    ParallelDescriptor::ReduceRealMax(end_total,
                                      ParallelDescriptor::IOProcessorNumber());
    if (maestro_verbose > 0) {
        Print() << "Time to regrid: " << end_total << '\n';
    }
}

// re-compute tag_array since the actual grid structure changed due to buffering
// this is required in order to compute numdisjointchunks, r_start_coord, r_end_coord
void Maestro::TagArray() {
    // this routine is not required for spherical
    if (spherical) {
        return;
    }

    // grids have not been initialized to tag yet.
    if (finest_level == 0) {
        return;
    }

    // timer for profiling
    BL_PROFILE_VAR("Maestro::TagArray()", TagArray);

    for (int lev = 1; lev <= finest_level; ++lev) {
        for (MFIter mfi(sold[lev], false); mfi.isValid(); ++mfi) {
            const Box& validBox = mfi.validbox();
            // re-compute tag_array since the actual grid structure changed due to buffering
            // this is required in order to compute numdisjointchunks, r_start_coord, r_end_coord
            RetagArray(validBox, lev);
        }
    }
    ParallelDescriptor::ReduceIntMax(
        tag_array.dataPtr(),
        (base_geom.max_radial_level + 1) * base_geom.nr_fine);
}

// tag all cells for refinement
// overrides the pure virtual function in AmrCore
void Maestro::ErrorEst(int lev, TagBoxArray& tags, Real time, [[maybe_unused]] int ng) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::ErrorEst()", ErrorEst);

    // reset the tag_array (marks radii for planar tagging)
    std::fill(tag_array.begin(), tag_array.end(), 0);

    // convert temperature to perturbation values
    if (use_tpert_in_tagging) {
        PutInPertForm(lev, sold, tempbar, Temp, Temp, bcs_s, true);
    }

    // if you add openMP here, make sure to collect tag_array across threads
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(sold[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        // tag cells for refinement
        // for planar problems, we keep track of when a cell at a particular
        // latitude is tagged using tag_array
        StateError(tags, sold[lev], mfi, lev, time);
    }

    // for planar refinement, we need to gather tagged entries in arrays
    // from all processors and then re-tag tileboxes across each tagged
    // height
    if (!spherical) {
        ParallelDescriptor::ReduceIntMax(
            tag_array.dataPtr(),
            (base_geom.max_radial_level + 1) * base_geom.nr_fine);

        // NOTE: adding OpenMP breaks the code - not exactly sure why
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(sold[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // tag all cells at a given height if any cells at that height were tagged
            TagBoxes(tags, mfi, lev, time);
        }
    }  // if (!spherical)

    // convert back to full temperature states
    if (use_tpert_in_tagging) {
        PutInPertForm(lev, sold, tempbar, Temp, Temp, bcs_s, false);
    }
}

// within a call to AmrCore::regrid, this function fills in data at a level
// that existed before, using pre-existing fine and interpolated coarse data
// overrides the pure virtual function in AmrCore
void Maestro::RemakeLevel(int lev, Real time, const BoxArray& ba,
                          const DistributionMapping& dm) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::RemakeLevel()", RemakeLevel);

    const int ng_snew = snew[lev].nGrow();
    const int ng_u = unew[lev].nGrow();
    const int ng_S = S_cc_new[lev].nGrow();
    const int ng_g = gpi[lev].nGrow();
    const int ng_d = dSdt[lev].nGrow();
    const int ng_w = w0_cart[lev].nGrow();
    const int ng_r = rhcc_for_nodalproj[lev].nGrow();
    const int ng_p = pi[lev].nGrow();
#ifdef SDC
    const int ng_i = intra[lev].nGrow();
#endif

    MultiFab snew_state(ba, dm, Nscal, ng_snew);
    MultiFab sold_state(ba, dm, Nscal, ng_snew);
    MultiFab unew_state(ba, dm, AMREX_SPACEDIM, ng_u);
    MultiFab uold_state(ba, dm, AMREX_SPACEDIM, ng_u);
    MultiFab S_cc_new_state(ba, dm, 1, ng_S);
    MultiFab S_cc_old_state(ba, dm, 1, ng_S);
    MultiFab gpi_state(ba, dm, AMREX_SPACEDIM, ng_g);
    MultiFab dSdt_state(ba, dm, 1, ng_d);
    MultiFab w0_cart_state(ba, dm, AMREX_SPACEDIM, ng_w);
    MultiFab rhcc_for_nodalproj_state(ba, dm, 1, ng_r);
    MultiFab pi_state(convert(ba, nodal_flag), dm, 1, ng_p);
#ifdef SDC
    MultiFab intra_state(ba, dm, Nscal, ng_i);
#endif

    FillPatch(lev, time, sold_state, sold, sold, 0, 0, Nscal, 0, bcs_s);
    std::swap(sold_state, sold[lev]);
    std::swap(snew_state, snew[lev]);

    FillPatch(lev, time, uold_state, uold, uold, 0, 0, AMREX_SPACEDIM, 0,
              bcs_u);
    std::swap(uold_state, uold[lev]);
    std::swap(unew_state, unew[lev]);

    FillPatch(lev, time, S_cc_old_state, S_cc_old, S_cc_old, 0, 0, 1, 0, bcs_f);
    std::swap(S_cc_old_state, S_cc_old[lev]);
    std::swap(S_cc_new_state, S_cc_new[lev]);

    FillPatch(lev, time, gpi_state, gpi, gpi, 0, 0, AMREX_SPACEDIM, 0, bcs_f);
    std::swap(gpi_state, gpi[lev]);

    FillPatch(lev, time, dSdt_state, dSdt, dSdt, 0, 0, 1, 0, bcs_f);
    std::swap(dSdt_state, dSdt[lev]);

    std::swap(w0_cart_state, w0_cart[lev]);
    std::swap(rhcc_for_nodalproj_state, rhcc_for_nodalproj[lev]);
    std::swap(pi_state, pi[lev]);
#ifdef SDC
    FillPatch(lev, time, intra_state, intra, intra, 0, 0, Nscal, 0, bcs_f);
    std::swap(intra_state, intra[lev]);
#endif

    if (spherical) {
        const int ng_n = normal[lev].nGrow();
        const int ng_c = cell_cc_to_r[lev].nGrow();
        MultiFab normal_state(ba, dm, 3, ng_n);
        iMultiFab cell_cc_to_r_state(ba, dm, 1, ng_c);
        std::swap(normal_state, normal[lev]);
        std::swap(cell_cc_to_r_state, cell_cc_to_r[lev]);
    }

    if (lev > 0 && reflux_type == 2) {
        flux_reg_s[lev] = std::make_unique<FluxRegister>(
            ba, dm, refRatio(lev - 1), lev, Nscal);
    }
}

// within a call to AmrCore::regrid, this function fills in data at a level
// that did NOT exist before, using interpolated coarse data
// overrides the pure virtual function in AmrCore
void Maestro::MakeNewLevelFromCoarse(int lev, Real time, const BoxArray& ba,
                                     const DistributionMapping& dm) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeNewLevelFromCoarse()", MakeNewLevelFromCoarse);

    sold[lev].define(ba, dm, Nscal, 0);
    snew[lev].define(ba, dm, Nscal, 0);
    uold[lev].define(ba, dm, AMREX_SPACEDIM, 0);
    unew[lev].define(ba, dm, AMREX_SPACEDIM, 0);
    S_cc_old[lev].define(ba, dm, 1, 0);
    S_cc_new[lev].define(ba, dm, 1, 0);
    gpi[lev].define(ba, dm, AMREX_SPACEDIM, 0);
    dSdt[lev].define(ba, dm, 1, 0);
    w0_cart[lev].define(ba, dm, AMREX_SPACEDIM, 2);
    rhcc_for_nodalproj[lev].define(ba, dm, 1, 1);

    pi[lev].define(convert(ba, nodal_flag), dm, 1, 0);  // nodal
#ifdef SDC
    intra[lev].define(ba, dm, Nscal, 0);
#endif

    if (spherical) {
        normal[lev].define(ba, dm, 3, 1);
        cell_cc_to_r[lev].define(ba, dm, 1, 0);
    }

    if (lev > 0 && reflux_type == 2) {
        flux_reg_s[lev] = std::make_unique<FluxRegister>(
            ba, dm, refRatio(lev - 1), lev, Nscal);
    }

    FillCoarsePatch(lev, time, sold[lev], sold, sold, 0, 0, Nscal, bcs_s);
    FillCoarsePatch(lev, time, uold[lev], uold, uold, 0, 0, AMREX_SPACEDIM,
                    bcs_u, 1);
    FillCoarsePatch(lev, time, S_cc_old[lev], S_cc_old, S_cc_old, 0, 0, 1,
                    bcs_f);
    FillCoarsePatch(lev, time, gpi[lev], gpi, gpi, 0, 0, AMREX_SPACEDIM, bcs_f);
    FillCoarsePatch(lev, time, dSdt[lev], dSdt, dSdt, 0, 0, 1, bcs_f);
#ifdef SDC
    FillCoarsePatch(lev, time, intra[lev], intra, intra, 0, 0, Nscal, bcs_f);
#endif
}

// within a call to AmrCore::regrid, this function deletes all data
// at a level of refinement that is no longer needed
// overrides the pure virtual function in AmrCore
void Maestro::ClearLevel(int lev) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::ClearLevel()", ClearLevel);

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
#ifdef SDC
    intra[lev].clear();
#endif
    if (spherical) {
        normal[lev].clear();
        cell_cc_to_r[lev].clear();
    }

    flux_reg_s[lev].reset(nullptr);
}

void Maestro::RegridBaseState(BaseState<Real>& base_s, const bool is_edge) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::RegridBaseState()", RegridBaseState);

    const int max_lev = base_geom.max_radial_level + 1;

    const int nrf = is_edge ? base_geom.nr_fine + 1 : base_geom.nr_fine;
    BaseState<Real> state_temp_s(max_lev, nrf);
    auto state_temp = state_temp_s.array();
    auto base = base_s.array();

    // copy the coarsest level of the real arrays into the
    // temp arrays
    ParallelFor(nrf,
                [=] AMREX_GPU_DEVICE(int r) { state_temp(0, r) = base(0, r); });
    Gpu::synchronize();

    // piecewise linear interpolation to fill the cc temp arrays
    for (auto n = 1; n < max_lev; ++n) {
        if (is_edge) {
            const auto nrn = base_geom.nr(n) + 1;
            ParallelFor(nrn, [=] AMREX_GPU_DEVICE(int r) {
                if (r % 2 == 0) {
                    state_temp(n, r) = state_temp(n - 1, r / 2);
                } else {
                    state_temp(n, r) =
                        0.5 * (state_temp(n - 1, r / 2) +
                               0.25 * state_temp(n - 1, r / 2 + 1));
                }
            });
        } else {
            const auto nrn = base_geom.nr(n);
            ParallelFor(nrn, [=] AMREX_GPU_DEVICE(int r) {
                if (r == 0 || r == nrn - 1) {
                    state_temp(n, r) = state_temp(n - 1, r / 2);
                } else {
                    if (r % 2 == 0) {
                        state_temp(n, r) = 0.75 * state_temp(n - 1, r / 2) +
                                           0.25 * state_temp(n - 1, r / 2 - 1);
                    } else {
                        state_temp(n, r) = 0.75 * state_temp(n - 1, r / 2) +
                                           0.25 * state_temp(n - 1, r / 2 + 1);
                    }
                }
            });
        }
        Gpu::synchronize();
    }

    // copy valid data into temp
    for (auto n = 1; n <= finest_level; ++n) {
        for (auto i = 1; i <= base_geom.numdisjointchunks(n); ++i) {
            const auto lo = base_geom.r_start_coord(n, i);
            const auto hi = is_edge ? base_geom.r_end_coord(n, i) + 1
                                    : base_geom.r_end_coord(n, i);
            ParallelFor(hi - lo + 1, [=] AMREX_GPU_DEVICE(int k) {
                int r = k + lo;
                state_temp(n, r) = base(n, r);
            });
            Gpu::synchronize();
        }
    }

    // copy temp array back into the real thing
    base_s.copy(state_temp_s);
}

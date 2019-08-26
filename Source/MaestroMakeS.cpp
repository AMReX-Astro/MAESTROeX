
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

// compute S at cell-centers
void
Maestro::Make_S_cc (Vector<MultiFab>& S_cc,
                    Vector<MultiFab>& delta_gamma1_term,
                    Vector<MultiFab>& delta_gamma1,
                    const Vector<MultiFab>& scal,
                    const Vector<MultiFab>& u,
                    const Vector<MultiFab>& rho_omegadot,
                    const Vector<MultiFab>& rho_Hnuc,
                    const Vector<MultiFab>& rho_Hext,
                    const Vector<MultiFab>& thermal,
                    const RealVector& p0,
                    const RealVector& gamma1bar,
                    RealVector& delta_gamma1_termbar,
                    const RealVector& psi_in)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Make_S_cc()", Make_S_cc);

#ifdef AMREX_USE_CUDA
    auto not_launched = Gpu::notInLaunchRegion();
    // turn on GPU
    if (not_launched) Gpu::setLaunchRegion(true);
#endif

    // put 1d base state quantities on cartestian grid for spherical case
    Vector<MultiFab> gamma1bar_cart(finest_level+1);
    Vector<MultiFab> p0_cart(finest_level+1);
    Vector<MultiFab> gradp0_cart(finest_level+1);
    Vector<MultiFab> psi_cart(finest_level+1);

    // calculate gradp0
    RealVector gradp0;

    if (spherical == 1) {
        gradp0.resize((max_radial_level+1)*nr_fine);
        if (use_delta_gamma1_term) {
            Real dr = r_cc_loc[1] - r_cc_loc[0];
            gradp0[0] = (p0[1] - p0[0]) / dr;

            dr = r_cc_loc[nr_fine-1] - r_cc_loc[nr_fine-2];
            gradp0[nr_fine-1] = (p0[nr_fine-1] - p0[nr_fine-2]) / dr;

            for (int r=1; r < nr_fine-1; r++) {
                dr = r_cc_loc[r+1] - r_cc_loc[r-1];
                gradp0[r] = (p0[r+1] - p0[r-1]) / dr;
            }
        }

        for (int lev=0; lev<=finest_level; ++lev) {
            gamma1bar_cart[lev].define(grids[lev], dmap[lev], 1, 0);
            p0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
            gradp0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
            psi_cart[lev].define(grids[lev], dmap[lev], 1, 0);

            gamma1bar_cart[lev].setVal(0.);
            p0_cart[lev].setVal(0.);
            gradp0_cart[lev].setVal(0.);
            psi_cart[lev].setVal(0.);
        }

        if (use_delta_gamma1_term) {
            Put1dArrayOnCart(gamma1bar,gamma1bar_cart,0,0,bcs_f,0);
            Put1dArrayOnCart(p0,p0_cart,0,0,bcs_f,0);
            Put1dArrayOnCart(gradp0,gradp0_cart,0,0,bcs_f,0);
            Put1dArrayOnCart(psi_in,psi_cart,0,0,bcs_f,0);
        }
    }


    for (int lev=0; lev<=finest_level; ++lev) {

        // Declare local storage now. This should be done outside the MFIter loop,
        // and then we will resize the Fabs in each MFIter loop iteration. Then,
        // we apply an Elixir to ensure that their memory is saved until it is no
        // longer needed (only relevant for the asynchronous case, usually on GPUs).

        // get references to the MultiFabs at level lev
        MultiFab& S_cc_mf = S_cc[lev];
        MultiFab& delta_gamma1_term_mf = delta_gamma1_term[lev];
        MultiFab& delta_gamma1_mf = delta_gamma1[lev];
        const MultiFab& scal_mf = scal[lev];
        const MultiFab& u_mf = u[lev];
        const MultiFab& rho_odot_mf = rho_omegadot[lev];
        const MultiFab& rho_Hnuc_mf = rho_Hnuc[lev];
        const MultiFab& rho_Hext_mf = rho_Hext[lev];
        const MultiFab& thermal_mf = thermal[lev];
        const MultiFab& normal_mf = normal[lev];
        const MultiFab& cc_to_r = cell_cc_to_r[lev];
        const Real* dx = geom[lev].CellSize();

        const MultiFab& p0_mf = p0_cart[lev];
        const MultiFab& gradp0_mf = gradp0_cart[lev];
        const MultiFab& gamma1bar_mf = gamma1bar_cart[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(S_cc_mf, true); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.

            if (spherical == 1) {

#pragma gpu box(tileBox)
                make_S_cc_sphr(AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),
                               BL_TO_FORTRAN_ANYD(S_cc_mf[mfi]),
                               BL_TO_FORTRAN_ANYD(delta_gamma1_term_mf[mfi]),
                               BL_TO_FORTRAN_ANYD(delta_gamma1_mf[mfi]),
                               BL_TO_FORTRAN_ANYD(scal_mf[mfi]),
                               BL_TO_FORTRAN_ANYD(u_mf[mfi]),
                               BL_TO_FORTRAN_ANYD(rho_odot_mf[mfi]),
                               BL_TO_FORTRAN_ANYD(rho_Hnuc_mf[mfi]),
                               BL_TO_FORTRAN_ANYD(rho_Hext_mf[mfi]),
                               BL_TO_FORTRAN_ANYD(thermal_mf[mfi]),
                               BL_TO_FORTRAN_ANYD(p0_mf[mfi]),
                               BL_TO_FORTRAN_ANYD(gradp0_mf[mfi]),
                               BL_TO_FORTRAN_ANYD(gamma1bar_mf[mfi]),
                               BL_TO_FORTRAN_ANYD(normal_mf[mfi]));

            } else {
#pragma gpu box(tileBox)
                make_S_cc(AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),
                          lev,
                          BL_TO_FORTRAN_ANYD(S_cc_mf[mfi]),
                          BL_TO_FORTRAN_ANYD(delta_gamma1_term_mf[mfi]),
                          BL_TO_FORTRAN_ANYD(delta_gamma1_mf[mfi]),
                          BL_TO_FORTRAN_ANYD(scal_mf[mfi]),
                          BL_TO_FORTRAN_ANYD(u_mf[mfi]),
                          BL_TO_FORTRAN_ANYD(rho_odot_mf[mfi]),
                          BL_TO_FORTRAN_ANYD(rho_Hnuc_mf[mfi]),
                          BL_TO_FORTRAN_ANYD(rho_Hext_mf[mfi]),
                          BL_TO_FORTRAN_ANYD(thermal_mf[mfi]),
                          p0.dataPtr(), gamma1bar.dataPtr(), AMREX_REAL_ANYD(dx));
            }
        }
    }

    // average fine data onto coarser cells
    AverageDown(S_cc,0,1);
    AverageDown(delta_gamma1_term,0,1);

    if (use_delta_gamma1_term) {

        // horizontal average of delta_gamma1_term
        Average(delta_gamma1_term,delta_gamma1_termbar,0);

        for (int lev=0; lev<=finest_level; ++lev) {

            // get references to the MultiFabs at level lev
            MultiFab& delta_gamma1_term_mf = delta_gamma1_term[lev];
            const MultiFab& delta_gamma1_mf = delta_gamma1[lev];
            const Real* dx = geom[lev].CellSize();

            const MultiFab& p0_mf = p0_cart[lev];
            const MultiFab& psi_mf = psi_cart[lev];
            const MultiFab& gamma1bar_mf = gamma1bar_cart[lev];

            // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
            for ( MFIter mfi(delta_gamma1_term_mf, true); mfi.isValid(); ++mfi ) {

                // Get the index space of the valid region
                const Box& tileBox = mfi.tilebox();

                // call fortran subroutine
                // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
                // lo/hi coordinates (including ghost cells), and/or the # of components
                // We will also pass "validBox", which specifies the "valid" region.
                if (spherical == 1) {
#pragma gpu box(tileBox)
                    create_correction_delta_gamma1_term_sphr(AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),
                                                             BL_TO_FORTRAN_ANYD(delta_gamma1_term_mf[mfi]),
                                                             BL_TO_FORTRAN_ANYD(delta_gamma1_mf[mfi]),
                                                             BL_TO_FORTRAN_ANYD(gamma1bar_mf[mfi]),
                                                             BL_TO_FORTRAN_ANYD(psi_mf[mfi]),
                                                             BL_TO_FORTRAN_ANYD(p0_mf[mfi]));
                } else {
#pragma gpu box(tileBox)
                    create_correction_delta_gamma1_term(AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),
                                                        lev,
                                                        BL_TO_FORTRAN_ANYD(delta_gamma1_term_mf[mfi]),
                                                        BL_TO_FORTRAN_ANYD(delta_gamma1_mf[mfi]),
                                                        gamma1bar.dataPtr(), psi_in.dataPtr(), p0.dataPtr());

                }
            }
        }

    }

#ifdef AMREX_USE_CUDA
    // turn off GPU
    if (not_launched) Gpu::setLaunchRegion(false);
#endif

}

// compute rhcc = beta0*(S_cc-Sbar) + beta0*delta_chi
void
Maestro::MakeRHCCforNodalProj (Vector<MultiFab>& rhcc,
                               const Vector<MultiFab>& S_cc,
                               const RealVector& Sbar,
                               const RealVector& beta0,
                               const Vector<MultiFab>& delta_gamma1_term)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeRHCCforNodalProj()",MakeRHCCforNodalProj);

#ifdef AMREX_USE_CUDA
    auto not_launched = Gpu::notInLaunchRegion();
    // turn on GPU
    if (not_launched) Gpu::setLaunchRegion(true);
#endif

    Vector<MultiFab> Sbar_cart(finest_level+1);
    Vector<MultiFab> beta0_cart(finest_level+1);

    if (spherical == 1) {
        for (int lev=0; lev<=finest_level; ++lev) {
            Sbar_cart [lev].define(grids[lev], dmap[lev], 1, 0);
            beta0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
            Sbar_cart [lev].setVal(0.);
            beta0_cart[lev].setVal(0.);
        }

        Put1dArrayOnCart(Sbar,Sbar_cart,0,0,bcs_f,0);
        Put1dArrayOnCart(beta0,beta0_cart,0,0,bcs_f,0);
    }

    for (int lev=0; lev<=finest_level; ++lev) {

        // fill rhcc
        // get references to the MultiFabs at level lev
        MultiFab& rhcc_mf = rhcc[lev];
        const MultiFab& S_cc_mf = S_cc[lev];
        const MultiFab& Sbar_mf = Sbar_cart[lev];
        const MultiFab& beta0_mf = beta0_cart[lev];

        const MultiFab& delta_gamma1_term_mf = delta_gamma1_term[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(S_cc_mf, true); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            if (spherical == 1) {

#pragma gpu box(tileBox)
                make_rhcc_for_nodalproj_sphr(AMREX_INT_ANYD(tileBox.loVect()),
                                             AMREX_INT_ANYD(tileBox.hiVect()),
                                             BL_TO_FORTRAN_ANYD(rhcc_mf[mfi]),
                                             BL_TO_FORTRAN_ANYD(S_cc_mf[mfi]),
                                             BL_TO_FORTRAN_ANYD(Sbar_mf[mfi]),
                                             BL_TO_FORTRAN_ANYD(beta0_mf[mfi]),
                                             BL_TO_FORTRAN_ANYD(delta_gamma1_term_mf[mfi]));
            } else {
#pragma gpu box(tileBox)
                make_rhcc_for_nodalproj(AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),
                                        lev,
                                        BL_TO_FORTRAN_ANYD(rhcc_mf[mfi]),
                                        BL_TO_FORTRAN_ANYD(S_cc_mf[mfi]),
                                        Sbar.dataPtr(), beta0.dataPtr(),
                                        BL_TO_FORTRAN_ANYD(delta_gamma1_term_mf[mfi]));
            }
        }
    }

    // averge down and fill ghost cells using first-order extrapolation
    AverageDown(rhcc,0,1);
    FillPatch(t_old, rhcc, rhcc, rhcc, 0, 0, 1, 0, bcs_f);

#ifdef AMREX_USE_CUDA
    // turn off GPU
    if (not_launched) Gpu::setLaunchRegion(false);
#endif
}

void
Maestro::CorrectRHCCforNodalProj(Vector<MultiFab>& rhcc,
                                 const RealVector& rho0,
                                 const RealVector& beta0,
                                 const RealVector& gamma1bar,
                                 const RealVector& p0,
                                 const Vector<MultiFab>& delta_p_term)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::CorrectRHCCforNodalProj()",CorrectRHCCforNodalProj);

#ifdef AMREX_USE_CUDA
    auto not_launched = Gpu::notInLaunchRegion();
    // turn on GPU
    if (not_launched) Gpu::setLaunchRegion(true);
#endif

    // Local variables
    Vector<MultiFab> correction_cc(finest_level+1);

    for (int lev=0; lev<=finest_level; ++lev) {
        correction_cc[lev].define(grids[lev], dmap[lev], 1, 1);
        correction_cc[lev].setVal(0.);
    }

    Vector<MultiFab> gamma1bar_cart(finest_level+1);
    Vector<MultiFab>        p0_cart(finest_level+1);
    Vector<MultiFab>     beta0_cart(finest_level+1);
    Vector<MultiFab>      rho0_cart(finest_level+1);

    if (spherical == 1) {
        for (int lev=0; lev<=finest_level; ++lev) {
            gamma1bar_cart[lev].define(grids[lev], dmap[lev], 1, 0);
            p0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
            beta0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
            rho0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
            gamma1bar_cart[lev].setVal(0.);
            p0_cart[lev].setVal(0.);
            beta0_cart[lev].setVal(0.);
            rho0_cart[lev].setVal(0.);
        }

        Put1dArrayOnCart(gamma1bar,gamma1bar_cart,0,0,bcs_f,0);
        Put1dArrayOnCart(p0,p0_cart,0,0,bcs_f,0);
        Put1dArrayOnCart(beta0,beta0_cart,0,0,bcs_f,0);
        Put1dArrayOnCart(rho0,rho0_cart,0,0,bcs_s,Rho);
    }

    for (int lev=0; lev<=finest_level; ++lev) {
        // get references to the MultiFabs at level lev
        const MultiFab& delta_p_mf = delta_p_term[lev];
        MultiFab& correction_cc_mf = correction_cc[lev];

        const MultiFab& gamma1bar_mf = gamma1bar_cart[lev];
        const MultiFab& p0_mf = p0_cart[lev];
        const MultiFab& beta0_mf = beta0_cart[lev];
        const MultiFab& rho0_mf = rho0_cart[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(correction_cc_mf, true); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            if (spherical == 1) {
#pragma gpu box(tileBox)
                create_correction_cc_sphr(AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),
                                          BL_TO_FORTRAN_ANYD(correction_cc_mf[mfi]),
                                          BL_TO_FORTRAN_ANYD(delta_p_mf[mfi]),
                                          BL_TO_FORTRAN_ANYD(beta0_mf[mfi]),
                                          BL_TO_FORTRAN_ANYD(gamma1bar_mf[mfi]),
                                          BL_TO_FORTRAN_ANYD(p0_mf[mfi]),
                                          BL_TO_FORTRAN_ANYD(rho0_mf[mfi]), dt);
            } else {
#pragma gpu box(tileBox)
                create_correction_cc(AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()), lev,
                                     BL_TO_FORTRAN_ANYD(correction_cc_mf[mfi]),
                                     BL_TO_FORTRAN_ANYD(delta_p_mf[mfi]),
                                     beta0.dataPtr(), gamma1bar.dataPtr(),
                                     p0.dataPtr(), dt);
            }
        }
    }

    // average down and fill ghost cells using first-order extrapolation
    AverageDown(correction_cc,0,1);
    FillPatch(t_old, correction_cc, correction_cc, correction_cc, 0, 0, 1, 0, bcs_f);

    // add correction term
    for (int lev=0; lev<=finest_level; ++lev) {
        MultiFab::Add(rhcc[lev],correction_cc[lev],0,0,1,1);
    }

#ifdef AMREX_USE_CUDA
    // turn off GPU
    if (not_launched) Gpu::setLaunchRegion(false);
#endif
}

// compute rhcc = beta0*(S_cc-Sbar) + beta0*delta_chi
void
Maestro::MakeRHCCforMacProj (Vector<MultiFab>& rhcc,
                             const RealVector& rho0,
                             const Vector<MultiFab>& S_cc,
                             const RealVector& Sbar,
                             const RealVector& beta0,
                             const Vector<MultiFab>& delta_gamma1_term,
                             const RealVector& gamma1bar,
                             const RealVector& p0,
                             const Vector<MultiFab>& delta_p_term,
                             Vector<MultiFab>& delta_chi,
                             int is_predictor)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeRHCCforMacProj()",MakeRHCCforMacProj);

#ifdef AMREX_USE_CUDA
    auto not_launched = Gpu::notInLaunchRegion();
    // turn on GPU
    if (not_launched) Gpu::setLaunchRegion(true);
#endif

    // put 1d base state quantities on cartestian grid for spherical case
    Vector<MultiFab> Sbar_cart(finest_level+1);
    Vector<MultiFab> beta0_cart(finest_level+1);
    Vector<MultiFab> gamma1bar_cart(finest_level+1);
    Vector<MultiFab> p0_cart(finest_level+1);
    Vector<MultiFab> rho0_cart(finest_level+1);

    if (spherical == 1) {
        for (int lev=0; lev<=finest_level; ++lev) {
            Sbar_cart [lev].define(grids[lev], dmap[lev], 1, 0);
            beta0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
            gamma1bar_cart[lev].define(grids[lev], dmap[lev], 1, 0);
            p0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
            rho0_cart[lev].define(grids[lev], dmap[lev], 1, 0);

            Sbar_cart [lev].setVal(0.);
            beta0_cart[lev].setVal(0.);
            gamma1bar_cart[lev].setVal(0.);
            p0_cart[lev].setVal(0.);
            rho0_cart[lev].setVal(0.);
        }

        Put1dArrayOnCart(Sbar,Sbar_cart,0,0,bcs_f,0);
        Put1dArrayOnCart(beta0,beta0_cart,0,0,bcs_f,0);

        if (dpdt_factor > 0.0) {
            Put1dArrayOnCart(gamma1bar,gamma1bar_cart,0,0,bcs_f,0);
            Put1dArrayOnCart(p0,p0_cart,0,0,bcs_f,0);
            Put1dArrayOnCart(rho0,rho0_cart,0,0,bcs_f,0);
        }
    }

    for (int lev=0; lev<=finest_level; ++lev) {

        // fill rhcc
        // get references to the MultiFabs at level lev
        MultiFab& rhcc_mf = rhcc[lev];
        const MultiFab& S_cc_mf = S_cc[lev];
        const MultiFab& delta_gamma1_term_mf = delta_gamma1_term[lev];
        const MultiFab& delta_p_mf = delta_p_term[lev];
        MultiFab& delta_chi_mf = delta_chi[lev];

        const MultiFab& Sbar_mf = Sbar_cart[lev];
        const MultiFab& beta0_mf = beta0_cart[lev];
        const MultiFab& gamma1bar_mf = gamma1bar_cart[lev];
        const MultiFab& p0_mf = p0_cart[lev];
        const MultiFab& rho0_mf = rho0_cart[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(S_cc_mf, true); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            if (spherical == 1) {
                const Real* dx = geom[lev].CellSize();
#pragma gpu box(tileBox)
                make_rhcc_for_macproj_sphr(AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),
                                           BL_TO_FORTRAN_ANYD(rhcc_mf[mfi]),
                                           BL_TO_FORTRAN_ANYD(S_cc_mf[mfi]),
                                           BL_TO_FORTRAN_ANYD(Sbar_mf[mfi]),
                                           BL_TO_FORTRAN_ANYD(beta0_mf[mfi]),
                                           BL_TO_FORTRAN_ANYD(rho0_mf[mfi]),
                                           BL_TO_FORTRAN_ANYD(delta_gamma1_term_mf[mfi]),
                                           BL_TO_FORTRAN_ANYD(gamma1bar_mf[mfi]),
                                           BL_TO_FORTRAN_ANYD(p0_mf[mfi]),
                                           BL_TO_FORTRAN_ANYD(delta_p_mf[mfi]),
                                           BL_TO_FORTRAN_ANYD(delta_chi_mf[mfi]),
                                           dt, is_predictor);
            } else {
#pragma gpu box(tileBox)
                make_rhcc_for_macproj(AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),lev,
                                      BL_TO_FORTRAN_ANYD(rhcc_mf[mfi]),
                                      BL_TO_FORTRAN_ANYD(S_cc_mf[mfi]),
                                      Sbar.dataPtr(), beta0.dataPtr(),
                                      BL_TO_FORTRAN_ANYD(delta_gamma1_term_mf[mfi]),
                                      gamma1bar.dataPtr(), p0.dataPtr(),
                                      BL_TO_FORTRAN_ANYD(delta_p_mf[mfi]),
                                      BL_TO_FORTRAN_ANYD(delta_chi_mf[mfi]),
                                      dt, is_predictor);
            }
        }
    }

#ifdef AMREX_USE_CUDA
    // turn off GPU
    if (not_launched) Gpu::setLaunchRegion(false);
#endif

}

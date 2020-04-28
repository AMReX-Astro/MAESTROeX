
#include <Maestro.H>
#include <Maestro_F.H>
#include <Maestro_F.H>

using namespace amrex;

#ifndef SDC
void Maestro::Burner(const Vector<MultiFab>& s_in,
                     Vector<MultiFab>& s_out,
                     const Vector<MultiFab>& rho_Hext,
                     Vector<MultiFab>& rho_omegadot,
                     Vector<MultiFab>& rho_Hnuc,
                     const RealVector& p0,
                     const Real dt_in,
                     const Real time_in)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Burner()", Burner);

    // Put tempbar_init on cart
    Vector<MultiFab> tempbar_init_cart(finest_level+1);

    if (spherical == 1) {
        for (int lev=0; lev<=finest_level; ++lev) {
            tempbar_init_cart[lev].define(grids[lev], dmap[lev], 1, 0);
            tempbar_init_cart[lev].setVal(0.);
        }

        if (drive_initial_convection == 1) {
            Put1dArrayOnCart(tempbar_init,tempbar_init_cart,0,0,bcs_f,0);
        }
    }

    for (int lev = 0; lev <= finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab&         s_in_mf =         s_in[lev];
        MultiFab&        s_out_mf =        s_out[lev];
        const MultiFab&     rho_Hext_mf =     rho_Hext[lev];
        MultiFab& rho_omegadot_mf = rho_omegadot[lev];
        MultiFab&     rho_Hnuc_mf =     rho_Hnuc[lev];
        const MultiFab& tempbar_cart_mf = tempbar_init_cart[lev];

        // create mask assuming refinement ratio = 2
        int finelev = lev+1;
        if (lev == finest_level) finelev = finest_level;
        auto max_lev = max_radial_level + 1;

        const BoxArray& fba = s_in[finelev].boxArray();
        const iMultiFab& mask = makeFineMask(s_in_mf, fba, IntVect(2));

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(s_in[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            bool use_mask = !(lev==finest_level);

            const Array4<const Real> s_in_arr = s_in[lev].array(mfi);
            const Array4<Real> s_out_arr = s_out[lev].array(mfi);
            const Array4<const Real> rho_Hext_arr = rho_Hext[lev].array(mfi);
            const Array4<Real> rho_omegadot_arr = rho_omegadot[lev].array(mfi);
            const Array4<Real> rho_Hnuc_arr = rho_Hnuc[lev].array(mfi);
            const Real * AMREX_RESTRICT tempbar_init_p = tempbar_init.dataPtr();
            const Array4<const int> mask_arr = mask.array(mfi);

            auto ispec_threshold = network_spec_index(burner_threshold_species);

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            if (spherical) {
#pragma gpu box(tileBox)
                burner_loop_sphr(AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),
                                 BL_TO_FORTRAN_ANYD(s_in_mf[mfi]),
                                 BL_TO_FORTRAN_ANYD(s_out_mf[mfi]),
                                 BL_TO_FORTRAN_ANYD(rho_Hext_mf[mfi]),
                                 BL_TO_FORTRAN_ANYD(rho_omegadot_mf[mfi]),
                                 BL_TO_FORTRAN_ANYD(rho_Hnuc_mf[mfi]),
                                 BL_TO_FORTRAN_ANYD(tempbar_cart_mf[mfi]), dt_in, time_in, 
                                 BL_TO_FORTRAN_ANYD(mask[mfi]), use_mask);
            } else {
                ParallelFor(tileBox, 
                [=] AMREX_GPU_DEVICE (int i, int j, int k) 
                {
                    if (use_mask && !mask_arr(i,j,k)) return; // cell is covered by finer cells

                    auto rho = s_in_arr(i,j,k,Rho);
                    Real x_in[NumSpec];
                    for (auto n = 0; n < NumSpec; ++n) {
                        x_in[n] = s_in_arr(i,j,k,FirstSpec+n) / s_in_arr(i,j,k,Rho);
                    }

                    Real T_in = 0.0;
                    if (drive_initial_convection) {
                        auto r = (AMREX_SPACEDIM == 2) ? j : k;
                        T_in = tempbar_init_p[lev+max_lev*r];
                    } else {
                        T_in = s_in_arr(i,j,k,Temp);
                    }

                    Real x_test = (ispec_threshold > 0) ? x_in[ispec_threshold] : 0.0;
                        
                    burn_t state_in;
                    burn_t state_out;

                    Real x_out[NumSpec];
                    Real rhowdot[NumSpec];
                    Real rhoH = 0.0;

                    // if the threshold species is not in the network, then we burn
                    // normally.  if it is in the network, make sure the mass
                    // fraction is above the cutoff.
                    if ((rho > burning_cutoff_density_lo && 
                        rho < burning_cutoff_density_hi) &&
                        ( ispec_threshold < 0 || 
                        (ispec_threshold > 0 && x_test > burner_threshold_cutoff) ) ) {
                        // Initialize burn state_in and state_out
                        state_in.e   = 0.0;
                        state_in.rho = rho;
                        state_in.T   = T_in;
                        for (auto n = 0; n < NumSpec; ++n) {
                            state_in.xn[n] = x_in[n];
                        }

                        // initialize state_out the same as state_in
                        state_out.e   = 0.0;
                        state_out.rho = rho;
                        state_out.T   = T_in;
                        for (auto n = 0; n < NumSpec; ++n) {
                            state_out.xn[n] = x_in[n];
                        }
                        
                        burner(state_out, dt_in);

                        for (auto n = 0; n < NumSpec; ++n) {
                            x_out[n] = state_out.xn[n];
                            rhowdot[n] = state_out.rho * \
                                (state_out.xn[n] - state_in.xn[n]) / dt_in;
                        }
                        rhoH = state_out.rho * (state_out.e - state_in.e) / dt_in;
                    } else {
                        for (auto n = 0; n < NumSpec; ++n) {
                            x_out[n] = x_in[n];
                            rhowdot[n] = 0.0;
                        }
                    }

                    // check if sum{X_k} = 1
                    Real sumX = 0.0;
                    for (auto n = 0; n < NumSpec; ++n) {
                        sumX += x_out[n];
                    }

                    if (fabs(sumX - 1.0) > reaction_sum_tol) {
    #ifndef AMREX_USE_GPU
                        Abort("ERROR: abundances do not sum to 1");
    #endif
                        for (auto n = 0; n < NumSpec; ++n) {
                            state_out.xn[n] /= sumX;
                        }
                    }

                    // pass the density and pi through
                    s_out_arr(i,j,k,Rho) = s_in_arr(i,j,k,Rho);
                    s_out_arr(i,j,k,Pi) = s_in_arr(i,j,k,Pi);

                    // update the species
                    for (auto n = 0; n < NumSpec; ++n) {
                        s_out_arr(i,j,k,FirstSpec+n) = x_out[n] * rho;
                    }

                    // store the energy generation and species create quantities
                    for (auto n = 0; n < NumSpec; ++n) {
                        rho_omegadot_arr(i,j,k,n) = rhowdot[n];
                    }
                    rho_Hnuc_arr(i,j,k) = rhoH;

                    // update the enthalpy -- include the change due to external heating
                    s_out_arr(i,j,k,RhoH) = s_in_arr(i,j,k,RhoH) \
                        + dt_in * rho_Hnuc_arr(i,j,k) + dt_in * rho_Hext_arr(i,j,k);
                });
            }
        }
    }
}

#else
// SDC burner
void Maestro::Burner(const Vector<MultiFab>& s_in,
                     Vector<MultiFab>& s_out,
                     const RealVector& p0,
                     const Real dt_in,
                     const Real time_in,
                     const Vector<MultiFab>& source)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::BurnerSDC()", BurnerSDC);

    // Put tempbar_init on cart
    Vector<MultiFab> p0_cart(finest_level+1);

    if (spherical == 1) {
        for (int lev=0; lev<=finest_level; ++lev) {
            p0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
            p0_cart[lev].setVal(0.);
        }

        if (drive_initial_convection == 1) {
            Put1dArrayOnCart(p0,p0_cart,0,0,bcs_f,0);
        }
    }

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab&    s_in_mf =    s_in[lev];
              MultiFab&   s_out_mf =   s_out[lev];
        const MultiFab& p0_cart_mf = p0_cart[lev];
        const MultiFab&  source_mf =  source[lev];
        
        // create mask assuming refinement ratio = 2
        int finelev = lev+1;
        if (lev == finest_level) finelev = finest_level;

        const BoxArray& fba = s_in[finelev].boxArray();
        const iMultiFab& mask = makeFineMask(s_in_mf, fba, IntVect(2));
        

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(s_in_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            int use_mask = !(lev==finest_level);

            // call fortran subroutine
            
            if (spherical == 1) {
#pragma gpu box(tileBox)
                burner_loop_sphr(AMREX_INT_ANYD(tileBox.loVect()), 
                    AMREX_INT_ANYD(tileBox.hiVect()),
                    BL_TO_FORTRAN_ANYD(s_in_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(s_out_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(source_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(p0_cart_mf[mfi]), dt_in, time_in,
                    BL_TO_FORTRAN_ANYD(mask[mfi]), use_mask);
            } else {
#pragma gpu box(tileBox)
                burner_loop(AMREX_INT_ANYD(tileBox.loVect()), 
                    AMREX_INT_ANYD(tileBox.hiVect()),
                    lev,
                    BL_TO_FORTRAN_ANYD(s_in_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(s_out_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(source_mf[mfi]), 
                    p0.dataPtr(), dt_in, time_in,
                    BL_TO_FORTRAN_ANYD(mask[mfi]), use_mask);
            }
        }
    }
}
#endif



#include <Maestro.H>

using namespace amrex;

void Maestro::Burner(const Vector<MultiFab>& s_in, Vector<MultiFab>& s_out,
                     const Vector<MultiFab>& rho_Hext,
                     Vector<MultiFab>& rho_omegadot, Vector<MultiFab>& rho_Hnuc,
                     [[maybe_unused]] const BaseState<Real>& p0, const Real dt_in,
                     [[maybe_unused]] const Real time_in) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Burner()", Burner);

    // Put tempbar_init on cart
    Vector<MultiFab> tempbar_init_cart(finest_level + 1);

    if (spherical) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            tempbar_init_cart[lev].define(grids[lev], dmap[lev], 1, 0);
            tempbar_init_cart[lev].setVal(0.);
        }

        if (drive_initial_convection) {
            Put1dArrayOnCart(tempbar_init, tempbar_init_cart, false, false,
                             bcs_f, 0);
        }
    }

    const auto ispec_threshold = network_spec_index(burner_threshold_species);

    for (int lev = 0; lev <= finest_level; ++lev) {
        // create mask assuming refinement ratio = 2
        int finelev = lev + 1;
        if (lev == finest_level) {
            finelev = finest_level;
        }

        const BoxArray& fba = s_in[finelev].boxArray();
        const iMultiFab& mask = makeFineMask(s_in[lev], fba, IntVect(2));

        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(s_in[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const bool use_mask = (lev != finest_level);

            const Array4<const Real> s_in_arr = s_in[lev].array(mfi);
            const Array4<Real> s_out_arr = s_out[lev].array(mfi);
            const Array4<const Real> rho_Hext_arr = rho_Hext[lev].array(mfi);
            const Array4<Real> rho_omegadot_arr = rho_omegadot[lev].array(mfi);
            const Array4<Real> rho_Hnuc_arr = rho_Hnuc[lev].array(mfi);
            const auto tempbar_init_arr = tempbar_init.const_array();

            // use a dummy value for non-spherical as tempbar_init_cart is not defined
            const Array4<const Real> tempbar_cart_arr =
                spherical ? tempbar_init_cart[lev].array(mfi)
                          : rho_Hext[lev].array(mfi);
            const Array4<const int> mask_arr = mask.array(mfi);


            reduce_op.eval(tileBox, reduce_data,
            [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) -> ReduceTuple
            {
                if (use_mask && mask_arr(i, j, k) == 1) {
                    return {0.0};  // cell is covered by finer cells
                }
                auto rho = s_in_arr(i, j, k, Rho);
                Real x_in[NumSpec];
                for (int n = 0; n < NumSpec; ++n) {
                    x_in[n] = s_in_arr(i, j, k, FirstSpec + n) / rho;
                }
#if NAUX_NET > 0
                Real aux_in[NumAux];
                for (int n = 0; n < NumAux; ++n) {
                    aux_in[n] = s_in_arr(i, j, k, FirstAux + n) / rho;
                }
#endif

                Real T_in = 0.0;
                if (drive_initial_convection) {
                    if (!spherical) {
                        auto r = (AMREX_SPACEDIM == 2) ? j : k;
                        T_in = tempbar_init_arr(lev, r);
                    } else {
                        T_in = tempbar_cart_arr(i, j, k);
                    }
                } else {
                    T_in = s_in_arr(i, j, k, Temp);
                }

                Real x_test =
                    (ispec_threshold > 0) ? x_in[ispec_threshold] : 0.0;

                burn_t state_in;
                burn_t state_out;


                Real x_out[NumSpec];
#if NAUX_NET > 0
                Real aux_out[NumAux];
#endif
                Real rhowdot[NumSpec];
                Real rhoH = 0.0;

                Real burn_failed = 0.0_rt;

                // if the threshold species is not in the network, then we burn
                // normally.  if it is in the network, make sure the mass
                // fraction is above the cutoff.
                if ((rho > burning_cutoff_density_lo &&
                     rho < burning_cutoff_density_hi) &&
                    (ispec_threshold < 0 ||
                     (ispec_threshold > 0 && x_test > burner_threshold_cutoff))) {

                    // Initialize burn state_in and state_out
                    state_in.e = 0.0;
                    state_in.rho = rho;
                    state_in.T = T_in;
                    for (int n = 0; n < NumSpec; ++n) {
                        state_in.xn[n] = x_in[n];
                    }
#if NAUX_NET > 0
                    for (int n = 0; n < NumAux; ++n) {
                        state_in.aux[n] = aux_in[n];
                    }
#endif

                    // initialize state_out the same as state_in
                    state_out.success = true;

                    state_out.e = 0.0;
                    state_out.rho = rho;
                    state_out.T = T_in;
                    for (int n = 0; n < NumSpec; ++n) {
                        state_out.xn[n] = x_in[n];
                    }
#if NAUX_NET > 0
                    for (int n = 0; n < NumAux; ++n) {
                        state_out.aux[n] = aux_in[n];
                    }
#endif

                    state_out.i = i;
                    state_out.j = j;
                    state_out.k = k;

                    burner(state_out, dt_in);

                    // if we were unsuccessful, update the failure count
                    if (!state_out.success) {
                        burn_failed = 1.0;
                    }

                    for (int n = 0; n < NumSpec; ++n) {
                        x_out[n] = state_out.xn[n];
                        rhowdot[n] = state_out.rho *
                                     (state_out.xn[n] - state_in.xn[n]) / dt_in;
                    }
#if NAUX_NET > 0
                    for (int n = 0; n < NumAux; ++n) {
                        aux_out[n] = state_out.aux[n];
                    }
#endif
                    rhoH = state_out.rho * (state_out.e - state_in.e) / dt_in;
                } else {
                    for (int n = 0; n < NumSpec; ++n) {
                        x_out[n] = x_in[n];
                        rhowdot[n] = 0.0;
                    }
#if NAUX_NET > 0
                    for (int n = 0; n < NumAux; ++n) {
                        aux_out[n] = aux_in[n];
                    }
#endif
                }

                // check if sum{X_k} = 1
                Real sumX = 0.0;
                for (auto X : x_out) {
                    sumX += X;
                }

                if (std::abs(sumX - 1.0) > reaction_sum_tol) {
#ifndef AMREX_USE_GPU
                    std::cout << amrex::Font::Bold << amrex::FGColor::Green
                              << "ERROR: abundances do not sum to 1" << amrex::ResetDisplay << std::endl;
#endif
                    burn_failed = 1.0_rt;
                }

                // pass the density and pi through
                s_out_arr(i, j, k, Rho) = s_in_arr(i, j, k, Rho);
                s_out_arr(i, j, k, Pi) = s_in_arr(i, j, k, Pi);

                // update the species
                for (int n = 0; n < NumSpec; ++n) {
                    s_out_arr(i, j, k, FirstSpec + n) = x_out[n] * rho;
                }

                // update the auxiliary variables
#if NAUX_NET > 0
                for (int n = 0; n < NumAux; ++n) {
                    s_out_arr(i, j, k, FirstAux + n) = aux_out[n] * rho;
                }
#endif

                // store the energy generation and species create quantities
                for (int n = 0; n < NumSpec; ++n) {
                    rho_omegadot_arr(i, j, k, n) = rhowdot[n];
                }
                rho_Hnuc_arr(i, j, k) = rhoH;

                // update the enthalpy -- include the change due to external heating
                s_out_arr(i, j, k, RhoH) = s_in_arr(i, j, k, RhoH) +
                                           dt_in * rho_Hnuc_arr(i, j, k) +
                                           dt_in * rho_Hext_arr(i, j, k);

                return {burn_failed};
            });
        }

        ReduceTuple hv = reduce_data.value();
        Real burn_failed = amrex::get<0>(hv);

        if (burn_failed != 0.0) {
            amrex::Abort("burning failed");
        }
    }
}


#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

#ifndef SDC

#ifndef ML
void Maestro::Burner(const Vector<MultiFab>& s_in, Vector<MultiFab>& s_out,
                     const Vector<MultiFab>& rho_Hext,
                     Vector<MultiFab>& rho_omegadot, Vector<MultiFab>& rho_Hnuc,
                     const BaseState<Real>& p0, const Real dt_in,
                     const Real time_in) {
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

    Vector<MultiFab> react_in(finest_level+1), react_out(finest_level+1);
    const bool save_react_data = istep > 0 && save_react_int > 0 && istep % save_react_int == 0;

    for (int lev = 0; lev <= finest_level; ++lev) {
        if (save_react_data) {
	    if (NumSpec == 3) {
		react_in[lev].define(grids[lev], dmap[lev], NumSpec+1, 0);  // X(C12), X(Mg24), rho, T -> burn
		react_out[lev].define(grids[lev], dmap[lev], 2*(NumSpec+0), 0); // burn -> X(C12), X(Mg24), Hnuc, dXdt, enucdot
	    } else {
		react_in[lev].define(grids[lev], dmap[lev], NumSpec+2, 0);  // X, rho, T -> burn
		react_out[lev].define(grids[lev], dmap[lev], 2*(NumSpec+1), 0); // burn -> X, Hnuc, dXdt, enucdot
	    }
	}

        // create mask assuming refinement ratio = 2
        int finelev = lev + 1;
        if (lev == finest_level) finelev = finest_level;

        const BoxArray& fba = s_in[finelev].boxArray();
        const iMultiFab& mask = makeFineMask(s_in[lev], fba, IntVect(2));

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

            Array4<Real> react_in_arr, react_out_arr;

            if (save_react_data) {
                react_in_arr = react_in[lev].array(mfi);
                react_out_arr = react_out[lev].array(mfi);
            }

            // use a dummy value for non-spherical as tempbar_init_cart is not defined
            const Array4<const Real> tempbar_cart_arr =
                spherical ? tempbar_init_cart[lev].array(mfi)
                          : rho_Hext[lev].array(mfi);
            const Array4<const int> mask_arr = mask.array(mfi);

            // If save_react_data, we add perturbations to the mass fractions to
            // get more varied data for training ML model
            if (save_react_data) {
                ParallelForRNG(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k,
                                                             amrex::RandomEngine const& engine) noexcept
                {
                    if (use_mask && mask_arr(i, j, k))
                        return;  // cell is covered by finer cells

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

		    // Perturb 75 percent of the data
		    Real rand_perc = 0.75;
		    if (amrex::Random(engine) < rand_perc) {
			// Loop over the species
			for (int n = 0; n < NumSpec; ++n) {
			    // aprox13 network (NumSpec == 13)
			    // Three main species are C12(1), O16(2), Mg24(4)
			    // Do not perturb the main species
			    if (NumSpec == 13 && n != 1 && n != 2 && n != 4) {
				// Set a random perturbation based on log(X_k)
				Real X_log = std::log10(x_in[n]);
				// We only want to perturb mass fractions below 10^-4
				// if (X_log <= -4.0) {
				    // Random number generated between [-X_log/15, X_log/15]
				    Real rand = (amrex::Random(engine) - 0.5) * 2*X_log/15.0;
				    Real perturb = pow(10.0, X_log + rand) - x_in[n];
				    x_in[n] += perturb;

				    // Pick the larger mass fraction of C12 and Mg24
				    // and subtract perturb since sum(X_k) = 1 
				    int spec_rand = (x_in[1] > x_in[4]) ? 1 : 4;
				    x_in[spec_rand] -= perturb;
				// }
			    }
			    // ignition_simple network (NumSpec == 3)
			    else if (NumSpec == 3 && n != 1) {
				int spec_other = (n == 0) ? 2 : 0;
				Real x_less = (x_in[n] < x_in[spec_other]) ? x_in[n] : x_in[spec_other];
				// Set a random perturbation within -/+ 10% of smaller mass fraction
				Real perturb = (amrex::Random(engine) - 0.5) * 0.1 * x_less;
				// avoid negative mass fractions
				if (x_in[n] < -perturb || x_in[spec_other] < perturb) {
				    perturb = 0.0; 
				}
				x_in[n] += perturb;

				// Subtract perturb from the other main species
				x_in[spec_other] -= perturb;
			    }
			}
			
			// Perturb density and temperature by small relative values
			// since they do not change much
			Real perturb_rel = (amrex::Random(engine) - 0.5) * 2.e-4;
			rho *= (1.0 + perturb_rel);
			perturb_rel = (amrex::Random(engine) - 0.5) * 6.e-4;
			T_in *= (1.0 + perturb_rel);
		    }

		    if (NumSpec == 3) {
			react_in_arr(i, j, k, 0) = x_in[0];
			react_in_arr(i, j, k, 1) = x_in[2];
			react_in_arr(i, j, k, 2) = rho;
			react_in_arr(i, j, k, 3) = T_in;

			react_out_arr(i, j, k, 2) = 0.0; // output generated energy
			for (int n = NumSpec; n < 2*(NumSpec); ++n) {
			    react_out_arr(i, j, k, n) = 0.0;
			}
		    } else {
			for (int n = 0; n < NumSpec; ++n) {
			    react_in_arr(i, j, k, n) = x_in[n];
			    react_out_arr(i, j, k, n) = x_in[n];
			}
			react_in_arr(i, j, k, NumSpec) = rho;    // input density
			react_in_arr(i, j, k, NumSpec+1) = T_in; // input temperature
			react_out_arr(i, j, k, NumSpec) = 0.0;   // output generated energy

			for (int n = NumSpec+1; n < 2*(NumSpec+1); ++n) {
			    react_out_arr(i, j, k, n) = 0.0;
			}
		    }
		    

                    burn_t state_in;
                    burn_t state_out;

                    Real x_out[NumSpec];
#if NAUX_NET > 0
                    Real aux_out[NumAux];
#endif
                    Real rhowdot[NumSpec];
                    Real rhoH = 0.0;

                    // if the threshold species is not in the network, then we burn
                    // normally.  if it is in the network, make sure the mass
                    // fraction is above the cutoff.
                    if ((rho > burning_cutoff_density_lo &&
                         rho < burning_cutoff_density_hi) &&
                        (ispec_threshold < 0 ||
                         (ispec_threshold > 0 &&
                          x_test > burner_threshold_cutoff))) {
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

                        burner(state_out, dt_in);

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

                        // use state_in and state_out to set the reaction outputs and RHS
                        // because the EOS calls require the actual internal energy e,
                        // this should come after we're otherwise finished with state in/out.
			if (NumSpec == 3) {
			    react_out_arr(i, j, k, 0) = state_out.xn[0];
			    react_out_arr(i, j, k, 1) = state_out.xn[2];
			    react_out_arr(i, j, k, 2) = state_out.e - state_in.e;
			} else {
			    for (int n = 0; n < NumSpec; ++n) {
				react_out_arr(i, j, k, n) = state_out.xn[n];
			    }
			    react_out_arr(i, j, k, NumSpec) = state_out.e - state_in.e;
			}
			
                        eos(eos_input_rt, state_in); // get initial e
                        state_out.e += state_in.e;
                        eos(eos_input_re, state_out); // using final e, get final T
                        Array1D<Real, 1, neqs> ydot; // this will hold the RHS
                        actual_rhs(state_out, ydot); // evaluate the RHS to get dYdt

			if (NumSpec == 3) {
			    react_out_arr(i, j, k, NumSpec) = ydot(1) * aion[0]; // save dXdt
			    react_out_arr(i, j, k, NumSpec+1) = ydot(3) * aion[2];
			    react_out_arr(i, j, k, 2*NumSpec-1) = ydot(net_ienuc); // save enucdot
			} else {
			    for (int n = 1; n <= NumSpec; ++n) {
				react_out_arr(i, j, k, NumSpec+n) = ydot(n) * aion[n-1]; // save dXdt
			    }
			    react_out_arr(i, j, k, 2*NumSpec+1) = ydot(net_ienuc); // save enucdot
			}
		    }
                });
            }

            // solve original input states
            ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                if (use_mask && mask_arr(i, j, k))
                    return;  // cell is covered by finer cells

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

                // if the threshold species is not in the network, then we burn
                // normally.  if it is in the network, make sure the mass
                // fraction is above the cutoff.
                if ((rho > burning_cutoff_density_lo &&
                     rho < burning_cutoff_density_hi) &&
                    (ispec_threshold < 0 ||
                     (ispec_threshold > 0 &&
                      x_test > burner_threshold_cutoff))) {
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

                    burner(state_out, dt_in);

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
                for (int n = 0; n < NumSpec; ++n) {
                    sumX += x_out[n];
                }

                if (fabs(sumX - 1.0) > reaction_sum_tol) {
#ifndef AMREX_USE_CUDA
                    Abort("ERROR: abundances do not sum to 1");
#endif
                    for (int n = 0; n < NumSpec; ++n) {
                        state_out.xn[n] /= sumX;
                    }
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
            });
        }
    }

    if (save_react_data) {
	Vector<std::string> react_in_varnames, react_out_varnames;
	for (int i = 0; i < NumSpec; ++i) {
	    std::string spec_string = "X(";
	    spec_string += short_spec_names_cxx[i];
	    spec_string += ')';
	    if (NumSpec == 3 && i != 1) {
		react_in_varnames.push_back(spec_string);
		react_out_varnames.push_back(spec_string);
	    }
	}
	react_in_varnames.push_back("rho");
	react_in_varnames.push_back("temp");
	react_out_varnames.push_back("enuc");
	for (int i = 0; i < NumSpec; ++i) {
	    std::string spec_string = "Xdot(";
	    spec_string += short_spec_names_cxx[i];
	    spec_string += ')';
	    if (NumSpec == 3 && i != 1) {
		react_out_varnames.push_back(spec_string);
	    }
	}
	react_out_varnames.push_back("enucdot");
    
        std::string react_in_name = "react_inputs_";
        std::string react_out_name = "react_outputs_";
        PlotFileName(istep, &react_in_name);
        PlotFileName(istep, &react_out_name);
        react_in_name = react_in_name + "_" + algo_step;
        react_out_name = react_out_name + "_" + algo_step;
        Vector<int> step_array(finest_level+1, istep);
        WriteMultiLevelPlotfile(react_in_name, finest_level + 1, GetVecOfConstPtrs(react_in), react_in_varnames,
                                Geom(), dt_in, step_array, refRatio());
        WriteMultiLevelPlotfile(react_out_name, finest_level + 1, GetVecOfConstPtrs(react_out), react_out_varnames,
                                Geom(), dt_in, step_array, refRatio());
    }
}

#else
// ifdef ML
void Maestro::Burner(const Vector<MultiFab>& s_in, Vector<MultiFab>& s_out,
                     const Vector<MultiFab>& rho_Hext,
                     Vector<MultiFab>& rho_omegadot, Vector<MultiFab>& rho_Hnuc,
                     const BaseState<Real>& p0, const Real dt_in,
                     const Real time_in) {
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

    const bool use_ml_const = use_ml;

    // match data type of pytorch tensor to multifab data
    auto dtype0 = torch::kFloat64;

    Vector<MultiFab> react_in(finest_level+1), react_out(finest_level+1);

    const bool save_react_data = istep > 0 && save_react_int > 0 && istep % save_react_int == 0;

    for (int lev = 0; lev <= finest_level; ++lev) {
        if (save_react_data) {
	    if (NumSpec == 3) {
		react_in[lev].define(grids[lev], dmap[lev], NumSpec+1, 0);  // X(C12), X(Mg24), rho, T -> burn
		react_out[lev].define(grids[lev], dmap[lev], 2*(NumSpec+0), 0); // burn -> X(C12), X(Mg24), Hnuc, dXdt, enucdot
	    } else {
		react_in[lev].define(grids[lev], dmap[lev], NumSpec+2, 0);  // X, rho, T -> burn
		react_out[lev].define(grids[lev], dmap[lev], 2*(NumSpec+1), 0); // burn -> X, Hnuc, dXdt, enucdot
	    }
	}

        // create mask assuming refinement ratio = 2
        int finelev = lev + 1;
        if (lev == finest_level) finelev = finest_level;

        const BoxArray& fba = s_in[finelev].boxArray();
        const iMultiFab& mask = makeFineMask(s_in[lev], fba, IntVect(2));

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

            Array4<Real> react_in_arr, react_out_arr;

            if (save_react_data) {
                react_in_arr = react_in[lev].array(mfi);
                react_out_arr = react_out[lev].array(mfi);
            }

            // use a dummy value for non-spherical as tempbar_init_cart is not defined
            const Array4<const Real> tempbar_cart_arr =
                spherical ? tempbar_init_cart[lev].array(mfi)
                          : rho_Hext[lev].array(mfi);
            const Array4<const int> mask_arr = mask.array(mfi);

	    // retrieve smallend and size of box
            const IntVect bx_lo = tileBox.smallEnd();
            const IntVect nbox = tileBox.size();

            // compute total cells in the box
            int ncell = AMREX_SPACEDIM == 2 ?
                nbox[0] * nbox[1] : nbox[0] * nbox[1] * nbox[2];

            // copy multifabs to input array
            // index ordering: (species, rho, temp)
	    int nextra = (NumSpec == 3) ? 1 : 2;
            const int NumInput = NumSpec + nextra;
            amrex::Gpu::ManagedVector<Real> state_temp(ncell*NumInput);
            Real* AMREX_RESTRICT temp_ptr = state_temp.dataPtr();

	    // copy input multifab to torch tensor
            amrex::ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                int ii = i - bx_lo[0];
                int jj = j - bx_lo[1];
                int index = jj*nbox[0] + ii;
#if AMREX_SPACEDIM == 3
                int kk = k - bx_lo[2];
                index += kk*nbox[0]*nbox[1];
#endif

                auto rho = s_in_arr(i, j, k, Rho);
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

                // array order is row-based [index][comp]
		if (NumSpec == 3) {
		    temp_ptr[index*NumInput] = s_in_arr(i, j, k, FirstSpec) / rho;
		    temp_ptr[index*NumInput + 1] = s_in_arr(i, j, k, FirstSpec + 2) / rho;
		    temp_ptr[index*NumInput + 2] = rho / dens_fac;
		    temp_ptr[index*NumInput + 3] = T_in / temp_fac;
		} else {
		    for (int n = 0; n < NumSpec; ++n) {
			temp_ptr[index*NumInput + n] = s_in_arr(i, j, k, FirstSpec + n) / rho;
		    }
		    temp_ptr[index*NumInput + NumSpec] = rho / dens_fac;
		    temp_ptr[index*NumInput + NumSpec + 1] = T_in / temp_fac;
		}
            });

	    at::Tensor outputs_torch = torch::zeros({ncell, NumSpec+1}, torch::TensorOptions().dtype(dtype0));
#ifdef AMREX_USE_CUDA
            outputs_torch.to(torch::kCUDA);
#endif

            if (use_ml) {
                // create torch tensor from array
#ifndef AMREX_USE_CUDA
                at::Tensor inputs_torch = torch::from_blob(temp_ptr, {ncell, NumInput},
                                                           torch::TensorOptions().dtype(dtype0));
#else
                at::Tensor inputs_torch = torch::from_blob(temp_ptr, {ncell, NumInput},
                                                           torch::TensorOptions().dtype(dtype0).device(torch::kCUDA));
#endif

                // evaluate torch model
                inputs_torch = inputs_torch.to(torch::kFloat32);
                outputs_torch = module.forward({inputs_torch}).toTensor();
                outputs_torch = outputs_torch.to(dtype0);

                if (/*DEBUG=*/0) {
                    Print() << "example input: \n"
                            << inputs_torch.slice(/*dim=*/0, /*start=*/0, /*end=*/5) << '\n';
                    Print() << "example output: \n"
                            << outputs_torch.slice(/*dim=*/0, /*start=*/0, /*end=*/5) << '\n';
                }
            }

            // get accessor to tensor (read-only)
#ifndef AMREX_USE_CUDA
            auto outputs_torch_acc = outputs_torch.accessor<Real,2>();
#else
            auto outputs_torch_acc = outputs_torch.packed_accessor64<Real,2>();
#endif
            auto use_ml_const_istep_gt_0 = use_ml_const && istep > 0; // cannot use use_ml_const and istep directly in CUDA kernel below.
	    
            ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                if (use_mask && mask_arr(i, j, k))
                    return;  // cell is covered by finer cells
		
		int ii = i - bx_lo[0];
                int jj = j - bx_lo[1];
                int index = jj*nbox[0] + ii;
#if AMREX_SPACEDIM == 3
                int kk = k - bx_lo[2];
                index += kk*nbox[0]*nbox[1];
#endif

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

		Real x_out[NumSpec];
#if NAUX_NET > 0
		Real aux_out[NumAux];
#endif
		Real rhowdot[NumSpec];
		Real rhoH = 0.0;

		if (save_react_data) {
		    if (NumSpec == 3) {
			react_in_arr(i, j, k, 0) = x_in[0];
			react_in_arr(i, j, k, 1) = x_in[2];
			react_in_arr(i, j, k, 2) = rho;
			react_in_arr(i, j, k, 3) = T_in;

			react_out_arr(i, j, k, 2) = 0.0; // output generated energy
			for (int n = NumSpec; n < 2*(NumSpec); ++n) {
			    react_out_arr(i, j, k, n) = 0.0;
			}
		    } else {
			for (int n = 0; n < NumSpec; ++n) {
			    react_in_arr(i, j, k, n) = x_in[n];
			    react_out_arr(i, j, k, n) = x_in[n];
			}
			react_in_arr(i, j, k, NumSpec) = rho;    // input density
			react_in_arr(i, j, k, NumSpec+1) = T_in; // input temperature
			react_out_arr(i, j, k, NumSpec) = 0.0;   // output generated energy
			
			for (int n = NumSpec+1; n < 2*(NumSpec+1); ++n) {
			    react_out_arr(i, j, k, n) = 0.0;
			}
		    }
		}

                // if the threshold species is not in the network, then we burn
                // normally.  if it is in the network, make sure the mass
                // fraction is above the cutoff.
                if ((rho > burning_cutoff_density_lo &&
                     rho < burning_cutoff_density_hi) &&
                    (ispec_threshold < 0 ||
                     (ispec_threshold > 0 &&
                      x_test > burner_threshold_cutoff))) {

		    if (use_ml_const_istep_gt_0) {
                        // copy output tensor to multifabs
                        // index ordering: (species, enuc)
			// check if X_k >= 0
			x_out[0] = (outputs_torch_acc[index][0] >= 0.0) ?
				    outputs_torch_acc[index][0] : 0.0;
			// X(O16) does not change
			x_out[1] = x_in[1];
			// check if X_k >= 0
			x_out[2] = (outputs_torch_acc[index][1] >= 0.0) ?
				    outputs_torch_acc[index][1] : 0.0;

                        // note enuc in output tensor is the normalized value
                        // of (state_out.e - state_in.e)
                        rhoH = rho * (outputs_torch_acc[index][2] * enuc_fac) / dt_in;
                    } else {
                        // need to use burner if no ML model was given or at step 0
			burn_t state_in;
			burn_t state_out;

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

			burner(state_out, dt_in);

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
		    }
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
                for (int n = 0; n < NumSpec; ++n) {
                    sumX += x_out[n];
                }

                if (fabs(sumX - 1.0) > reaction_sum_tol) {
		    if (NumSpec == 3) {
			Real sumXmO16 = sumX - x_out[1];
			x_out[0] *= (1.0 - x_out[1]) / sumXmO16;
			x_out[2] *= (1.0 - x_out[1]) / sumXmO16;
		    } else {
			for (int n = 0; n < NumSpec; ++n) {
			    x_out[n] /= sumX;
			}
		    }
                }

		if (use_ml_const) {
		    for (int n = 0; n < NumSpec; ++n) {
			rhowdot[n] = rho * (x_out[n] - x_in[n]) / dt_in;
		    }
		}

		if (save_react_data) {
		    if (NumSpec == 3) {
			react_out_arr(i, j, k, 0) = x_out[0];
			react_out_arr(i, j, k, 1) = x_out[2];
			react_out_arr(i, j, k, 2) = rhoH * dt_in;
			
			react_out_arr(i, j, k, NumSpec) = rhowdot[0]; // save dXdt
			react_out_arr(i, j, k, NumSpec+1) = rhowdot[2];
			react_out_arr(i, j, k, 2*NumSpec-1) = rhoH; // save enucdot
		    } else {
			for (int n = 0; n < NumSpec; ++n) {
			    react_out_arr(i, j, k, n) = x_out[n];
			}
			react_out_arr(i, j, k, NumSpec) = rhoH * dt_in;

			for (int n = 0; n < NumSpec; ++n) {
			    react_out_arr(i, j, k, NumSpec+n) = rhowdot[n]; // save dXdt
			}
			react_out_arr(i, j, k, 2*NumSpec+1) = rhoH; // save enucdot
		    }
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
            });
        }
    }

    if (save_react_data) {
	Vector<std::string> react_in_varnames, react_out_varnames;
	for (int i = 0; i < NumSpec; ++i) {
	    std::string spec_string = "X(";
	    spec_string += short_spec_names_cxx[i];
	    spec_string += ')';
	    if (NumSpec == 3 && i != 1) {
		react_in_varnames.push_back(spec_string);
		react_out_varnames.push_back(spec_string);
	    }
	}
	react_in_varnames.push_back("rho");
	react_in_varnames.push_back("temp");
	react_out_varnames.push_back("enuc");
	for (int i = 0; i < NumSpec; ++i) {
	    std::string spec_string = "Xdot(";
	    spec_string += short_spec_names_cxx[i];
	    spec_string += ')';
	    if (NumSpec == 3 && i != 1) {
		react_out_varnames.push_back(spec_string);
	    }
	}
	react_out_varnames.push_back("enucdot");
    
        std::string react_in_name = "react_inputs_";
        std::string react_out_name = "react_outputs_";
        PlotFileName(istep, &react_in_name);
        PlotFileName(istep, &react_out_name);
        react_in_name = react_in_name + "_" + algo_step;
        react_out_name = react_out_name + "_" + algo_step;
        Vector<int> step_array(finest_level+1, istep);
        WriteMultiLevelPlotfile(react_in_name, finest_level + 1, GetVecOfConstPtrs(react_in), react_in_varnames,
                                Geom(), dt_in, step_array, refRatio());
        WriteMultiLevelPlotfile(react_out_name, finest_level + 1, GetVecOfConstPtrs(react_out), react_out_varnames,
                                Geom(), dt_in, step_array, refRatio());
    }
}
#endif

#else
// SDC burner
void Maestro::Burner(const Vector<MultiFab>& s_in, Vector<MultiFab>& s_out,
                     const BaseState<Real>& p0, const Real dt_in,
                     const Real time_in, const Vector<MultiFab>& source) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::BurnerSDC()", BurnerSDC);

    // Put tempbar_init on cart
    Vector<MultiFab> p0_cart(finest_level + 1);
    const auto ispec_threshold = network_spec_index(burner_threshold_species);

    for (int lev = 0; lev <= finest_level; ++lev) {
        p0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
        p0_cart[lev].setVal(0.);
    }
    if (spherical) {
        Put1dArrayOnCart(p0, p0_cart, false, false, bcs_f, 0);
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
        // create mask assuming refinement ratio = 2
        int finelev = lev + 1;
        if (lev == finest_level) finelev = finest_level;

        const BoxArray& fba = s_in[finelev].boxArray();
        const iMultiFab& mask = makeFineMask(s_in[lev], fba, IntVect(2));

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(s_in[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const bool use_mask = !(lev == finest_level);

            const Array4<const Real> s_in_arr = s_in[lev].array(mfi);
            const Array4<Real> s_out_arr = s_out[lev].array(mfi);
            const Array4<const Real> p0_cart_arr = p0_cart[lev].array(mfi);
            const Array4<const Real> source_arr = source[lev].array(mfi);
            const Array4<const int> mask_arr = mask.array(mfi);

            const auto p0_arr = p0.const_array();

            ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                if (use_mask && mask_arr(i, j, k))
                    return;  // cell is covered by finer cells

                auto r = (AMREX_SPACEDIM == 2) ? j : k;

                Real sdc_rhoX[NumSpec];
                for (int n = 0; n < NumSpec; ++n) {
                    sdc_rhoX[n] = source_arr(i, j, k, FirstSpec + n);
                }
                auto sdc_rhoh = source_arr(i, j, k, RhoH);
                auto sdc_p0 = spherical ? p0_cart_arr(i, j, k) : p0_arr(lev, r);

                auto rho_in = s_in_arr(i, j, k, Rho);
                Real rhoX_in[NumSpec];
                for (int n = 0; n < NumSpec; ++n) {
                    rhoX_in[n] = s_in_arr(i, j, k, FirstSpec + n);
                }
#if NAUX_NET > 0
                Real aux_in[NumAux];
                for (int n = 0; n < NumAux; ++n) {
                    aux_in[n] = s_in_arr(i, j, k, FirstAux + n);
                }
#endif
                auto rhoh_in = s_in_arr(i, j, k, RhoH);

                Real x_test = (ispec_threshold > 0)
                                  ? rhoX_in[ispec_threshold] / rho_in
                                  : 0.0;

                burn_t state_in;
                burn_t state_out;

                Real rhoX_out[NumSpec];
#if NAUX_NET > 0
                Real aux_out[NumAux];
#endif
                Real rho_out = 0.0;
                Real rhoh_out = 0.0;

                // if the threshold species is not in the network, then we burn
                // normally.  if it is in the network, make sure the mass
                // fraction is above the cutoff.
                if ((rho_in > burning_cutoff_density_lo &&
                     rho_in < burning_cutoff_density_hi) &&
                    (ispec_threshold < 0 ||
                     (ispec_threshold > 0 &&
                      x_test > burner_threshold_cutoff))) {
                    state_in.p0 = sdc_p0;
                    state_in.rho = rho_in;
                    for (int n = 0; n < NumSpec; ++n) {
                        state_in.y[n] = rhoX_in[n];
                    }
                    state_in.y[SENTH] = rhoh_in;
                    for (int n = 0; n < NumSpec; ++n) {
                        state_in.ydot_a[n] = sdc_rhoX[n];
                    }
                    state_in.ydot_a[SENTH] = sdc_rhoh;
#if NAUX_NET > 0
                    for (int n = 0; n < NumAux; ++n) {
                        state_in.aux[n] = rhoaux_in[n] / rho_in;
                    }
#endif
                    state_in.success = true;

                    // initialize state_out the same
                    state_out.p0 = sdc_p0;
                    state_out.rho = rho_in;
                    for (int n = 0; n < NumSpec; ++n) {
                        state_out.y[n] = rhoX_in[n];
                    }
                    state_out.y[NumSpec] = rhoh_in;
                    for (int n = 0; n < NumSpec; ++n) {
                        state_out.ydot_a[n] = sdc_rhoX[n];
                    }
                    state_out.ydot_a[NumSpec] = sdc_rhoh;
#if NAUX_NET > 0
                    for (int n = 0; n < NumAux; ++n) {
                        state_out.aux[n] = rhoaux_in[n] / rho_in;
                    }
#endif
                    state_out.success = true;

                    integrator(state_out, dt_in);

                    for (int n = 0; n < NumSpec; ++n) {
                        rho_out += state_out.y[n];
                        rhoX_out[n] = state_out.y[n];
                    }
                    rhoh_out = state_out.y[NumSpec];
#if NAUX_NET > 0
                    for (int n = 0; n < NumAux; ++n) {
                        rhoaux_out[n] = rho_out * state_out.aux[n];
                    }
#endif
                } else {
                    rho_out = rho_in;
                    for (int n = 0; n < NumSpec; ++n) {
                        rho_out += sdc_rhoX[n] * dt_in;
                        rhoX_out[n] = rhoX_in[n] + sdc_rhoX[n] * dt_in;
                    }
                    rhoh_out = rhoh_in + sdc_rhoh * dt_in;
#if NAUX_NET > 0
                    for (int n = 0; n < NumAux; ++n) {
                        rhoaux_out[n] = rho_out / rho_in * rhoaux_in[n];
                    }
#endif
                }

                // update the density
                s_out_arr(i, j, k, Rho) = rho_out;

                // update the species
                for (int n = 0; n < NumSpec; ++n) {
                    s_out_arr(i, j, k, FirstSpec + n) = rhoX_out[n];
                }
#if NAUX_NET > 0
                for (int n = 0; n < NumAux; ++n) {
                    s_out_arr(i, j, k, FirstAux + n) = rhoaux_out[n];
                }
#endif

                // update the enthalpy -- include the change due to external heating
                s_out_arr(i, j, k, RhoH) = rhoh_out;

                // pass the tracers through (currently not implemented)
            });
        }
    }
}
#endif

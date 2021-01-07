
#include <Maestro.H>
#include <Maestro_F.H>
#include <Problem_F.H>

using namespace amrex;

// advance a single level for a single time step, updates flux registers
void Maestro::AdvanceTimeStep(bool is_initIter) {
    // // timer for profiling
    BL_PROFILE_VAR("Maestro::AdvanceTimeStep()", AdvanceTimeStep);

    Print() << "\nTimestep " << istep << " starts with TIME = " << t_old
            << " DT = " << dt << std::endl
            << std::endl;

    const auto nr_fine = base_geom.nr_fine;
    const auto max_radial_level = base_geom.max_radial_level;

    // vectors store the multilevel 1D states as one very long array
    // these are cell-centered
    BaseState<Real> p0_minus_peosbar(max_radial_level + 1, nr_fine);
    BaseState<Real> w0_force(max_radial_level + 1, nr_fine);
    BaseState<Real> Sbar_old(max_radial_level + 1, nr_fine);
    BaseState<Real> Sbar_new(max_radial_level + 1, nr_fine);
    BaseState<Real> Sbar_nph(max_radial_level + 1, nr_fine);
    BaseState<Real> delta_chi_w0(max_radial_level + 1, nr_fine);
    BaseState<Real> Hext_bar(max_radial_level + 1, nr_fine);
    BaseState<Real> tempbar_new(max_radial_level + 1, nr_fine);
    BaseState<Real> rhoX0_old(max_radial_level + 1, nr_fine, NumSpec);
    BaseState<Real> rhoX0_new(max_radial_level + 1, nr_fine, NumSpec);

    // vectors store the multilevel 1D states as one very long array
    // these are edge-centered
    BaseState<Real> w0_old(max_radial_level + 1, nr_fine + 1);
    BaseState<Real> rho0_predicted_edge(max_radial_level + 1, nr_fine + 1);

    bool is_predictor;

    p0_minus_peosbar.setVal(0.);
    delta_chi_w0.setVal(0.);

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // ! compute initial gamma1bar_old
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    auto rho0_old_arr = rho0_old.array();
    auto p0_old_arr = p0_old.array();
    auto rhoX0_old_arr = rhoX0_old.array();
    auto tempbar_arr = tempbar.array();
    auto gamma1bar_old_arr = gamma1bar_old.array();

    for (auto l = 0; l <= base_geom.max_radial_level; ++l) {
        for (auto r = 0; r < base_geom.nr(l); ++r) {
            eos_t eos_state;

            eos_state.rho = rho0_old_arr(l, r);
            eos_state.p = p0_old_arr(l, r);
            for (auto n = 0; n < NumSpec; ++n) {
                eos_state.xn[n] = rhoX0_old_arr(l, r, n) / rho0_old_arr(l, r);
            }
            eos_state.T = tempbar_arr(l, r);

            eos(eos_input_rp, eos_state);

            gamma1bar_old_arr(l, r) = eos_state.gam1;
        }
    }

    // copy rhoX0_old
    auto s0_init_arr = s0_init.array();

    for (auto n = 0; n < base_geom.max_radial_level; ++n) {
        for (auto r = 0; r < base_geom.nr(n); ++r) {
            for (auto comp = 0; comp < NumSpec; ++comp) {
                rhoX0_old_arr(n, r, comp) = s0_init_arr(n, r, FirstSpec + comp);
            }
        }
    }

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // ! compute the heating term and Sbar
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    auto Hext_bar_arr = Hext_bar.array();

    Hext_bar.setVal(0.0);

    if (prob_type == 1) {
        if (t_old <= heating_time) {
            Real fac;
            if ((t_old + dt) > heating_time) {
                fac = (heating_time - t_old) / dt;
            } else {
                fac = 1.0;
            }

            for (int n = 0; n <= max_radial_level; ++n) {
                for (int r = 0; r < base_geom.nr(n); ++r) {
                    if (!spherical) {
                        // plane-parallel -- do the heating term in paper II (section 4)
                        Hext_bar_arr(n, r) =
                            fac * heating_peak *
                            std::exp(
                                -((base_geom.r_cc_loc(n, r) - heating_rad) *
                                  (base_geom.r_cc_loc(n, r) - heating_rad)) /
                                heating_sigma);
                    } else {
                        // spherical -- lower amplitude heating term
                        Hext_bar_arr(n, r) =
                            fac * heating_peak *
                            std::exp(
                                -((base_geom.r_cc_loc(n, r) - heating_rad) *
                                  (base_geom.r_cc_loc(n, r) - heating_rad)) /
                                heating_sigma);
                    }
                }
            }
        }

    } else if (prob_type == 2) {
        // analytic heating modeling CNO cycle

        const auto h1_comp = network_spec_index("hydrogen-1");
        const auto c12_comp = network_spec_index("carbon-12");
        const auto n14_comp = network_spec_index("nitrogen-14");
        const auto o16_comp = network_spec_index("oxygen-16");

        for (int n = 0; n <= max_radial_level; ++n) {
            for (int r = 0; r < base_geom.nr(n); ++r) {
                Real rho = rho0_old_arr(n, r);
                Real T_6_third =
                    std::pow((tempbar_arr(n, r) / 1.0e6), 1.0 / 3.0);
                Real tmp1 = rhoX0_old_arr(n, r, c12_comp);
                Real tmp2 = rhoX0_old_arr(n, r, n14_comp);
                Real tmp3 = rhoX0_old_arr(n, r, o16_comp);
                Real X_CNO = (tmp1 + tmp2 + tmp3) / rho;
                Real X_1 = rhoX0_old_arr(n, r, h1_comp) / rho;

                tmp1 = 2.7e-3 * T_6_third;
                tmp2 = -7.78e-3 * T_6_third * T_6_third;
                tmp3 = -1.49e-4 * T_6_third * T_6_third * T_6_third;
                Real g14 = 1.0 + tmp1 + tmp2 + tmp3;

                tmp1 =
                    8.67e27 * g14 * X_CNO * X_1 * rho / (T_6_third * T_6_third);
                tmp2 = std::exp(-1.5228e2 / T_6_third);
                Hext_bar_arr(n, r) = tmp1 * tmp2;
            }
        }

    } else if (prob_type == 3) {
        const auto he4_comp = network_spec_index("helium-4");

        // off-center heating for sub_chandra
        if (t_old <= heating_time) {
            Real fac;
            if ((t_old + dt) > heating_time) {
                fac = (heating_time - t_old) / dt;
            } else {
                fac = 1.0;
            }

            for (int n = 0; n <= max_radial_level; ++n) {
                for (int r = 0; r < base_geom.nr(n); ++r) {
                    if (!spherical) {
                        // Abort("ERROR: heating not supported")

                        Hext_bar_arr(n, r) =
                            fac * heating_peak *
                            std::exp(
                                -((base_geom.r_cc_loc(n, r) - heating_rad) *
                                  (base_geom.r_cc_loc(n, r) - heating_rad)) /
                                heating_sigma);
                    } else {
                        // spherical -- lower amplitude heating term
                        Hext_bar_arr(n, r) =
                            fac * heating_peak *
                            std::exp(
                                -((base_geom.r_cc_loc(n, r) - heating_rad) *
                                  (base_geom.r_cc_loc(n, r) - heating_rad)) /
                                (heating_sigma * heating_sigma));

                        // only heat if there is He-4
                        Hext_bar_arr(n, r) *=
                            rhoX0_old_arr(n, r, he4_comp) / rho0_old_arr(n, r);
                    }
                }
            }
        }

    } else if (prob_type == 4) {
        // Apply both heating and cooling for an Urca process

        if (t_old <= heating_time) {
            Real fac;
            if ((t_old + dt) > heating_time) {
                fac = (heating_time - t_old) / dt;
            } else {
                fac = 1.0;
            }

            for (int n = 0; n <= max_radial_level; ++n) {
                for (int r = 0; r < base_geom.nr(n); ++r) {
                    if (!spherical) {
                        // plane-parallel -- do the heating term in paper II (section 4)
                        // plus a similar cooling term for Urca
                        Hext_bar_arr(n, r) =
                            fac * (heating_peak *
                                       std::exp(-((base_geom.r_cc_loc(n, r) -
                                                   heating_rad) *
                                                  (base_geom.r_cc_loc(n, r) -
                                                   heating_rad)) /
                                                heating_sigma) +
                                   cooling_peak *
                                       std::exp(-((base_geom.r_cc_loc(n, r) -
                                                   cooling_rad) *
                                                  (base_geom.r_cc_loc(n, r) -
                                                   cooling_rad)) /
                                                cooling_sigma));

                    } else {
                        // spherical -- lower amplitude heating/cooling term
                        Hext_bar_arr(n, r) =
                            fac * (heating_peak *
                                       std::exp(-((base_geom.r_cc_loc(n, r) -
                                                   heating_rad) *
                                                  (base_geom.r_cc_loc(n, r) -
                                                   heating_rad)) /
                                                heating_sigma) +
                                   cooling_peak *
                                       std::exp(-((base_geom.r_cc_loc(n, r) -
                                                   cooling_rad) *
                                                  (base_geom.r_cc_loc(n, r) -
                                                   cooling_rad)) /
                                                cooling_sigma));
                    }
                }
            }
        }

    } else {
        Print() << "ERROR: " << prob_type << " prob_type not yet supported."
                << std::endl;
        Abort();
    }

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // ! make Sbar
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    auto Sbar_old_arr = Sbar_old.array();

    for (auto l = 0; l <= base_geom.max_radial_level; ++l) {
        for (auto r = 0; r < base_geom.nr(l); ++r) {
            eos_t eos_state;

            eos_state.rho = rho0_old_arr(l, r);
            eos_state.T = tempbar_arr(l, r);
            for (auto n = 0; n < NumSpec; ++n) {
                eos_state.xn[n] = rhoX0_old_arr(l, r, n) / rho0_old_arr(l, r);
            }

            eos(eos_input_rt, eos_state);

            Sbar_old_arr(l, r) =
                Hext_bar_arr(l, r) * eos_state.dpdT /
                (eos_state.rho * eos_state.cp * eos_state.dpdr);
        }
    }

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // ! compute w_0
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    w0_old.copy(w0);

    // compute w0, w0_force, and delta_chi_w0
    is_predictor = true;
    Makew0(w0_old, w0_force, Sbar_old, rho0_old, rho0_old, p0_old, p0_old,
           gamma1bar_old, gamma1bar_old, p0_minus_peosbar, dt, dtold,
           is_predictor);

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // ! update density and compute rho0_predicted_edge
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    AdvectBaseDens(rho0_predicted_edge);

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // ! recompute cutoff coordinates now that rho0 has changed
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    auto rho0_new_arr = rho0_new.array();
    auto p0_new_arr = p0_new.array();
    auto rhoX0_new_arr = rhoX0_new.array();
    auto gamma1bar_new_arr = gamma1bar_new.array();
    auto tempbar_new_arr = tempbar_new.array();

    ComputeCutoffCoords(rho0_new);
    base_geom.ComputeCutoffCoords(rho0_new.array());

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // ! compute gravity
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    MakeGravCell(grav_cell_new, rho0_new);

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // ! update species
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    UpdateSpecies(rho0_old, rho0_predicted_edge, rhoX0_old, rhoX0_new);

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // ! update pressure
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    // set new p0 through HSE
    p0_new.copy(p0_old);

    EnforceHSE(rho0_new, p0_new, grav_cell_new);

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // ! compute gamma1bar_new
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    for (auto l = 0; l <= base_geom.max_radial_level; ++l) {
        for (auto r = 0; r < base_geom.nr(l); ++r) {
            eos_t eos_state;

            eos_state.rho = rho0_new_arr(l, r);
            eos_state.p = p0_new_arr(l, r);
            for (auto n = 0; n < NumSpec; ++n) {
                eos_state.xn[n] = rhoX0_new_arr(l, r, n) / rho0_new_arr(l, r);
            }
            eos_state.T = tempbar_arr(l, r);

            eos(eos_input_rp, eos_state);

            gamma1bar_new_arr(l, r) = eos_state.gam1;
        }
    }

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // ! update temperature
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    for (auto l = 0; l <= base_geom.max_radial_level; ++l) {
        for (auto r = 0; r < base_geom.nr(l); ++r) {
            eos_t eos_state;

            eos_state.rho = rho0_new_arr(l, r);
            eos_state.p = p0_new_arr(l, r);
            for (auto n = 0; n < NumSpec; ++n) {
                eos_state.xn[n] = rhoX0_new_arr(l, r, n) / rho0_new_arr(l, r);
            }
            eos_state.T = tempbar_arr(l, r);

            eos(eos_input_rp, eos_state);

            tempbar_new_arr(l, r) = eos_state.T;
        }
    }

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // ! reset cutoff coordinates
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ComputeCutoffCoords(rho0_old);
    base_geom.ComputeCutoffCoords(rho0_old.array());

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // ! make Sbar
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    auto Sbar_new_arr = Sbar_new.array();

    for (auto l = 0; l <= base_geom.max_radial_level; ++l) {
        for (auto r = 0; r < base_geom.nr(l); ++r) {
            eos_t eos_state;

            eos_state.rho = rho0_new_arr(l, r);
            eos_state.T = tempbar_arr(l, r);
            for (auto n = 0; n < NumSpec; ++n) {
                eos_state.xn[n] = rhoX0_new_arr(l, r, n) / rho0_new_arr(l, r);
            }

            eos(eos_input_rt, eos_state);

            Sbar_new_arr(l, r) =
                Hext_bar_arr(l, r) * eos_state.dpdT /
                (eos_state.rho * eos_state.cp * eos_state.dpdr);
        }
    }

    Sbar_nph.copy(0.5 * (Sbar_old + Sbar_new));

    p0_minus_peosbar.setVal(0.);

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // ! compute w_0
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    w0_old.copy(w0);

    is_predictor = false;
    Makew0(w0_old, w0_force, Sbar_nph, rho0_old, rho0_new, p0_old, p0_new,
           gamma1bar_old, gamma1bar_new, p0_minus_peosbar, dt, dtold,
           is_predictor);

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // ! update density and compute rho0_predicted_edge
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    AdvectBaseDens(rho0_predicted_edge);

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // ! recompute cutoff coordinates now that rho0 has changed
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ComputeCutoffCoords(rho0_new);
    base_geom.ComputeCutoffCoords(rho0_new.array());

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // ! compute gravity
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    MakeGravCell(grav_cell_new, rho0_new);

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // ! update species
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    UpdateSpecies(rho0_old, rho0_predicted_edge, rhoX0_old, rhoX0_new);

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // ! update pressure
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    p0_new.copy(p0_old);

    EnforceHSE(rho0_new, p0_new, grav_cell_new);

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // ! compute gamma1bar_new
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    for (auto l = 0; l <= base_geom.max_radial_level; ++l) {
        for (auto r = 0; r < base_geom.nr(l); ++r) {
            eos_t eos_state;

            eos_state.rho = rho0_new_arr(l, r);
            eos_state.p = p0_new_arr(l, r);
            for (auto n = 0; n < NumSpec; ++n) {
                eos_state.xn[n] = rhoX0_new_arr(l, r, n) / rho0_new_arr(l, r);
            }
            eos_state.T = tempbar_arr(l, r);

            eos(eos_input_rp, eos_state);

            gamma1bar_new_arr(l, r) = eos_state.gam1;
        }
    }

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // ! update temperature
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    for (auto l = 0; l <= base_geom.max_radial_level; ++l) {
        for (auto r = 0; r < base_geom.nr(l); ++r) {
            eos_t eos_state;

            eos_state.rho = rho0_new_arr(l, r);
            eos_state.p = p0_new_arr(l, r);
            for (auto n = 0; n < NumSpec; ++n) {
                eos_state.xn[n] = rhoX0_new_arr(l, r, n) / rho0_new_arr(l, r);
            }
            eos_state.T = tempbar_arr(l, r);

            eos(eos_input_rp, eos_state);

            tempbar_new_arr(l, r) = eos_state.T;
        }
    }

    // rhoX0_old.swap(rhoX0_new);
    tempbar.swap(tempbar_new);

    // copy rhoX0_new into s0_init_arr
    for (auto n = 0; n < base_geom.max_radial_level; ++n) {
        for (auto r = 0; r < base_geom.nr(n); ++r) {
            for (auto comp = 0; comp < NumSpec; ++comp) {
                s0_init_arr(n, r, FirstSpec + comp) = rhoX0_new_arr(n, r, comp);
            }
        }
    }
}

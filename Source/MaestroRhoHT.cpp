
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void Maestro::TfromRhoH(Vector<MultiFab>& scal, const BaseState<Real>& p0) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::TfromRhoH()", TfromRhoH);

    Vector<MultiFab> p0_cart(finest_level + 1);

    for (int lev = 0; lev <= finest_level; ++lev) {
        p0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
        p0_cart[lev].setVal(0.);
    }
    Put1dArrayOnCart(p0, p0_cart, false, false, bcs_f, 0);

    const auto use_eos_e_instead_of_h_loc = use_eos_e_instead_of_h;

    for (int lev = 0; lev <= finest_level; ++lev) {
        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(scal[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const Array4<Real> state = scal[lev].array(mfi);
            const Array4<const Real> p0_arr = p0_cart[lev].array(mfi);

            if (use_eos_e_instead_of_h_loc) {
                // (rho, (h->e)) --> T, p
                ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    eos_t eos_state;

                    eos_state.rho = state(i, j, k, Rho);
                    eos_state.T = state(i, j, k, Temp);
                    for (auto n = 0; n < NumSpec; ++n) {
                        eos_state.xn[n] =
                            state(i, j, k, FirstSpec + n) / eos_state.rho;
                    }
#if NAUX_NET > 0
                    for (auto n = 0; n < NumAux; ++n) {
                        eos_state.aux[n] =
                            state(i, j, k, FirstAux + n) / eos_state.rho;
                    }
#endif

                    // e = (rhoh - p)/rho
                    eos_state.e = (state(i, j, k, RhoH) - p0_arr(i, j, k)) /
                                  state(i, j, k, Rho);

                    eos(eos_input_re, eos_state);

                    state(i, j, k, Temp) = eos_state.T;
                });
            } else {
                // (rho, h) --> T, p
                ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    eos_t eos_state;

                    eos_state.rho = state(i, j, k, Rho);
                    eos_state.T = state(i, j, k, Temp);
                    for (auto n = 0; n < NumSpec; ++n) {
                        eos_state.xn[n] =
                            state(i, j, k, FirstSpec + n) / eos_state.rho;
                    }
#if NAUX_NET > 0
                    for (auto n = 0; n < NumAux; ++n) {
                        eos_state.aux[n] =
                            state(i, j, k, FirstAux + n) / eos_state.rho;
                    }
#endif

                    eos_state.h = state(i, j, k, RhoH) / state(i, j, k, Rho);

                    eos(eos_input_rh, eos_state);

                    state(i, j, k, Temp) = eos_state.T;
                });
            }
        }
    }

    // average down and fill ghost cells
    AverageDown(scal, Temp, 1);
    FillPatch(t_old, scal, scal, scal, Temp, Temp, 1, Temp, bcs_s);
}

void Maestro::TfromRhoP(Vector<MultiFab>& scal, const BaseState<Real>& p0,
                        const bool updateRhoH) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::TfromRhoP()", TfromRhoP);

    Vector<MultiFab> p0_cart(finest_level + 1);

    for (int lev = 0; lev <= finest_level; ++lev) {
        p0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
        p0_cart[lev].setVal(0.);
    }
    Put1dArrayOnCart(p0, p0_cart, false, false, bcs_f, 0);

    const auto use_pprime_in_tfromp_loc = use_pprime_in_tfromp;

    for (int lev = 0; lev <= finest_level; ++lev) {
        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(scal[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const Array4<Real> state = scal[lev].array(mfi);
            const Array4<const Real> p0_arr = p0_cart[lev].array(mfi);

            // (rho, p) --> T
            ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                eos_t eos_state;

                eos_state.rho = state(i, j, k, Rho);
                eos_state.T = state(i, j, k, Temp);

                if (use_pprime_in_tfromp_loc) {
                    eos_state.p = p0_arr(i, j, k) + state(i, j, k, Pi);
                } else {
                    eos_state.p = p0_arr(i, j, k);
                }

                for (auto n = 0; n < NumSpec; ++n) {
                    eos_state.xn[n] =
                        state(i, j, k, FirstSpec + n) / eos_state.rho;
                }
#if NAUX_NET > 0
                for (auto n = 0; n < NumAux; ++n) {
                    eos_state.aux[n] =
                        state(i, j, k, FirstAux + n) / eos_state.rho;
                }
#endif

                eos(eos_input_rp, eos_state);

                state(i, j, k, Temp) = eos_state.T;

                if (updateRhoH) {
                    state(i, j, k, RhoH) = eos_state.rho * eos_state.h;
                }
            });
        }
    }

    // average down and fill ghost cells (Temperature)
    AverageDown(scal, Temp, 1);
    FillPatch(t_old, scal, scal, scal, Temp, Temp, 1, Temp, bcs_s);

    // average down and fill ghost cells (Enthalpy)
    if (updateRhoH) {
        AverageDown(scal, RhoH, 1);
        FillPatch(t_old, scal, scal, scal, RhoH, RhoH, 1, RhoH, bcs_s);
    }
}

void Maestro::PfromRhoH(const Vector<MultiFab>& state,
                        const Vector<MultiFab>& s_old, Vector<MultiFab>& peos) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::PfromRhoH()", PfromRhoH);

    for (int lev = 0; lev <= finest_level; ++lev) {
        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(state[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const Array4<const Real> state_arr = state[lev].array(mfi);
            const Array4<const Real> temp_old = s_old[lev].array(mfi, Temp);
            const Array4<Real> peos_arr = peos[lev].array(mfi);

            // (rho, H) --> T, p
            ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                eos_t eos_state;

                eos_state.rho = state_arr(i, j, k, Rho);
                eos_state.T = temp_old(i, j, k);

                for (auto n = 0; n < NumSpec; ++n) {
                    eos_state.xn[n] =
                        state_arr(i, j, k, FirstSpec + n) / eos_state.rho;
                }
#if NAUX_NET > 0
                for (auto n = 0; n < NumAux; ++n) {
                    eos_state.aux[n] =
                        state_arr(i, j, k, FirstAux + n) / eos_state.rho;
                }
#endif

                eos_state.h =
                    state_arr(i, j, k, RhoH) / state_arr(i, j, k, Rho);

                eos(eos_input_rh, eos_state);

                peos_arr(i, j, k) = eos_state.p;
            });
        }
    }

    // average down and fill ghost cells
    AverageDown(peos, 0, 1);
    FillPatch(t_old, peos, peos, peos, 0, 0, 1, 0, bcs_f);
}

void Maestro::MachfromRhoH(const Vector<MultiFab>& scal,
                           const Vector<MultiFab>& vel,
                           const BaseState<Real>& p0,
                           const Vector<MultiFab>& w0cart,
                           Vector<MultiFab>& mach) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MachfromRhoH()", MachfromRhoH);

    Vector<MultiFab> p0_cart(finest_level + 1);

    for (int lev = 0; lev <= finest_level; ++lev) {
        p0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
        p0_cart[lev].setVal(0.);
    }
    Put1dArrayOnCart(p0, p0_cart, false, false, bcs_f, 0);

    const auto use_eos_e_instead_of_h_loc = use_eos_e_instead_of_h;

    for (int lev = 0; lev <= finest_level; ++lev) {
        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(scal[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const Array4<const Real> state = scal[lev].array(mfi);
            const Array4<const Real> u = vel[lev].array(mfi);
            const Array4<const Real> p0_arr = p0_cart[lev].array(mfi);
            const Array4<const Real> w0_arr = w0cart[lev].array(mfi);
            const Array4<Real> mach_arr = mach[lev].array(mfi);

            ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

            // vel is the magnitude of the velocity, including w0
#if (AMREX_SPACEDIM == 2)
                Real velocity = sqrt(u(i, j, k, 0) * u(i, j, k, 0) +
                                     (u(i, j, k, 1) + w0_arr(i, j, k)) *
                                         (u(i, j, k, 1) + w0_arr(i, j, k)));
#else
                Real velocity = sqrt(u(i,j,k,0)*u(i,j,k,0) + 
                    u(i,j,k,1)*u(i,j,k,1) + 
                    (u(i,j,k,2) + w0_arr(i,j,k))*(u(i,j,k,2) + w0_arr(i,j,k)));
#endif

                eos_t eos_state;

                eos_state.rho = state(i, j, k, Rho);
                eos_state.T = state(i, j, k, Temp);

                for (auto n = 0; n < NumSpec; ++n) {
                    eos_state.xn[n] =
                        state(i, j, k, FirstSpec + n) / eos_state.rho;
                }
#if NAUX_NET > 0
                for (auto n = 0; n < NumAux; ++n) {
                    eos_state.aux[n] =
                        state(i, j, k, FirstAux + n) / eos_state.rho;
                }
#endif

                if (use_eos_e_instead_of_h_loc) {
                    // e = h - p/rho
                    eos_state.e = (state(i, j, k, RhoH) - p0_arr(i, j, k)) /
                                  state(i, j, k, Rho);

                    eos(eos_input_re, eos_state);
                } else {
                    eos_state.h = state(i, j, k, RhoH) / state(i, j, k, Rho);

                    eos(eos_input_rh, eos_state);
                }

                mach_arr(i, j, k) = velocity / eos_state.cs;
            });
        }
    }

    // average down and fill ghost cells
    AverageDown(mach, 0, 1);
    FillPatch(t_old, mach, mach, mach, 0, 0, 1, 0, bcs_f);
}

void Maestro::CsfromRhoH(const Vector<MultiFab>& scal,
                         const Vector<MultiFab>& p0_cart,
                         Vector<MultiFab>& cs) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::CsfromRhoH()", CsfromRhoH);

    const auto use_eos_e_instead_of_h_loc = use_eos_e_instead_of_h;

    for (int lev = 0; lev <= finest_level; ++lev) {
        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(scal[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const Array4<const Real> state = scal[lev].array(mfi);
            const Array4<const Real> p0_arr = p0_cart[lev].array(mfi);
            const Array4<Real> cs_arr = cs[lev].array(mfi);

            ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                eos_t eos_state;

                eos_state.rho = state(i, j, k, Rho);
                eos_state.T = state(i, j, k, Temp);

                for (auto n = 0; n < NumSpec; ++n) {
                    eos_state.xn[n] =
                        state(i, j, k, FirstSpec + n) / eos_state.rho;
                }
#if NAUX_NET > 0
                for (auto n = 0; n < NumAux; ++n) {
                    eos_state.aux[n] =
                        state(i, j, k, FirstAux + n) / eos_state.rho;
                }
#endif

                if (use_eos_e_instead_of_h_loc) {
                    // e = h - p/rho
                    eos_state.e = (state(i, j, k, RhoH) - p0_arr(i, j, k)) /
                                  state(i, j, k, Rho);

                    eos(eos_input_re, eos_state);
                } else {
                    eos_state.h = state(i, j, k, RhoH) / state(i, j, k, Rho);

                    eos(eos_input_rh, eos_state);
                }

                cs_arr(i, j, k) = eos_state.cs;
            });
        }
    }

    // average down and fill ghost cells
    AverageDown(cs, 0, 1);
    FillPatch(t_old, cs, cs, cs, 0, 0, 1, 0, bcs_f);
}

void Maestro::HfromRhoTedge(
    Vector<std::array<MultiFab, AMREX_SPACEDIM> >& sedge,
    const BaseState<Real>& rho0_edge_old, const BaseState<Real>& rhoh0_edge_old,
    const BaseState<Real>& rho0_edge_new,
    const BaseState<Real>& rhoh0_edge_new) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::HfromRhoTedge()", HfromRhoTedge);

    Vector<MultiFab> rho0_cart(finest_level + 1);
    Vector<MultiFab> rhoh0_cart(finest_level + 1);
    Vector<MultiFab> tempbar_cart(finest_level + 1);

    for (int lev = 0; lev <= finest_level; ++lev) {
        rho0_cart[lev].define(grids[lev], dmap[lev], 1, 2);
        rhoh0_cart[lev].define(grids[lev], dmap[lev], 1, 2);
        tempbar_cart[lev].define(grids[lev], dmap[lev], 1, 2);
        rho0_cart[lev].setVal(0.);
        rhoh0_cart[lev].setVal(0.);
        tempbar_cart[lev].setVal(0.);
    }

    BaseState<Real> rho0_halftime(base_geom.max_radial_level + 1,
                                  base_geom.nr_fine);
    BaseState<Real> rhoh0_halftime(base_geom.max_radial_level + 1,
                                   base_geom.nr_fine);

    rho0_halftime.copy(0.5 * (rho0_old + rho0_new));
    rhoh0_halftime.copy(0.5 * (rhoh0_old + rhoh0_new));

    Put1dArrayOnCart(rho0_halftime, rho0_cart, false, false, bcs_s, Rho);
    Put1dArrayOnCart(rhoh0_halftime, rhoh0_cart, false, false, bcs_s, RhoH);
    Put1dArrayOnCart(tempbar, tempbar_cart, false, false, bcs_s, Temp);

    // edge variables
    BaseState<Real> rho0_edge(base_geom.max_radial_level + 1,
                              base_geom.nr_fine + 1);
    BaseState<Real> rhoh0_edge(base_geom.max_radial_level + 1,
                               base_geom.nr_fine + 1);
    BaseState<Real> tempbar_edge(base_geom.max_radial_level + 1,
                                 base_geom.nr_fine + 1);

    if (!spherical) {
        CelltoEdge(tempbar, tempbar_edge);
        rho0_edge.copy(0.5 * (rho0_edge_old + rho0_edge_new));
        rhoh0_edge.copy(0.5 * (rhoh0_edge_old + rhoh0_edge_new));
    }

    Vector<MultiFab> rho0_edge_cart(finest_level + 1);
    Vector<MultiFab> rhoh0_edge_cart(finest_level + 1);
    Vector<MultiFab> tempbar_edge_cart(finest_level + 1);

    for (int lev = 0; lev <= finest_level; ++lev) {
        rho0_edge_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        rhoh0_edge_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        tempbar_edge_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        rho0_edge_cart[lev].setVal(0.);
        rhoh0_edge_cart[lev].setVal(0.);
        tempbar_edge_cart[lev].setVal(0.);
    }

    if (!spherical) {
        Put1dArrayOnCart(rho0_edge, rho0_edge_cart, true, false, bcs_s, Rho);
        Put1dArrayOnCart(rhoh0_edge, rhoh0_edge_cart, true, false, bcs_s, RhoH);
        Put1dArrayOnCart(tempbar_edge, tempbar_edge_cart, true, false, bcs_s,
                         Temp);
    }

    // make a lot of local copies
    const auto enthalpy_pred_type_loc = enthalpy_pred_type;
    const auto species_pred_type_loc = species_pred_type;
    const auto predict_Tprime_then_h_loc = predict_Tprime_then_h;
    const auto predict_rhoprime_and_X_loc = predict_rhoprime_and_X;
    const auto predict_rhoX_loc = predict_rhoX;
    const auto predict_rho_and_X_loc = predict_rho_and_X;
    const auto predict_T_then_h_loc = predict_T_then_h;
    const auto predict_T_then_rhohprime_loc = predict_T_then_rhohprime;
    const auto small_temp_loc = small_temp;

    for (int lev = 0; lev <= finest_level; ++lev) {
        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(sold[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const Box& xbx = amrex::growHi(tileBox, 0, 1);
            const Box& ybx = amrex::growHi(tileBox, 1, 1);
#if (AMREX_SPACEDIM == 3)
            const Box& zbx = amrex::growHi(tileBox, 2, 1);
#endif
            const Array4<const Real> rho0_arr = rho0_cart[lev].array(mfi);
            const Array4<const Real> rhoh0_arr = rhoh0_cart[lev].array(mfi);
            const Array4<const Real> tempbar_arr = tempbar_cart[lev].array(mfi);

            const Array4<Real> sedgex = sedge[lev][0].array(mfi);
            const Array4<Real> sedgey = sedge[lev][1].array(mfi);
#if (AMREX_SPACEDIM == 3)
            const Array4<Real> sedgez = sedge[lev][2].array(mfi);
#endif

            if (!spherical) {
                const Array4<const Real> rho0_edge_arr =
                    rho0_edge_cart[lev].array(mfi);
                const Array4<const Real> rhoh0_edge_arr =
                    rhoh0_edge_cart[lev].array(mfi);
                const Array4<const Real> tempbar_edge_arr =
                    tempbar_edge_cart[lev].array(mfi);
                // x-edge
                ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    eos_t eos_state;

                    // get edge-centered temperature
                    if (enthalpy_pred_type_loc == predict_Tprime_then_h_loc) {
                        eos_state.T = amrex::max(
                            sedgex(i, j, k, Temp) + tempbar_arr(i, j, k),
                            small_temp_loc);
                    } else {
                        eos_state.T =
                            amrex::max(sedgex(i, j, k, Temp), small_temp_loc);
                    }

                    // get edge-centered density and species
                    if (species_pred_type_loc == predict_rhoprime_and_X_loc) {
                        // interface states are rho' and X
                        eos_state.rho =
                            sedgex(i, j, k, Rho) + rho0_arr(i, j, k);

                        for (auto n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = sedgex(i, j, k, FirstSpec + n);
                        }
#if NAUX_NET > 0
                        for (auto n = 0; n < NumAux; ++n) {
                            eos_state.aux[n] = sedgex(i, j, k, FirstAux + n);
                        }
#endif

                    } else if (species_pred_type_loc == predict_rhoX_loc) {
                        // interface states are rho and (rho X)
                        eos_state.rho = sedgex(i, j, k, Rho);

                        for (auto n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] =
                                sedgex(i, j, k, FirstSpec + n) / eos_state.rho;
                        }
#if NAUX_NET > 0
                        for (auto n = 0; n < NumAux; ++n) {
                            eos_state.aux[n] =
                                sedgex(i, j, k, FirstAux + n) / eos_state.rho;
                        }
#endif

                    } else if (species_pred_type_loc == predict_rho_and_X_loc) {
                        // interface states are rho and X
                        eos_state.rho = sedgex(i, j, k, Rho);

                        for (auto n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = sedgex(i, j, k, FirstSpec + n);
                        }
#if NAUX_NET > 0
                        for (auto n = 0; n < NumAux; ++n) {
                            eos_state.aux[n] = sedgex(i, j, k, FirstAux + n);
                        }
#endif
                    }

                    eos(eos_input_rt, eos_state);

                    if (enthalpy_pred_type_loc == predict_T_then_h_loc ||
                        enthalpy_pred_type_loc == predict_Tprime_then_h_loc) {
                        sedgex(i, j, k, RhoH) = eos_state.h;
                    } else if (enthalpy_pred_type_loc ==
                               predict_T_then_rhohprime_loc) {
                        sedgex(i, j, k, RhoH) =
                            eos_state.rho * eos_state.h - rhoh0_arr(i, j, k);
                    }
                });

                // y-edge
                ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    eos_t eos_state;

                    // get edge-centered temperature
                    if (enthalpy_pred_type_loc == predict_Tprime_then_h_loc) {
#if (AMREX_SPACEDIM == 2)
                        eos_state.T = amrex::max(
                            sedgey(i, j, k, Temp) + tempbar_edge_arr(i, j, k),
                            small_temp_loc);
#else
                        eos_state.T = amrex::max(sedgey(i,j,k,Temp) + tempbar_arr(i,j,k), small_temp_loc);
#endif
                    } else {
                        eos_state.T =
                            amrex::max(sedgey(i, j, k, Temp), small_temp_loc);
                    }

                    // get edge-centered density and species
                    if (species_pred_type_loc == predict_rhoprime_and_X_loc) {
                        // interface states are rho' and X
#if (AMREX_SPACEDIM == 2)
                        eos_state.rho =
                            sedgey(i, j, k, Rho) + rho0_edge_arr(i, j, k);
#else
                        eos_state.rho = sedgey(i,j,k,Rho) + rho0_arr(i,j,k);
#endif
                        for (auto n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = sedgey(i, j, k, FirstSpec + n);
                        }
#if NAUX_NET > 0
                        for (auto n = 0; n < NumAux; ++n) {
                            eos_state.aux[n] = sedgey(i, j, k, FirstAux + n);
                        }
#endif

                    } else if (species_pred_type_loc == predict_rhoX_loc) {
                        // interface states are rho and (rho X)
                        eos_state.rho = sedgey(i, j, k, Rho);

                        for (auto n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] =
                                sedgey(i, j, k, FirstSpec + n) / eos_state.rho;
                        }
#if NAUX_NET > 0
                        for (auto n = 0; n < NumAux; ++n) {
                            eos_state.aux[n] =
                                sedgey(i, j, k, FirstAux + n) / eos_state.rho;
                        }
#endif

                    } else if (species_pred_type_loc == predict_rho_and_X_loc) {
                        // interface states are rho and X
                        eos_state.rho = sedgey(i, j, k, Rho);

                        for (auto n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = sedgey(i, j, k, FirstSpec + n);
                        }
#if NAUX_NET > 0
                        for (auto n = 0; n < NumAux; ++n) {
                            eos_state.aux[n] = sedgey(i, j, k, FirstAux + n);
                        }
#endif
                    }

                    eos(eos_input_rt, eos_state);

                    if (enthalpy_pred_type_loc == predict_T_then_h_loc ||
                        enthalpy_pred_type_loc == predict_Tprime_then_h_loc) {
                        sedgey(i, j, k, RhoH) = eos_state.h;
                    } else if (enthalpy_pred_type_loc ==
                               predict_T_then_rhohprime_loc) {
#if (AMREX_SPACEDIM == 2)
                        sedgey(i, j, k, RhoH) = eos_state.rho * eos_state.h -
                                                rhoh0_edge_arr(i, j, k);
#else
                        sedgey(i,j,k,RhoH) = eos_state.rho * eos_state.h - rhoh0_arr(i,j,k);
#endif
                    }
                });

#if (AMREX_SPACEDIM == 3)
                // z-edge
                ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    eos_t eos_state;

                    // get edge-centered temperature
                    if (enthalpy_pred_type_loc == predict_Tprime_then_h_loc) {
                        eos_state.T = amrex::max(
                            sedgez(i, j, k, Temp) + tempbar_edge_arr(i, j, k),
                            small_temp_loc);
                    } else {
                        eos_state.T =
                            amrex::max(sedgez(i, j, k, Temp), small_temp_loc);
                    }

                    // get edge-centered density and species
                    if (species_pred_type_loc == predict_rhoprime_and_X_loc) {
                        // interface states are rho' and X
                        eos_state.rho =
                            sedgez(i, j, k, Rho) + rho0_edge_arr(i, j, k);

                        for (auto n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = sedgez(i, j, k, FirstSpec + n);
                        }
#if NAUX_NET > 0
                        for (auto n = 0; n < NumAux; ++n) {
                            eos_state.aux[n] = sedgez(i, j, k, FirstAux + n);
                        }
#endif

                    } else if (species_pred_type_loc == predict_rhoX_loc) {
                        // interface states are rho and (rho X)
                        eos_state.rho = sedgez(i, j, k, Rho);

                        for (auto n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] =
                                sedgez(i, j, k, FirstSpec + n) / eos_state.rho;
                        }
#if NAUX_NET > 0
                        for (auto n = 0; n < NumAux; ++n) {
                            eos_state.aux[n] =
                                sedgez(i, j, k, FirstAux + n) / eos_state.rho;
                        }
#endif

                    } else if (species_pred_type_loc == predict_rho_and_X_loc) {
                        // interface states are rho and X
                        eos_state.rho = sedgez(i, j, k, Rho);

                        for (auto n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = sedgez(i, j, k, FirstSpec + n);
                        }
#if NAUX_NET > 0
                        for (auto n = 0; n < NumAux; ++n) {
                            eos_state.aux[n] = sedgez(i, j, k, FirstAux + n);
                        }
#endif
                    }

                    eos(eos_input_rt, eos_state);

                    if (enthalpy_pred_type_loc == predict_T_then_h_loc ||
                        enthalpy_pred_type_loc == predict_Tprime_then_h_loc) {
                        sedgez(i, j, k, RhoH) = eos_state.h;
                    } else if (enthalpy_pred_type_loc ==
                               predict_T_then_rhohprime_loc) {
                        sedgez(i, j, k, RhoH) = eos_state.rho * eos_state.h -
                                                rhoh0_edge_arr(i, j, k);
                    }
                });
#endif
            } else {
#if (AMREX_SPACEDIM == 3)
                // x-edge
                ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    eos_t eos_state;

                    // get edge-centered temperature
                    if (enthalpy_pred_type_loc == predict_Tprime_then_h_loc) {
                        Real tempbar_edge_l = 0.5 * (tempbar_arr(i - 1, j, k) +
                                                     tempbar_arr(i, j, k));
                        eos_state.T =
                            amrex::max(sedgex(i, j, k, Temp) + tempbar_edge_l,
                                       small_temp_loc);
                    } else {
                        eos_state.T =
                            amrex::max(sedgex(i, j, k, Temp), small_temp_loc);
                    }

                    // get edge-centered density and species
                    if (species_pred_type_loc == predict_rhoprime_and_X_loc) {
                        // interface states are rho' and X
                        Real rho0_edge_loc =
                            0.5 * (rho0_arr(i - 1, j, k) + rho0_arr(i, j, k));
                        eos_state.rho = sedgex(i, j, k, Rho) + rho0_edge_loc;

                        for (auto n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = sedgex(i, j, k, FirstSpec + n);
                        }
#if NAUX_NET > 0
                        for (auto n = 0; n < NumAux; ++n) {
                            eos_state.aux[n] = sedgex(i, j, k, FirstAux + n);
                        }
#endif
                    } else if (species_pred_type_loc == predict_rhoX_loc) {
                        // interface states are rho and (rho X)
                        eos_state.rho = sedgex(i, j, k, Rho);

                        for (auto n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] =
                                sedgex(i, j, k, FirstSpec + n) / eos_state.rho;
                        }
#if NAUX_NET > 0
                        for (auto n = 0; n < NumAux; ++n) {
                            eos_state.aux[n] =
                                sedgex(i, j, k, FirstAux + n) / eos_state.rho;
                        }
#endif
                    } else if (species_pred_type_loc == predict_rho_and_X_loc) {
                        // interface states are rho and X
                        eos_state.rho = sedgex(i, j, k, Rho);

                        for (auto n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = sedgex(i, j, k, FirstSpec + n);
                        }
#if NAUX_NET > 0
                        for (auto n = 0; n < NumAux; ++n) {
                            eos_state.aux[n] = sedgex(i, j, k, FirstAux + n);
                        }
#endif
                    }

                    eos(eos_input_rt, eos_state);

                    if (enthalpy_pred_type_loc == predict_T_then_h_loc ||
                        enthalpy_pred_type_loc == predict_Tprime_then_h_loc) {
                        sedgex(i, j, k, RhoH) = eos_state.h;
                    } else if (enthalpy_pred_type_loc ==
                               predict_T_then_rhohprime_loc) {
                        Real rhoh0_edge_l =
                            0.5 * (rhoh0_arr(i - 1, j, k) + rhoh0_arr(i, j, k));
                        sedgex(i, j, k, RhoH) =
                            eos_state.rho * eos_state.h - rhoh0_edge_l;
                    }
                });

                // y-edge
                ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    eos_t eos_state;

                    // get edge-centered temperature
                    if (enthalpy_pred_type_loc == predict_Tprime_then_h_loc) {
                        Real tempbar_edge_l = 0.5 * (tempbar_arr(i, j - 1, k) +
                                                     tempbar_arr(i, j, k));
                        eos_state.T =
                            amrex::max(sedgey(i, j, k, Temp) + tempbar_edge_l,
                                       small_temp_loc);
                    } else {
                        eos_state.T =
                            amrex::max(sedgey(i, j, k, Temp), small_temp_loc);
                    }

                    // get edge-centered density and species
                    if (species_pred_type_loc == predict_rhoprime_and_X_loc) {
                        // interface states are rho' and X
                        Real rho0_edge_l =
                            0.5 * (rho0_arr(i, j - 1, k) + rho0_arr(i, j, k));
                        eos_state.rho = sedgey(i, j, k, Rho) + rho0_edge_l;

                        for (auto n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = sedgey(i, j, k, FirstSpec + n);
                        }
#if NAUX_NET > 0
                        for (auto n = 0; n < NumAux; ++n) {
                            eos_state.aux[n] = sedgey(i, j, k, FirstAux + n);
                        }
#endif
                    } else if (species_pred_type_loc == predict_rhoX_loc) {
                        // interface states are rho and (rho X)
                        eos_state.rho = sedgey(i, j, k, Rho);

                        for (auto n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] =
                                sedgey(i, j, k, FirstSpec + n) / eos_state.rho;
                        }
#if NAUX_NET > 0
                        for (auto n = 0; n < NumAux; ++n) {
                            eos_state.aux[n] =
                                sedgey(i, j, k, FirstAux + n) / eos_state.rho;
                        }
#endif
                    } else if (species_pred_type_loc == predict_rho_and_X_loc) {
                        // interface states are rho and X
                        eos_state.rho = sedgey(i, j, k, Rho);

                        for (auto n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = sedgey(i, j, k, FirstSpec + n);
                        }
#if NAUX_NET > 0
                        for (auto n = 0; n < NumAux; ++n) {
                            eos_state.aux[n] = sedgey(i, j, k, FirstAux + n);
                        }
#endif
                    }

                    eos(eos_input_rt, eos_state);

                    if (enthalpy_pred_type_loc == predict_T_then_h_loc ||
                        enthalpy_pred_type_loc == predict_Tprime_then_h_loc) {
                        sedgey(i, j, k, RhoH) = eos_state.h;
                    } else if (enthalpy_pred_type_loc ==
                               predict_T_then_rhohprime_loc) {
                        Real rhoh0_edge_l =
                            0.5 * (rhoh0_arr(i, j - 1, k) + rhoh0_arr(i, j, k));
                        sedgey(i, j, k, RhoH) =
                            eos_state.rho * eos_state.h - rhoh0_edge_l;
                    }
                });

                // z-edge
                ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    eos_t eos_state;

                    // get edge-centered temperature
                    if (enthalpy_pred_type_loc == predict_Tprime_then_h_loc) {
                        Real tempbar_edge_l = 0.5 * (tempbar_arr(i, j, k - 1) +
                                                     tempbar_arr(i, j, k));
                        eos_state.T =
                            amrex::max(sedgez(i, j, k, Temp) + tempbar_edge_l,
                                       small_temp_loc);
                    } else {
                        eos_state.T =
                            amrex::max(sedgez(i, j, k, Temp), small_temp_loc);
                    }

                    // get edge-centered density and species
                    if (species_pred_type_loc == predict_rhoprime_and_X_loc) {
                        // interface states are rho' and X
                        Real rho0_edge_loc =
                            0.5 * (rho0_arr(i, j, k - 1) + rho0_arr(i, j, k));
                        eos_state.rho = sedgez(i, j, k, Rho) + rho0_edge_loc;

                        for (auto n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = sedgez(i, j, k, FirstSpec + n);
                        }
#if NAUX_NET > 0
                        for (auto n = 0; n < NumAux; ++n) {
                            eos_state.aux[n] = sedgez(i, j, k, FirstAux + n);
                        }
#endif
                    } else if (species_pred_type_loc == predict_rhoX_loc) {
                        // interface states are rho and (rho X)
                        eos_state.rho = sedgez(i, j, k, Rho);

                        for (auto n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] =
                                sedgez(i, j, k, FirstSpec + n) / eos_state.rho;
                        }
#if NAUX_NET > 0
                        for (auto n = 0; n < NumAux; ++n) {
                            eos_state.aux[n] =
                                sedgez(i, j, k, FirstAux + n) / eos_state.rho;
                        }
#endif
                    } else if (species_pred_type_loc == predict_rho_and_X_loc) {
                        // interface states are rho and X
                        eos_state.rho = sedgez(i, j, k, Rho);

                        for (auto n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = sedgez(i, j, k, FirstSpec + n);
                        }
#if NAUX_NET > 0
                        for (auto n = 0; n < NumAux; ++n) {
                            eos_state.aux[n] = sedgez(i, j, k, FirstAux + n);
                        }
#endif
                    }

                    eos(eos_input_rt, eos_state);

                    if (enthalpy_pred_type_loc == predict_T_then_h_loc ||
                        enthalpy_pred_type_loc == predict_Tprime_then_h_loc) {
                        sedgez(i, j, k, RhoH) = eos_state.h;
                    } else if (enthalpy_pred_type_loc ==
                               predict_T_then_rhohprime_loc) {
                        Real rhoh0_edge_l =
                            0.5 * (rhoh0_arr(i, j, k - 1) + rhoh0_arr(i, j, k));
                        sedgez(i, j, k, RhoH) =
                            eos_state.rho * eos_state.h - rhoh0_edge_l;
                    }
                });
#endif
            }
        }
    }
}

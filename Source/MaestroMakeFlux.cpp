
#include <Maestro.H>

using namespace amrex;

const int pred_rhoh = 0;
const int pred_rhohprime = 1;
const int pred_h = 2;
const int pred_T_then_rhohprime = 3;
const int pred_T_then_h = 4;
const int pred_hprime = 5;
const int pred_Tprime_then_h = 6;

const int pred_rhoprime_and_X = 1;
const int pred_rhoX = 2;
const int pred_rho_and_X = 3;

void Maestro::MakeRhoXFlux(
    const Vector<MultiFab>& state,
    Vector<std::array<MultiFab, AMREX_SPACEDIM> >& sflux,
    Vector<MultiFab>& etarhoflux,
    Vector<std::array<MultiFab, AMREX_SPACEDIM> >& sedge,
    const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& umac,
    const BaseState<Real>& rho0_old_in,
    const BaseState<Real>& rho0_edge_old_state,
    [[maybe_unused]] const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& r0mac_old,
    const BaseState<Real>& rho0_new_in,
    const BaseState<Real>& rho0_edge_new_state,
    [[maybe_unused]] const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& r0mac_new,
    const BaseState<Real>& rho0_predicted_edge_state, int start_comp,
    int num_comp) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeRhoXFlux()", MakeRhoXFlux);

    const int species_pred_type_loc = species_pred_type;
    const bool use_exact_base_state_loc = use_exact_base_state;
    const bool evolve_base_state_loc = evolve_base_state;

    const auto rho0_old_arr = rho0_old_in.const_array();
    const auto rho0_new_arr = rho0_new_in.const_array();

    const auto rho0_edge_old = rho0_edge_old_state.const_array();
    const auto rho0_edge_new = rho0_edge_new_state.const_array();
    const auto rho0_predicted_edge = rho0_predicted_edge_state.const_array();

    // reset density flux
    for (int lev = 0; lev <= finest_level; ++lev) {
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            sflux[lev][i].setVal(0., Rho, 1, 0);
        }
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
#if (AMREX_SPACEDIM == 3)
        MultiFab rho0mac_edgex, rho0mac_edgey, rho0mac_edgez;

        if (spherical) {
            rho0mac_edgex.define(convert(grids[lev], nodal_flag_x), dmap[lev],
                                 1, 1);
            rho0mac_edgey.define(convert(grids[lev], nodal_flag_y), dmap[lev],
                                 1, 1);
            rho0mac_edgez.define(convert(grids[lev], nodal_flag_z), dmap[lev],
                                 1, 1);
            MultiFab::LinComb(rho0mac_edgex, 0.5, r0mac_old[lev][0], 0, 0.5,
                              r0mac_new[lev][0], 0, 0, 1, 1);
            MultiFab::LinComb(rho0mac_edgey, 0.5, r0mac_old[lev][1], 0, 0.5,
                              r0mac_new[lev][1], 0, 0, 1, 1);
            MultiFab::LinComb(rho0mac_edgez, 0.5, r0mac_old[lev][2], 0, 0.5,
                              r0mac_new[lev][2], 0, 0, 1, 1);
        }
#endif

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(state[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region
            const Box& xbx = mfi.nodaltilebox(0);
            const Box& ybx = mfi.nodaltilebox(1);
#if (AMREX_SPACEDIM == 3)
            const Box& zbx = mfi.nodaltilebox(2);
#endif

            const Array4<Real> sedgex = sedge[lev][0].array(mfi);
            const Array4<Real> sfluxx = sflux[lev][0].array(mfi);
            const Array4<Real> etarhoflux_arr = etarhoflux[lev].array(mfi);
            const Array4<const Real> umacx = umac[lev][0].array(mfi);
            const Array4<Real> sedgey = sedge[lev][1].array(mfi);
            const Array4<Real> sfluxy = sflux[lev][1].array(mfi);
            const Array4<const Real> vmac = umac[lev][1].array(mfi);
#if (AMREX_SPACEDIM == 3)
            const Array4<Real> sedgez = sedge[lev][2].array(mfi);
            const Array4<Real> sfluxz = sflux[lev][2].array(mfi);
            const Array4<const Real> wmac = umac[lev][2].array(mfi);
#endif

            const auto w0_arr = w0.const_array();

#if (AMREX_SPACEDIM == 2)

            // x-direction
            ParallelFor(
                xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    // loop over components cannot be part of the ParallelFor
                    // due to race condition on sflux for the Rho component
                    for (int comp = start_comp; comp < num_comp + start_comp;
                         ++comp) {
                        Real rho0_edge =
                            0.5 * (rho0_old_arr(lev, j) + rho0_new_arr(lev, j));

                        if (species_pred_type_loc == pred_rhoprime_and_X) {
                            // edge states are rho' and X.  To make the (rho X) flux,
                            // we need the edge state of rho0
                            sfluxx(i, j, k, comp) =
                                umacx(i, j, k) *
                                (rho0_edge + sedgex(i, j, k, Rho)) *
                                sedgex(i, j, k, comp);

                        } else if (species_pred_type_loc == pred_rhoX) {
                            // edge states are (rho X)
                            sfluxx(i, j, k, comp) =
                                umacx(i, j, k) * sedgex(i, j, k, comp);

                        } else if (species_pred_type_loc == pred_rho_and_X) {
                            // edge states are rho and X
                            sfluxx(i, j, k, comp) = umacx(i, j, k) *
                                                    sedgex(i, j, k, Rho) *
                                                    sedgex(i, j, k, comp);
                        }

                        // compute the density fluxes by summing the species fluxes
                        sfluxx(i, j, k, Rho) += sfluxx(i, j, k, comp);
                    }
                });

            // y-direction
            ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j,
                                                  int k) noexcept {
                for (int comp = start_comp; comp < num_comp + start_comp;
                     ++comp) {
                    Real rho0_edge =
                        0.5 * (rho0_edge_old(lev, j) + rho0_edge_new(lev, j));

                    if (species_pred_type_loc == pred_rhoprime_and_X) {
                        // edge states are rho' and X.  To make the (rho X) flux,
                        // we need the edge state of rho0
                        sfluxy(i, j, k, comp) =
                            vmac(i, j, k) * (rho0_edge + sedgey(i, j, k, Rho)) *
                            sedgey(i, j, k, comp);

                    } else if (species_pred_type_loc == pred_rhoX) {
                        // edge states are (rho X)
                        sfluxy(i, j, k, comp) =
                            vmac(i, j, k) * sedgey(i, j, k, comp);

                    } else if (species_pred_type_loc == pred_rho_and_X) {
                        // edge state are rho and X
                        sfluxy(i, j, k, comp) = vmac(i, j, k) *
                                                sedgey(i, j, k, Rho) *
                                                sedgey(i, j, k, comp);
                    }

                    if (evolve_base_state_loc && !use_exact_base_state_loc) {
                        if (comp >= FirstSpec &&
                            comp <= FirstSpec + NumSpec - 1) {
                            etarhoflux_arr(i, j, k) += sfluxy(i, j, k, comp);
                        }

                        if (comp == FirstSpec + NumSpec - 1) {
                            etarhoflux_arr(i, j, k) -=
                                w0_arr(lev, j) * rho0_predicted_edge(lev, j);
                        }
                    }

                    // compute the density fluxes by summing the species fluxes
                    sfluxy(i, j, k, Rho) += sfluxy(i, j, k, comp);
                }
            });

#elif (AMREX_SPACEDIM == 3)

            if (!spherical) {
                // x-direction
                ParallelFor(
                    xbx, ybx, zbx,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        // loop over components cannot be part of the ParallelFor
                        // due to race condition on sflux for the Rho component
                        for (int comp = start_comp;
                             comp < num_comp + start_comp; ++comp) {
                            Real rho0_edge = 0.5 * (rho0_old_arr(lev, k) +
                                                    rho0_new_arr(lev, k));

                            if (species_pred_type_loc == pred_rhoprime_and_X) {
                                // edge states are rho' and X.  To make the (rho X) flux,
                                // we need the edge state of rho0
                                sfluxx(i, j, k, comp) =
                                    umacx(i, j, k) *
                                    (rho0_edge + sedgex(i, j, k, Rho)) *
                                    sedgex(i, j, k, comp);

                            } else if (species_pred_type_loc == pred_rhoX) {
                                // edge states are (rho X)
                                sfluxx(i, j, k, comp) =
                                    umacx(i, j, k) * sedgex(i, j, k, comp);

                            } else if (species_pred_type_loc ==
                                       pred_rho_and_X) {
                                // edge states are rho and X
                                sfluxx(i, j, k, comp) = umacx(i, j, k) *
                                                        sedgex(i, j, k, Rho) *
                                                        sedgex(i, j, k, comp);
                            }

                            // compute the density fluxes by summing the species fluxes
                            sfluxx(i, j, k, Rho) += sfluxx(i, j, k, comp);
                        }
                    },

                    // y-direction
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        for (int comp = start_comp;
                             comp < num_comp + start_comp; ++comp) {
                            Real rho0_edge = 0.5 * (rho0_old_arr(lev, k) +
                                                    rho0_new_arr(lev, k));

                            if (species_pred_type_loc == pred_rhoprime_and_X) {
                                // edge states are rho' and X.  To make the (rho X) flux,
                                // we need the edge state of rho0
                                sfluxy(i, j, k, comp) =
                                    vmac(i, j, k) *
                                    (rho0_edge + sedgey(i, j, k, Rho)) *
                                    sedgey(i, j, k, comp);

                            } else if (species_pred_type_loc == pred_rhoX) {
                                // edge states are (rho X)
                                sfluxy(i, j, k, comp) =
                                    vmac(i, j, k) * sedgey(i, j, k, comp);

                            } else if (species_pred_type_loc ==
                                       pred_rho_and_X) {
                                // edge states are rho and X
                                sfluxy(i, j, k, comp) = vmac(i, j, k) *
                                                        sedgey(i, j, k, Rho) *
                                                        sedgey(i, j, k, comp);
                            }

                            // compute the density fluxes by summing the species fluxes
                            sfluxy(i, j, k, Rho) += sfluxy(i, j, k, comp);
                        }
                    },

                    // z-direction
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        for (int comp = start_comp;
                             comp < num_comp + start_comp; ++comp) {
                            Real rho0_edge = 0.5 * (rho0_edge_old(lev, k) +
                                                    rho0_edge_new(lev, k));

                            if (species_pred_type_loc == pred_rhoprime_and_X) {
                                // edge states are rho' and X.  To make the (rho X) flux,
                                // we need the edge state of rho0
                                sfluxz(i, j, k, comp) =
                                    wmac(i, j, k) *
                                    (rho0_edge + sedgez(i, j, k, Rho)) *
                                    sedgez(i, j, k, comp);

                            } else if (species_pred_type_loc == pred_rhoX) {
                                // edge states are (rho X)
                                sfluxz(i, j, k, comp) =
                                    wmac(i, j, k) * sedgez(i, j, k, comp);

                            } else if (species_pred_type_loc ==
                                       pred_rho_and_X) {
                                // edge state are rho and X
                                sfluxz(i, j, k, comp) = wmac(i, j, k) *
                                                        sedgez(i, j, k, Rho) *
                                                        sedgez(i, j, k, comp);
                            }

                            if (evolve_base_state_loc &&
                                !use_exact_base_state_loc) {
                                if (comp >= FirstSpec &&
                                    comp <= FirstSpec + NumSpec - 1) {
                                    etarhoflux_arr(i, j, k) +=
                                        sfluxz(i, j, k, comp);
                                }

                                if (comp == FirstSpec + NumSpec - 1) {
                                    etarhoflux_arr(i, j, k) -=
                                        w0_arr(lev, k) *
                                        rho0_predicted_edge(lev, k);
                                }
                            }

                            // compute the density fluxes by summing the species fluxes
                            sfluxz(i, j, k, Rho) += sfluxz(i, j, k, comp);
                        }
                    });
            } else {
                // spherical case

                const Array4<const Real> rho0_edgex = rho0mac_edgex.array(mfi);
                const Array4<const Real> rho0_edgey = rho0mac_edgey.array(mfi);
                const Array4<const Real> rho0_edgez = rho0mac_edgez.array(mfi);

                // x-direction
                amrex::ParallelFor(
                    xbx, ybx, zbx,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        // loop over components cannot be part of the ParallelFor
                        // due to race condition on sflux for the Rho component
                        for (int comp = start_comp;
                             comp < num_comp + start_comp; ++comp) {
                            if (species_pred_type_loc == pred_rhoprime_and_X) {
                                // edge states are rho' and X.  To make the (rho X) flux,
                                // we need the edge state of rho0
                                sfluxx(i, j, k, comp) = umacx(i, j, k) *
                                                        (rho0_edgex(i, j, k) +
                                                         sedgex(i, j, k, Rho)) *
                                                        sedgex(i, j, k, comp);

                            } else if (species_pred_type_loc == pred_rhoX) {
                                // edge states are (rho X)
                                sfluxx(i, j, k, comp) =
                                    umacx(i, j, k) * sedgex(i, j, k, comp);

                            } else if (species_pred_type_loc ==
                                       pred_rho_and_X) {
                                // edge states are rho and X
                                sfluxx(i, j, k, comp) = umacx(i, j, k) *
                                                        sedgex(i, j, k, Rho) *
                                                        sedgex(i, j, k, comp);
                            }

                            // compute the density fluxes by summing the species fluxes
                            sfluxx(i, j, k, Rho) += sfluxx(i, j, k, comp);
                        }
                    },

                    // y-direction
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        for (int comp = start_comp;
                             comp < num_comp + start_comp; ++comp) {
                            if (species_pred_type_loc == pred_rhoprime_and_X) {
                                // edge states are rho' and X.  To make the (rho X) flux,
                                // we need the edge state of rho0
                                sfluxy(i, j, k, comp) = vmac(i, j, k) *
                                                        (rho0_edgey(i, j, k) +
                                                         sedgey(i, j, k, Rho)) *
                                                        sedgey(i, j, k, comp);

                            } else if (species_pred_type_loc == pred_rhoX) {
                                // edge states are (rho X)
                                sfluxy(i, j, k, comp) =
                                    vmac(i, j, k) * sedgey(i, j, k, comp);

                            } else if (species_pred_type_loc ==
                                       pred_rho_and_X) {
                                // edge states are rho and X
                                sfluxy(i, j, k, comp) = vmac(i, j, k) *
                                                        sedgey(i, j, k, Rho) *
                                                        sedgey(i, j, k, comp);
                            }

                            // compute the density fluxes by summing the species fluxes
                            sfluxy(i, j, k, Rho) += sfluxy(i, j, k, comp);
                        }
                    },

                    // z-direction
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        for (int comp = start_comp;
                             comp < num_comp + start_comp; ++comp) {
                            if (species_pred_type_loc == pred_rhoprime_and_X) {
                                // edge states are rho' and X.  To make the (rho X) flux,
                                // we need the edge state of rho0
                                sfluxz(i, j, k, comp) = wmac(i, j, k) *
                                                        (rho0_edgez(i, j, k) +
                                                         sedgez(i, j, k, Rho)) *
                                                        sedgez(i, j, k, comp);

                            } else if (species_pred_type_loc == pred_rhoX) {
                                // edge states are (rho X)
                                sfluxz(i, j, k, comp) =
                                    wmac(i, j, k) * sedgez(i, j, k, comp);

                            } else if (species_pred_type_loc ==
                                       pred_rho_and_X) {
                                // edge state are rho and X
                                sfluxz(i, j, k, comp) = wmac(i, j, k) *
                                                        sedgez(i, j, k, Rho) *
                                                        sedgez(i, j, k, comp);
                            }

                            // compute the density fluxes by summing the species fluxes
                            sfluxz(i, j, k, Rho) += sfluxz(i, j, k, comp);
                        }
                    });
            }  // end spherical
#endif
        }  // end MFIter loop

        // increment or decrement the flux registers by area and time-weighted fluxes
        // Note that the fluxes need to be scaled by dt and area
        // In this example we are solving s_t = -div(+F)
        // The fluxes contain, e.g., F_{i+1/2,j} = (s*u)_{i+1/2,j}
        // Keep this in mind when considering the different sign convention for updating
        // the flux registers from the coarse or fine grid perspective
        // NOTE: the flux register associated with flux_reg_s[lev] is associated
        // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
        if (reflux_type == 2) {
            // Get the grid size
            const Real* dx = geom[lev].CellSize();
            // NOTE: areas are different in DIM=2 and DIM=3
#if (AMREX_SPACEDIM == 3)
            const Real area[3] = {dx[1] * dx[2], dx[0] * dx[2], dx[0] * dx[1]};
#else
            const Real area[2] = {dx[1], dx[0]};
#endif

            if (flux_reg_s[lev + 1]) {
                for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                    // update the lev+1/lev flux register (index lev+1)
                    flux_reg_s[lev + 1]->CrseInit(sflux[lev][i], i, start_comp,
                                                  start_comp, num_comp,
                                                  -1.0 * dt * area[i]);
                    // also include density flux
                    flux_reg_s[lev + 1]->CrseInit(sflux[lev][i], i, Rho, Rho, 1,
                                                  -1.0 * dt * area[i]);
                }
            }
            if (flux_reg_s[lev]) {
                for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                    // update the lev/lev-1 flux register (index lev)
                    flux_reg_s[lev]->FineAdd(sflux[lev][i], i, start_comp,
                                             start_comp, num_comp,
                                             1.0 * dt * area[i]);
                    // also include density flux
                    flux_reg_s[lev]->FineAdd(sflux[lev][i], i, Rho, Rho, 1,
                                             1.0 * dt * area[i]);
                }
            }

            if (!spherical) {
                // need edge_restrict for etarhoflux
            }
        }
    }  // end loop over levels

    // average down fluxes
    if (reflux_type == 1) {
        AverageDownFaces(sflux);
    }

    // Something analogous to edge_restriction is done in UpdateScal()
}

void Maestro::MakeRhoHFlux(
    const Vector<MultiFab>& state,
    Vector<std::array<MultiFab, AMREX_SPACEDIM> >& sflux,
    Vector<std::array<MultiFab, AMREX_SPACEDIM> >& sedge,
    const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& umac,
    const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& w0mac,
    const BaseState<Real>& rho0_old_in, const BaseState<Real>& rho0_edge_old,
    const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& r0mac_old,
    const BaseState<Real>& rho0_new_in, const BaseState<Real>& rho0_edge_new,
    const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& r0mac_new,
    const BaseState<Real>& rhoh0_old_in, const BaseState<Real>& rhoh0_edge_old,
    const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& rh0mac_old,
    const BaseState<Real>& rhoh0_new_in, const BaseState<Real>& rhoh0_edge_new,
    const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& rh0mac_new,
    const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& h0mac_old,
    const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& h0mac_new) {

    amrex::ignore_unused(w0mac);
    amrex::ignore_unused(r0mac_old);
    amrex::ignore_unused(r0mac_new);
    amrex::ignore_unused(rh0mac_old);
    amrex::ignore_unused(rh0mac_new);
    amrex::ignore_unused(h0mac_old);
    amrex::ignore_unused(h0mac_new);

    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeRhoHFlux()", MakeRhoHFlux);

    const bool have_h = enthalpy_pred_type == pred_h ||
                        enthalpy_pred_type == pred_T_then_h ||
                        enthalpy_pred_type == pred_Tprime_then_h;
#if (AMREX_SPACEDIM == 3)
    const bool have_hprime = enthalpy_pred_type == pred_hprime;
#endif
    const bool have_rhoh = enthalpy_pred_type == pred_rhoh;

    const int species_pred_type_loc = species_pred_type;
    const int enthalpy_pred_type_loc = enthalpy_pred_type;

    for (int lev = 0; lev <= finest_level; ++lev) {
#if (AMREX_SPACEDIM == 3)
        MultiFab rho0mac_edgex, rho0mac_edgey, rho0mac_edgez;
        MultiFab h0mac_edgex, h0mac_edgey, h0mac_edgez;
        MultiFab rhoh0mac_edgex, rhoh0mac_edgey, rhoh0mac_edgez;

        rho0mac_edgex.define(convert(grids[lev], nodal_flag_x), dmap[lev], 1,
                             0);
        rho0mac_edgey.define(convert(grids[lev], nodal_flag_y), dmap[lev], 1,
                             0);
        rho0mac_edgez.define(convert(grids[lev], nodal_flag_z), dmap[lev], 1,
                             0);
        h0mac_edgex.define(convert(grids[lev], nodal_flag_x), dmap[lev], 1, 0);
        h0mac_edgey.define(convert(grids[lev], nodal_flag_y), dmap[lev], 1, 0);
        h0mac_edgez.define(convert(grids[lev], nodal_flag_z), dmap[lev], 1, 0);
        rhoh0mac_edgex.define(convert(grids[lev], nodal_flag_x), dmap[lev], 1,
                              0);
        rhoh0mac_edgey.define(convert(grids[lev], nodal_flag_y), dmap[lev], 1,
                              0);
        rhoh0mac_edgez.define(convert(grids[lev], nodal_flag_z), dmap[lev], 1,
                              0);

        rho0mac_edgex.setVal(0.);
        rho0mac_edgey.setVal(0.);
        rho0mac_edgez.setVal(0.);

        h0mac_edgex.setVal(0.);
        h0mac_edgey.setVal(0.);
        h0mac_edgez.setVal(0.);

        rhoh0mac_edgex.setVal(0.);
        rhoh0mac_edgey.setVal(0.);
        rhoh0mac_edgez.setVal(0.);

        if (spherical) {
            if (use_exact_base_state) {
                MultiFab::LinComb(rhoh0mac_edgex, 0.5, rh0mac_old[lev][0], 0,
                                  0.5, rh0mac_new[lev][0], 0, 0, 1, 0);
                MultiFab::LinComb(rhoh0mac_edgey, 0.5, rh0mac_old[lev][1], 0,
                                  0.5, rh0mac_new[lev][1], 0, 0, 1, 0);
                MultiFab::LinComb(rhoh0mac_edgez, 0.5, rh0mac_old[lev][2], 0,
                                  0.5, rh0mac_new[lev][2], 0, 0, 1, 0);
            } else {
                MultiFab::LinComb(rho0mac_edgex, 0.5, r0mac_old[lev][0], 0, 0.5,
                                  r0mac_new[lev][0], 0, 0, 1, 0);
                MultiFab::LinComb(rho0mac_edgey, 0.5, r0mac_old[lev][1], 0, 0.5,
                                  r0mac_new[lev][1], 0, 0, 1, 0);
                MultiFab::LinComb(rho0mac_edgez, 0.5, r0mac_old[lev][2], 0, 0.5,
                                  r0mac_new[lev][2], 0, 0, 1, 0);
                MultiFab::LinComb(h0mac_edgex, 0.5, h0mac_old[lev][0], 0, 0.5,
                                  h0mac_new[lev][0], 0, 0, 1, 0);
                MultiFab::LinComb(h0mac_edgey, 0.5, h0mac_old[lev][1], 0, 0.5,
                                  h0mac_new[lev][1], 0, 0, 1, 0);
                MultiFab::LinComb(h0mac_edgez, 0.5, h0mac_old[lev][2], 0, 0.5,
                                  h0mac_new[lev][2], 0, 0, 1, 0);
            }
        }
#endif

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(state[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region
            const Box& xbx = mfi.nodaltilebox(0);
            const Box& ybx = mfi.nodaltilebox(1);
#if (AMREX_SPACEDIM == 3)
            const Box& zbx = mfi.nodaltilebox(2);
#endif
            const Array4<Real> sedgex = sedge[lev][0].array(mfi);
            const Array4<Real> sfluxx = sflux[lev][0].array(mfi);
            const Array4<const Real> umacx = umac[lev][0].array(mfi);
            const Array4<Real> sedgey = sedge[lev][1].array(mfi);
            const Array4<Real> sfluxy = sflux[lev][1].array(mfi);
            const Array4<const Real> vmac = umac[lev][1].array(mfi);
#if (AMREX_SPACEDIM == 3)
            const Array4<Real> sedgez = sedge[lev][2].array(mfi);
            const Array4<Real> sfluxz = sflux[lev][2].array(mfi);
            const Array4<const Real> wmac = umac[lev][2].array(mfi);
#endif

            const auto rho0_old_arr = rho0_old_in.const_array();
            const auto rho0_new_arr = rho0_new_in.const_array();

            const auto rho0_edge_old_arr = rho0_edge_old.const_array();
            const auto rho0_edge_new_arr = rho0_edge_new.const_array();

            const auto rhoh0_old_arr = rhoh0_old_in.const_array();
            const auto rhoh0_new_arr = rhoh0_new_in.const_array();
            const auto rhoh0_edge_old_arr = rhoh0_edge_old.const_array();
            const auto rhoh0_edge_new_arr = rhoh0_edge_new.const_array();

#if (AMREX_SPACEDIM == 2)
            ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                // create x-fluxes
                if (have_h) {
                    // enthalpy edge state is h
                    if (species_pred_type_loc == pred_rhoprime_and_X) {
                        // density edge state is rho'
                        Real rho0_edge =
                            0.5 * (rho0_old_arr(lev, j) + rho0_new_arr(lev, j));

                        sfluxx(i, j, k, RhoH) =
                            umacx(i, j, k) *
                            (rho0_edge + sedgex(i, j, k, Rho)) *
                            sedgex(i, j, k, RhoH);

                    } else if (species_pred_type_loc == pred_rho_and_X ||
                               species_pred_type_loc == pred_rhoX) {
                        // density edge state is rho
                        sfluxx(i, j, k, RhoH) = umacx(i, j, k) *
                                                sedgex(i, j, k, Rho) *
                                                sedgex(i, j, k, RhoH);
                    }
                } else if (have_rhoh) {
                    sfluxx(i, j, k, RhoH) =
                        umacx(i, j, k) * sedgex(i, j, k, RhoH);

                } else if (enthalpy_pred_type_loc == pred_rhohprime ||
                           enthalpy_pred_type_loc == pred_T_then_rhohprime) {
                    // enthalpy edge state is (rho h)'
                    Real rhoh0_edge =
                        0.5 * (rhoh0_old_arr(lev, j) + rhoh0_new_arr(lev, j));

                    sfluxx(i, j, k, RhoH) =
                        umacx(i, j, k) * (rhoh0_edge + sedgex(i, j, k, RhoH));
                }
            });

            ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                // create y-fluxes
                if (have_h) {
                    // enthalpy edge state is h
                    if (species_pred_type_loc == pred_rhoprime_and_X) {
                        // density edge state is rho'
                        Real rho0_edge = 0.5 * (rho0_edge_old_arr(lev, j) +
                                                rho0_edge_new_arr(lev, j));

                        sfluxy(i, j, k, RhoH) =
                            vmac(i, j, k) * (rho0_edge + sedgey(i, j, k, Rho)) *
                            sedgey(i, j, k, RhoH);

                    } else if (species_pred_type_loc == pred_rho_and_X ||
                               species_pred_type_loc == pred_rhoX) {
                        // density edge state is rho
                        sfluxy(i, j, k, RhoH) = vmac(i, j, k) *
                                                sedgey(i, j, k, Rho) *
                                                sedgey(i, j, k, RhoH);
                    }
                } else if (have_rhoh) {
                    sfluxy(i, j, k, RhoH) =
                        vmac(i, j, k) * sedgey(i, j, k, RhoH);

                } else if (enthalpy_pred_type_loc == pred_rhohprime ||
                           enthalpy_pred_type_loc == pred_T_then_rhohprime) {
                    // enthalpy edge state is (rho h)'
                    Real rhoh0_edge = 0.5 * (rhoh0_edge_old_arr(lev, j) +
                                             rhoh0_edge_new_arr(lev, j));

                    sfluxy(i, j, k, RhoH) =
                        vmac(i, j, k) * (sedgey(i, j, k, RhoH) + rhoh0_edge);
                }
            });

#elif (AMREX_SPACEDIM == 3)

            if (!spherical) {
                ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    // create x-fluxes
                    if (have_h) {
                        // enthalpy edge state is h
                        if (species_pred_type_loc == pred_rhoprime_and_X) {
                            // density edge state is rho'
                            Real rho0_edge = 0.5 * (rho0_old_arr(lev, k) +
                                                    rho0_new_arr(lev, k));

                            sfluxx(i, j, k, RhoH) =
                                umacx(i, j, k) *
                                (rho0_edge + sedgex(i, j, k, Rho)) *
                                sedgex(i, j, k, RhoH);

                        } else if (species_pred_type_loc == pred_rho_and_X ||
                                   species_pred_type_loc == pred_rhoX) {
                            // density edge state is rho
                            sfluxx(i, j, k, RhoH) = umacx(i, j, k) *
                                                    sedgex(i, j, k, Rho) *
                                                    sedgex(i, j, k, RhoH);
                        }
                    } else if (have_rhoh) {
                        sfluxx(i, j, k, RhoH) =
                            umacx(i, j, k) * sedgex(i, j, k, RhoH);

                    } else if (enthalpy_pred_type_loc == pred_rhohprime ||
                               enthalpy_pred_type_loc ==
                                   pred_T_then_rhohprime) {
                        // enthalpy edge state is (rho h)'
                        Real rhoh0_edge = 0.5 * (rhoh0_old_arr(lev, k) +
                                                 rhoh0_new_arr(lev, k));

                        sfluxx(i, j, k, RhoH) =
                            umacx(i, j, k) *
                            (rhoh0_edge + sedgex(i, j, k, RhoH));
                    }
                });

                ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    // create y-fluxes
                    if (have_h) {
                        // enthalpy edge state is h
                        if (species_pred_type_loc == pred_rhoprime_and_X) {
                            // density edge state is rho'
                            Real rho0_edge = 0.5 * (rho0_old_arr(lev, k) +
                                                    rho0_new_arr(lev, k));

                            sfluxy(i, j, k, RhoH) =
                                vmac(i, j, k) *
                                (rho0_edge + sedgey(i, j, k, Rho)) *
                                sedgey(i, j, k, RhoH);

                        } else if (species_pred_type_loc == pred_rho_and_X ||
                                   species_pred_type_loc == pred_rhoX) {
                            // density edge state is rho
                            sfluxy(i, j, k, RhoH) = vmac(i, j, k) *
                                                    sedgey(i, j, k, Rho) *
                                                    sedgey(i, j, k, RhoH);
                        }
                    } else if (have_rhoh) {
                        sfluxy(i, j, k, RhoH) =
                            vmac(i, j, k) * sedgey(i, j, k, RhoH);

                    } else if (enthalpy_pred_type_loc == pred_rhohprime ||
                               enthalpy_pred_type_loc ==
                                   pred_T_then_rhohprime) {
                        // enthalpy edge state is (rho h)'
                        Real rhoh0_edge = 0.5 * (rhoh0_old_arr(lev, k) +
                                                 rhoh0_new_arr(lev, k));

                        sfluxy(i, j, k, RhoH) =
                            vmac(i, j, k) *
                            (rhoh0_edge + sedgey(i, j, k, RhoH));
                    }
                });

                ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    // create z-fluxes
                    if (have_h) {
                        // enthalpy edge state is h
                        if (species_pred_type_loc == pred_rhoprime_and_X) {
                            // density edge state is rho'
                            Real rho0_edge = 0.5 * (rho0_edge_old_arr(lev, k) +
                                                    rho0_edge_new_arr(lev, k));

                            sfluxz(i, j, k, RhoH) =
                                wmac(i, j, k) *
                                (rho0_edge + sedgez(i, j, k, Rho)) *
                                sedgez(i, j, k, RhoH);

                        } else if (species_pred_type_loc == pred_rho_and_X ||
                                   species_pred_type_loc == pred_rhoX) {
                            // density edge state is rho
                            sfluxz(i, j, k, RhoH) = wmac(i, j, k) *
                                                    sedgez(i, j, k, Rho) *
                                                    sedgez(i, j, k, RhoH);
                        }
                    } else if (have_rhoh) {
                        sfluxz(i, j, k, RhoH) =
                            wmac(i, j, k) * sedgez(i, j, k, RhoH);

                    } else if (enthalpy_pred_type_loc == pred_rhohprime ||
                               enthalpy_pred_type_loc ==
                                   pred_T_then_rhohprime) {
                        // enthalpy edge state is (rho h)'
                        Real rhoh0_edge = 0.5 * (rhoh0_edge_old_arr(lev, k) +
                                                 rhoh0_edge_new_arr(lev, k));

                        sfluxz(i, j, k, RhoH) =
                            wmac(i, j, k) *
                            (sedgez(i, j, k, RhoH) + rhoh0_edge);
                    }
                });
            } else {
                if (use_exact_base_state) {
                    const Array4<const Real> rhoh0_edgex =
                        rhoh0mac_edgex.array(mfi);
                    const Array4<const Real> rhoh0_edgey =
                        rhoh0mac_edgey.array(mfi);
                    const Array4<const Real> rhoh0_edgez =
                        rhoh0mac_edgez.array(mfi);

                    ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        if (have_h) {
                            // enthalpy edge state is h
                            // this is not supported on irregular-spaced base state
                        } else if (have_hprime) {
                            // enthalpy edge state is h'
                            // this is not supported on irregular-spaced base state
                        } else if (have_rhoh) {
                            sfluxx(i, j, k, RhoH) =
                                umacx(i, j, k) * sedgex(i, j, k, RhoH);
                        } else {
                            // enthalpy edge state is (rho h)'

                            // Average (rho h) onto edges by averaging rho and h
                            // separately onto edges.
                            //  (rho h)_edge = (rho h)' + rhoh_0

                            sfluxx(i, j, k, RhoH) =
                                umacx(i, j, k) *
                                (rhoh0_edgex(i, j, k) + sedgex(i, j, k, RhoH));
                        }
                    });

                    ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        if (have_h) {
                            // enthalpy edge state is h
                            // this is not supported on irregular-spaced base state
                        } else if (have_hprime) {
                            // enthalpy edge state is h'
                            // this is not supported on irregular-spaced base state
                        } else if (have_rhoh) {
                            sfluxy(i, j, k, RhoH) =
                                vmac(i, j, k) * sedgey(i, j, k, RhoH);
                        } else {
                            // enthalpy edge state is (rho h)'

                            // Average (rho h) onto edges by averaging rho and h
                            // separately onto edges.
                            //  (rho h)_edge = (rho h)' + rhoh_0
                            sfluxy(i, j, k, RhoH) =
                                vmac(i, j, k) *
                                (rhoh0_edgey(i, j, k) + sedgey(i, j, k, RhoH));
                        }
                    });

                    ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        if (have_h) {
                            // enthalpy edge state is h
                            // this is not supported on irregular-spaced base state
                        } else if (have_hprime) {
                            // enthalpy edge state is h'
                            // this is not supported on irregular-spaced base state
                        } else if (have_rhoh) {
                            sfluxz(i, j, k, RhoH) =
                                wmac(i, j, k) * sedgez(i, j, k, RhoH);
                        } else {
                            // enthalpy edge state is (rho h)'

                            // Average (rho h) onto edges by averaging rho and h
                            // separately onto edges.
                            //  (rho h)_edge = (rho h)' + rhoh_0
                            sfluxz(i, j, k, RhoH) =
                                wmac(i, j, k) *
                                (rhoh0_edgez(i, j, k) + sedgez(i, j, k, RhoH));
                        }
                    });
                } else {
                    const Array4<const Real> rho0_edgex =
                        rho0mac_edgex.array(mfi);
                    const Array4<const Real> rho0_edgey =
                        rho0mac_edgey.array(mfi);
                    const Array4<const Real> rho0_edgez =
                        rho0mac_edgez.array(mfi);

                    const Array4<const Real> h0_edgex = h0mac_edgex.array(mfi);
                    const Array4<const Real> h0_edgey = h0mac_edgey.array(mfi);
                    const Array4<const Real> h0_edgez = h0mac_edgez.array(mfi);

                    ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        if (have_h) {
                            // enthalpy edge state is h
                            if (species_pred_type_loc == pred_rhoprime_and_X) {
                                // density edge state is rho'
                                sfluxx(i, j, k, RhoH) = umacx(i, j, k) *
                                                        (rho0_edgex(i, j, k) +
                                                         sedgex(i, j, k, Rho)) *
                                                        sedgex(i, j, k, RhoH);

                            } else if (species_pred_type_loc ==
                                           pred_rho_and_X ||
                                       species_pred_type_loc == pred_rhoX) {
                                // density edge state is rho
                                sfluxx(i, j, k, RhoH) = umacx(i, j, k) *
                                                        sedgex(i, j, k, Rho) *
                                                        sedgex(i, j, k, RhoH);
                            }
                        } else if (have_hprime) {
                            // enthalpy edge state is h'
                            if (species_pred_type_loc == pred_rhoprime_and_X) {
                                // density edge state is rho'

                                // (rho h)_edge = (h' + h_0) * (rho' + rho_0) where h0 is
                                // computed from (rho h)_0 / rho_0
                                // sfluxx = (umac(i,j,k)+w0macx(i,j,k)) * (rho h)_edge
                                sfluxx(i, j, k, RhoH) =
                                    umacx(i, j, k) *
                                    (sedgex(i, j, k, Rho) +
                                     rho0_edgex(i, j, k)) *
                                    (sedgex(i, j, k, RhoH) + h0_edgex(i, j, k));
                            }

                        } else if (have_rhoh) {
                            sfluxx(i, j, k, RhoH) =
                                umacx(i, j, k) * sedgex(i, j, k, RhoH);

                        } else {
                            // enthalpy edge state is (rho h)'

                            // Average (rho h) onto edges by averaging rho and h
                            // separately onto edges.
                            //  (rho h)_edge = (rho h)' + (rho_0 * h_0)
                            // where h_0 is computed from (rho h)_0 / rho_0

                            sfluxx(i, j, k, RhoH) =
                                umacx(i, j, k) *
                                (rho0_edgex(i, j, k) * h0_edgex(i, j, k) +
                                 sedgex(i, j, k, RhoH));
                        }
                    });

                    ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        if (have_h) {
                            // enthalpy edge state is h
                            if (species_pred_type_loc == pred_rhoprime_and_X) {
                                // density edge state is rho'
                                sfluxy(i, j, k, RhoH) = vmac(i, j, k) *
                                                        (rho0_edgey(i, j, k) +
                                                         sedgey(i, j, k, Rho)) *
                                                        sedgey(i, j, k, RhoH);

                            } else if (species_pred_type_loc ==
                                           pred_rho_and_X ||
                                       species_pred_type_loc == pred_rhoX) {
                                // density edge state is rho
                                sfluxy(i, j, k, RhoH) = vmac(i, j, k) *
                                                        sedgey(i, j, k, Rho) *
                                                        sedgey(i, j, k, RhoH);
                            }
                        } else if (have_hprime) {
                            // enthalpy edge state is h'
                            if (species_pred_type_loc == pred_rhoprime_and_X) {
                                // density edge state is rho'

                                // (rho h)_edge = (h' + h_0) * (rho' + rho_0) where h0 is
                                // computed from (rho h)_0 / rho_0
                                // sfluxx = (umac(i,j,k)+w0macx(i,j,k)) * (rho h)_edge
                                sfluxy(i, j, k, RhoH) =
                                    vmac(i, j, k) *
                                    (sedgey(i, j, k, Rho) +
                                     rho0_edgey(i, j, k)) *
                                    (sedgey(i, j, k, RhoH) + h0_edgey(i, j, k));
                            }
                        } else if (have_rhoh) {
                            sfluxy(i, j, k, RhoH) =
                                vmac(i, j, k) * sedgey(i, j, k, RhoH);

                        } else {
                            // enthalpy edge state is (rho h)'

                            // Average (rho h) onto edges by averaging rho and h
                            // separately onto edges.
                            //  (rho h)_edge = (rho h)' + (rho_0 * h_0)
                            // where h_0 is computed from (rho h)_0 / rho_0
                            sfluxy(i, j, k, RhoH) =
                                vmac(i, j, k) *
                                (rho0_edgey(i, j, k) * h0_edgey(i, j, k) +
                                 sedgey(i, j, k, RhoH));
                        }
                    });

                    ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        if (have_h) {
                            // enthalpy edge state is h
                            if (species_pred_type_loc == pred_rhoprime_and_X) {
                                // density edge state is rho'
                                sfluxz(i, j, k, RhoH) = wmac(i, j, k) *
                                                        (rho0_edgez(i, j, k) +
                                                         sedgez(i, j, k, Rho)) *
                                                        sedgez(i, j, k, RhoH);

                            } else if (species_pred_type_loc ==
                                           pred_rho_and_X ||
                                       species_pred_type_loc == pred_rhoX) {
                                // density edge state is rho
                                sfluxz(i, j, k, RhoH) = wmac(i, j, k) *
                                                        sedgez(i, j, k, Rho) *
                                                        sedgez(i, j, k, RhoH);
                            }
                        } else if (have_hprime) {
                            // enthalpy edge state is h'
                            if (species_pred_type_loc == pred_rhoprime_and_X) {
                                // density edge state is rho'

                                // (rho h)_edge = (h' + h_0) * (rho' + rho_0) where h0 is
                                // computed from (rho h)_0 / rho_0
                                // sfluxx = (umac(i,j,k)+w0macx(i,j,k)) * (rho h)_edge
                                sfluxz(i, j, k, RhoH) =
                                    wmac(i, j, k) *
                                    (sedgez(i, j, k, Rho) +
                                     rho0_edgez(i, j, k)) *
                                    (sedgez(i, j, k, RhoH) + h0_edgez(i, j, k));
                            }
                        } else if (have_rhoh) {
                            sfluxz(i, j, k, RhoH) =
                                wmac(i, j, k) * sedgez(i, j, k, RhoH);

                        } else {
                            // enthalpy edge state is (rho h)'

                            // Average (rho h) onto edges by averaging rho and h
                            // separately onto edges.
                            //  (rho h)_edge = (rho h)' + (rho_0 * h_0)
                            // where h_0 is computed from (rho h)_0 / rho_0
                            sfluxz(i, j, k, RhoH) =
                                wmac(i, j, k) *
                                (rho0_edgez(i, j, k) * h0_edgez(i, j, k) +
                                 sedgez(i, j, k, RhoH));
                        }
                    });
                }
            }
#endif
        }  // end MFIter loop

        // increment or decrement the flux registers by area and time-weighted fluxes
        // Note that the fluxes need to be scaled by dt and area
        // In this example we are solving s_t = -div(+F)
        // The fluxes contain, e.g., F_{i+1/2,j} = (s*u)_{i+1/2,j}
        // Keep this in mind when considering the different sign convention for updating
        // the flux registers from the coarse or fine grid perspective
        // NOTE: the flux register associated with flux_reg_s[lev] is associated
        // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
        if (reflux_type == 2) {
            // Get the grid size
            const Real* dx = geom[lev].CellSize();
            // NOTE: areas are different in DIM=2 and DIM=3
#if (AMREX_SPACEDIM == 3)
            const Real area[3] = {dx[1] * dx[2], dx[0] * dx[2], dx[0] * dx[1]};
#else
            const Real area[2] = {dx[1], dx[0]};
#endif

            if (flux_reg_s[lev + 1]) {
                for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                    // update the lev+1/lev flux register (index lev+1)
                    flux_reg_s[lev + 1]->CrseInit(sflux[lev][i], i, RhoH, RhoH,
                                                  1, -1.0 * dt * area[i]);
                }
            }
            if (flux_reg_s[lev]) {
                for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                    // update the lev/lev-1 flux register (index lev)
                    flux_reg_s[lev]->FineAdd(sflux[lev][i], i, RhoH, RhoH, 1,
                                             1.0 * dt * area[i]);
                }
            }
        }
    }  // end loop over levels

    if (reflux_type == 1) {
        AverageDownFaces(sflux);
    }

    // Something analogous to edge_restriction is done in UpdateScal()
}

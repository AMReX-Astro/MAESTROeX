
#include <Maestro.H>
#include <MaestroBCThreads.H>
#include <Maestro_F.H>

using namespace amrex;

// compute unprojected mac velocities
void Maestro::AdvancePremac(
    Vector<std::array<MultiFab, AMREX_SPACEDIM> >& umac,
    const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& w0mac,
    const Vector<MultiFab>& w0_force_cart) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::AdvancePremac()", AdvancePremac);

    // create a uold with filled ghost cells
    Vector<MultiFab> utilde(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        utilde[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, ng_adv);
        utilde[lev].setVal(0.);
    }

    FillPatch(t_new, utilde, uold, uold, 0, 0, AMREX_SPACEDIM, 0, bcs_u, 1);

    // create a MultiFab to hold uold + w0
    Vector<MultiFab> ufull(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        ufull[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, ng_adv);
        // needed to avoid NaNs in filling corner ghost cells with 2 physical boundaries
        ufull[lev].setVal(0.);
    }

    // create ufull = uold + w0
    for (int lev = 0; lev <= finest_level; ++lev) {
        MultiFab::Copy(ufull[lev], w0_cart[lev], 0, 0, AMREX_SPACEDIM, 0);
    }
    // fill ufull ghost cells
    FillPatch(t_old, ufull, ufull, ufull, 0, 0, AMREX_SPACEDIM, 0, bcs_u, 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        MultiFab::Add(ufull[lev], utilde[lev], 0, 0, AMREX_SPACEDIM, ng_adv);
    }

    // create a face-centered MultiFab to hold utrans
    Vector<std::array<MultiFab, AMREX_SPACEDIM> > utrans(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        utrans[lev][0].define(convert(grids[lev], nodal_flag_x), dmap[lev], 1,
                              1);
        utrans[lev][1].define(convert(grids[lev], nodal_flag_y), dmap[lev], 1,
                              1);
#if (AMREX_SPACEDIM == 3)
        utrans[lev][2].define(convert(grids[lev], nodal_flag_z), dmap[lev], 1,
                              1);
#endif
        for (int j = 0; j < AMREX_SPACEDIM; j++) {
            utrans[lev][j].setVal(0.);
        }
    }

    // create utrans
    MakeUtrans(utilde, ufull, utrans, w0mac);

    // create a MultiFab to hold the velocity forcing
    Vector<MultiFab> vel_force(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        if (ppm_trace_forces == 0) {
            vel_force[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        } else {
            // tracing needs more ghost cells
            vel_force[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, ng_s);
        }
        vel_force[lev].setVal(0.);
    }

    int do_add_utilde_force = 1;
    MakeVelForce(vel_force, utrans, sold, rho0_old, grav_cell_old,
                 w0_force_cart,
#ifdef ROTATION
                 w0mac, false,
#endif
                 do_add_utilde_force);

    // add w0 to trans velocities
    Addw0(utrans, w0mac, 1.);

    VelPred(utilde, ufull, utrans, umac, w0mac, vel_force);
}

void Maestro::UpdateScal(
    const Vector<MultiFab>& stateold, Vector<MultiFab>& statenew,
    const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& sflux,
    const Vector<MultiFab>& force, int start_comp, int num_comp,
    const Vector<MultiFab>& p0_cart) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::UpdateScal()", UpdateScal);

    const Real dt_loc = dt;

    for (int lev = 0; lev <= finest_level; ++lev) {
        const auto dx = geom[lev].CellSizeArray();

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(stateold[lev], TilingIfNotGPU()); mfi.isValid();
             ++mfi) {
            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const Array4<const Real> sold_arr = stateold[lev].array(mfi);
            const Array4<Real> snew_arr = statenew[lev].array(mfi);
            const Array4<const Real> sfluxx = sflux[lev][0].array(mfi);
            const Array4<const Real> sfluxy = sflux[lev][1].array(mfi);
#if (AMREX_SPACEDIM == 3)
            const Array4<const Real> sfluxz = sflux[lev][2].array(mfi);
#endif
            const Array4<const Real> force_arr = force[lev].array(mfi);

            if (start_comp == RhoH) {
                // Enthalpy update
                const Array4<const Real> p0_arr = p0_cart[lev].array(mfi);

                ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    Real divterm =
                        ((sfluxx(i + 1, j, k, RhoH) - sfluxx(i, j, k, RhoH)) /
                             dx[0] +
                         (sfluxy(i, j + 1, k, RhoH) - sfluxy(i, j, k, RhoH)) /
                             dx[1]
#if (AMREX_SPACEDIM == 3)
                         + (sfluxz(i, j, k + 1, RhoH) - sfluxz(i, j, k, RhoH)) /
                               dx[2]
#endif
                        );

                    snew_arr(i, j, k, RhoH) =
                        sold_arr(i, j, k, RhoH) +
                        dt_loc * (-divterm + force_arr(i, j, k, RhoH));

                    if (do_eos_h_above_cutoff &&
                        snew_arr(i, j, k, Rho) <= base_cutoff_density) {
                        eos_t eos_state;

                        eos_state.rho = snew_arr(i, j, k, Rho);
                        eos_state.T = sold_arr(i, j, k, Temp);
                        eos_state.p = p0_arr(i, j, k);
                        for (auto n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = snew_arr(i, j, k, FirstSpec + n) /
                                              eos_state.rho;
                        }
#if NAUX_NET > 0
                        for (auto n = 0; n < NumAux; ++n) {
                            eos_state.aux[n] =
                                snew_arr(i, j, k, FirstAux + n) / eos_state.rho;
                        }
#endif

                        eos(eos_input_rp, eos_state);

                        snew_arr(i, j, k, RhoH) =
                            snew_arr(i, j, k, Rho) * eos_state.h;
                    }
                });

            } else if (start_comp == FirstSpec) {
                // RhoX update

                ParallelFor(tileBox, NumSpec,
                            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                                int comp = FirstSpec + n;

                                Real divterm = (sfluxx(i + 1, j, k, comp) -
                                                sfluxx(i, j, k, comp)) /
                                               dx[0];
                                divterm += (sfluxy(i, j + 1, k, comp) -
                                            sfluxy(i, j, k, comp)) /
                                           dx[1];
#if (AMREX_SPACEDIM == 3)
                                divterm += (sfluxz(i, j, k + 1, comp) -
                                            sfluxz(i, j, k, comp)) /
                                           dx[2];
#endif
                                snew_arr(i, j, k, comp) =
                                    sold_arr(i, j, k, comp) +
                                    dt_loc *
                                        (-divterm + force_arr(i, j, k, comp));
                            });

                ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    // update density
                    snew_arr(i, j, k, Rho) = sold_arr(i, j, k, Rho);

                    bool has_negative_species = false;

                    // define the update to rho as the sum of the updates to (rho X)_i
                    for (int comp = start_comp; comp < start_comp + NumSpec;
                         ++comp) {
                        snew_arr(i, j, k, Rho) +=
                            snew_arr(i, j, k, comp) - sold_arr(i, j, k, comp);
                        if (snew_arr(i, j, k, comp) < 0.0) {
                            has_negative_species = true;
                        }
                    }

// update auxiliary variables
#if NAUX_NET > 0
                    for (int comp = FirstAux; comp < FirstAux + NumAux;
                         ++comp) {
                        snew_arr(i, j, k, comp) = sold_arr(i, j, k, comp) *
                                                  snew_arr(i, j, k, Rho) /
                                                  sold_arr(i, j, k, Rho);
                    }
#endif

                    // enforce a density floor
                    if (snew_arr(i, j, k, Rho) < 0.5 * base_cutoff_density) {
                        for (int comp = start_comp; comp < start_comp + NumSpec;
                             ++comp) {
                            snew_arr(i, j, k, comp) *= 0.5 *
                                                       base_cutoff_density /
                                                       snew_arr(i, j, k, Rho);
                        }
                        snew_arr(i, j, k, Rho) = 0.5 * base_cutoff_density;
                    }

                    // do not allow the species to leave here negative.
                    if (has_negative_species) {
                        for (int comp = start_comp; comp < start_comp + NumSpec;
                             ++comp) {
                            if (snew_arr(i, j, k, comp) < 0.0) {
                                Real delta = -snew_arr(i, j, k, comp);
                                Real sumX = 0.0;
                                for (int comp2 = start_comp;
                                     comp2 < start_comp + NumSpec; ++comp2) {
                                    if (comp2 != comp &&
                                        snew_arr(i, j, k, comp2) >= 0.0) {
                                        sumX += snew_arr(i, j, k, comp2);
                                    }
                                }
                                for (int comp2 = start_comp;
                                     comp2 < start_comp + NumSpec; ++comp2) {
                                    if (comp2 != comp &&
                                        snew_arr(i, j, k, comp2) >= 0.0) {
                                        Real frac =
                                            snew_arr(i, j, k, comp2) / sumX;
                                        snew_arr(i, j, k, comp2) -=
                                            frac * delta;
                                    }
                                }
                                snew_arr(i, j, k, comp) = 0.0;
                            }
                        }
                    }
                });
            } else {
                Abort("Invalid scalar in UpdateScal().");
            }  // }
        }      // end MFIter loop
    }          // end loop over levels

    // synchronize by refluxing and averaging down, starting from the finest_level-1/finest_level pair
    if (reflux_type == 2) {
        for (int lev = finest_level - 1; lev >= 0; --lev) {
            // update lev based on coarse-fine flux mismatch
            flux_reg_s[lev + 1]->Reflux(statenew[lev], 1.0, start_comp,
                                        start_comp, num_comp, geom[lev]);
            if (start_comp == FirstSpec) {
                // do the same for density if we updated the species
                flux_reg_s[lev + 1]->Reflux(statenew[lev], 1.0, Rho, Rho, 1,
                                            geom[lev]);

// and the aux variables
#if NAUX_NET > 0
                for (int comp = 0; comp < NumAux; ++comp) {
                    flux_reg_s[lev + 1]->Reflux(statenew[lev], 1.0, FirstAux,
                                                FirstAux, NumAux, geom[lev]);
                }
#endif
            }
        }
    }

    // average fine data onto coarser cells
    // fill ghost cells
    AverageDown(statenew, start_comp, num_comp);
    FillPatch(t_old, statenew, statenew, statenew, start_comp, start_comp,
              num_comp, start_comp, bcs_s);

    // do the same for density and aux if we updated the species
    if (start_comp == FirstSpec) {
        AverageDown(statenew, Rho, 1);
        FillPatch(t_old, statenew, statenew, statenew, Rho, Rho, 1, Rho, bcs_s);

#if NAUX_NET > 0
        AverageDown(statenew, FirstAux, NumAux);
        FillPatch(t_old, statenew, statenew, statenew, FirstAux, FirstAux,
                  NumAux, FirstAux, bcs_s);
#endif
    }
}

void Maestro::UpdateVel(
    const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& umac,
    const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& uedge,
    const Vector<MultiFab>& force, const Vector<MultiFab>& sponge,
    [[maybe_unused]] const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& w0mac) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::UpdateVel()", UpdateVel);

    // 1) Subtract (Utilde dot grad) Utilde term from old Utilde
    // 2) Add forcing term to new Utilde

    const Real dt_loc = dt;

    for (int lev = 0; lev <= finest_level; ++lev) {
        const auto dx = geom[lev].CellSizeArray();

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(force[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const Array4<const Real> uold_arr = uold[lev].array(mfi);
            const Array4<Real> unew_arr = unew[lev].array(mfi);
            const Array4<const Real> umacx = umac[lev][0].array(mfi);
            const Array4<const Real> uedgex = uedge[lev][0].array(mfi);
            const Array4<const Real> vmac = umac[lev][1].array(mfi);
            const Array4<const Real> uedgey = uedge[lev][1].array(mfi);
#if (AMREX_SPACEDIM == 3)
            const Array4<const Real> wmac = umac[lev][2].array(mfi);
            const Array4<const Real> uedgez = uedge[lev][2].array(mfi);
#endif
            const Array4<const Real> force_arr = force[lev].array(mfi);
            const Array4<const Real> sponge_arr = sponge[lev].array(mfi);
            const Array4<const Real> w0_arr = w0_cart[lev].array(mfi);

            if (!spherical) {
                ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    // create cell-centered Utilde
                    Real ubar = 0.5 * (umacx(i, j, k) + umacx(i + 1, j, k));
                    Real vbar = 0.5 * (vmac(i, j, k) + vmac(i, j + 1, k));
#if (AMREX_SPACEDIM == 3)
                    Real wbar = 0.5 * (wmac(i, j, k) + wmac(i, j, k + 1));
#endif

                    // create (Utilde dot grad) Utilde
                    Real ugradu =
                        (ubar * (uedgex(i + 1, j, k, 0) - uedgex(i, j, k, 0)) /
                             dx[0] +
                         vbar * (uedgey(i, j + 1, k, 0) - uedgey(i, j, k, 0)) /
                             dx[1]
#if (AMREX_SPACEDIM == 3)
                         +
                         wbar * (uedgez(i, j, k + 1, 0) - uedgez(i, j, k, 0)) /
                             dx[2]
#endif
                        );

                    Real ugradv =
                        (ubar * (uedgex(i + 1, j, k, 1) - uedgex(i, j, k, 1)) /
                             dx[0] +
                         vbar * (uedgey(i, j + 1, k, 1) - uedgey(i, j, k, 1)) /
                             dx[1]
#if (AMREX_SPACEDIM == 3)
                         +
                         wbar * (uedgez(i, j, k + 1, 1) - uedgez(i, j, k, 1)) /
                             dx[2]
#endif
                        );

#if (AMREX_SPACEDIM == 3)
                    Real ugradw =
                        ubar * (uedgex(i + 1, j, k, 2) - uedgex(i, j, k, 2)) /
                            dx[0] +
                        vbar * (uedgey(i, j + 1, k, 2) - uedgey(i, j, k, 2)) /
                            dx[1] +
                        wbar * (uedgez(i, j, k + 1, 2) - uedgez(i, j, k, 2)) /
                            dx[2];
#endif

                    // update with (Utilde dot grad) Utilde and force
                    unew_arr(i, j, k, 0) = uold_arr(i, j, k, 0) -
                                           dt_loc * ugradu +
                                           dt_loc * force_arr(i, j, k, 0);
                    unew_arr(i, j, k, 1) = uold_arr(i, j, k, 1) -
                                           dt_loc * ugradv +
                                           dt_loc * force_arr(i, j, k, 1);
#if (AMREX_SPACEDIM == 3)
                    unew_arr(i, j, k, 2) = uold_arr(i, j, k, 2) -
                                           dt_loc * ugradw +
                                           dt_loc * force_arr(i, j, k, 2);
#endif

                    // subtract (w0 dot grad) Utilde term
#if (AMREX_SPACEDIM == 2)
                    Real w0bar =
                        0.5 * (w0_arr(i, j, k, AMREX_SPACEDIM - 1) +
                               w0_arr(i, j + 1, k, AMREX_SPACEDIM - 1));
#else 
                    Real w0bar = 0.5*(w0_arr(i,j,k,AMREX_SPACEDIM-1) + w0_arr(i,j,k+1,AMREX_SPACEDIM-1));
#endif

                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        unew_arr(i, j, k, n) -=
                            dt_loc * w0bar *
#if (AMREX_SPACEDIM == 2)
                            (uedgey(i, j + 1, k, n) - uedgey(i, j, k, n)) /
                            dx[1];
#else
                        (uedgez(i,j,k+1,n) - uedgez(i,j,k,n))/dx[2];
#endif
                        // Add the sponge
                        if (do_sponge) {
                            unew_arr(i, j, k, n) *= sponge_arr(i, j, k);
                        }
                    }
                });
            } else {
#if (AMREX_SPACEDIM == 3)
                const Array4<const Real> w0macx = w0mac[lev][0].array(mfi);
                const Array4<const Real> w0macy = w0mac[lev][1].array(mfi);
                const Array4<const Real> w0macz = w0mac[lev][2].array(mfi);

                ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    // create cell-centered Utilde
                    Real ubar = 0.5 * (umacx(i, j, k) + umacx(i + 1, j, k));
                    Real vbar = 0.5 * (vmac(i, j, k) + vmac(i, j + 1, k));
                    Real wbar = 0.5 * (wmac(i, j, k) + wmac(i, j, k + 1));

                    // create (Utilde dot grad) Utilde
                    Real ugradu =
                        ubar * (uedgex(i + 1, j, k, 0) - uedgex(i, j, k, 0)) /
                            dx[0] +
                        vbar * (uedgey(i, j + 1, k, 0) - uedgey(i, j, k, 0)) /
                            dx[1] +
                        wbar * (uedgez(i, j, k + 1, 0) - uedgez(i, j, k, 0)) /
                            dx[2];

                    Real ugradv =
                        ubar * (uedgex(i + 1, j, k, 1) - uedgex(i, j, k, 1)) /
                            dx[0] +
                        vbar * (uedgey(i, j + 1, k, 1) - uedgey(i, j, k, 1)) /
                            dx[1] +
                        wbar * (uedgez(i, j, k + 1, 1) - uedgez(i, j, k, 1)) /
                            dx[2];

                    Real ugradw =
                        ubar * (uedgex(i + 1, j, k, 2) - uedgex(i, j, k, 2)) /
                            dx[0] +
                        vbar * (uedgey(i, j + 1, k, 2) - uedgey(i, j, k, 2)) /
                            dx[1] +
                        wbar * (uedgez(i, j, k + 1, 2) - uedgez(i, j, k, 2)) /
                            dx[2];

                    // update with (Utilde dot grad) Utilde and force
                    unew_arr(i, j, k, 0) = uold_arr(i, j, k, 0) -
                                           dt_loc * ugradu +
                                           dt_loc * force_arr(i, j, k, 0);
                    unew_arr(i, j, k, 1) = uold_arr(i, j, k, 1) -
                                           dt_loc * ugradv +
                                           dt_loc * force_arr(i, j, k, 1);
                    unew_arr(i, j, k, 2) = uold_arr(i, j, k, 2) -
                                           dt_loc * ugradw +
                                           dt_loc * force_arr(i, j, k, 2);

                    // Subtract (w0 dot grad) Utilde term from new Utilde
                    Real gradux =
                        (uedgex(i + 1, j, k, 0) - uedgex(i, j, k, 0)) / dx[0];
                    Real gradvx =
                        (uedgex(i + 1, j, k, 1) - uedgex(i, j, k, 1)) / dx[0];
                    Real gradwx =
                        (uedgex(i + 1, j, k, 2) - uedgex(i, j, k, 2)) / dx[0];

                    Real graduy =
                        (uedgey(i, j + 1, k, 0) - uedgey(i, j, k, 0)) / dx[1];
                    Real gradvy =
                        (uedgey(i, j + 1, k, 1) - uedgey(i, j, k, 1)) / dx[1];
                    Real gradwy =
                        (uedgey(i, j + 1, k, 2) - uedgey(i, j, k, 2)) / dx[1];

                    Real graduz =
                        (uedgez(i, j, k + 1, 0) - uedgez(i, j, k, 0)) / dx[2];
                    Real gradvz =
                        (uedgez(i, j, k + 1, 1) - uedgez(i, j, k, 1)) / dx[2];
                    Real gradwz =
                        (uedgez(i, j, k + 1, 2) - uedgez(i, j, k, 2)) / dx[2];

                    Real w0_gradur =
                        gradux * 0.5 * (w0macx(i, j, k) + w0macx(i + 1, j, k)) +
                        graduy * 0.5 * (w0macy(i, j, k) + w0macy(i, j + 1, k)) +
                        graduz * 0.5 * (w0macz(i, j, k) + w0macz(i, j, k + 1));

                    Real w0_gradvr =
                        gradvx * 0.5 * (w0macx(i, j, k) + w0macx(i + 1, j, k)) +
                        gradvy * 0.5 * (w0macy(i, j, k) + w0macy(i, j + 1, k)) +
                        gradvz * 0.5 * (w0macz(i, j, k) + w0macz(i, j, k + 1));

                    Real w0_gradwr =
                        gradwx * 0.5 * (w0macx(i, j, k) + w0macx(i + 1, j, k)) +
                        gradwy * 0.5 * (w0macy(i, j, k) + w0macy(i, j + 1, k)) +
                        gradwz * 0.5 * (w0macz(i, j, k) + w0macz(i, j, k + 1));

                    unew_arr(i, j, k, 0) -= dt_loc * w0_gradur;
                    unew_arr(i, j, k, 1) -= dt_loc * w0_gradvr;
                    unew_arr(i, j, k, 2) -= dt_loc * w0_gradwr;

                    // Add the sponge
                    if (do_sponge) {
                        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                            unew_arr(i, j, k, n) *= sponge_arr(i, j, k);
                        }
                    }
                });
#else
                Abort("UpdateVel: Spherical is not valid for DIM < 3");
#endif
            }
        }  // end MFIter loop
    }      // end loop over levels

    // average fine data onto coarser cells
    AverageDown(unew, 0, AMREX_SPACEDIM);

    // fill ghost cells
    FillPatch(t_old, unew, unew, unew, 0, 0, AMREX_SPACEDIM, 0, bcs_u, 1);
}


#include <Maestro.H>
#include <Maestro_F.H>

#include <AMReX_MLMG.H>
#include <AMReX_MLNodeLaplacian.H>

using namespace amrex;

// Perform a nodal projection.
// Given a cell-centered velocity field Vproj (assembled in CreateUvecForProj),
// Vproj can decomposed into Vproj = Utilde + sig grad phi,
// where Utilde satisfies the constraint div(beta0*Utilde) = beta0*(S-Sbar)
// Depending on proj_type we use different values of Vproj, sig, and beta0
// to solve for phi we use div sigma grad phi = div(beta*Vproj) - beta0*(S-Sbar)
// where sigma = beta0 or beta0/rho depending on proj_type
// then solve for Utilde, pi, and grad(pi) based on proj_type.
// rhcc should enter as beta0*(S-Sbar) so we need to multiply by -1 in this routine
// the projection (done below)
void Maestro::NodalProj(int proj_type, Vector<MultiFab>& rhcc,
                        int istep_divu_iter, bool sdc_off, bool is_predictor) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::NodalProj()", NodalProj);

    AMREX_ASSERT(rhcc[0].nGrow() == 1);

    if (!(proj_type == initial_projection_comp ||
          proj_type == divu_iters_comp || proj_type == pressure_iters_comp ||
          proj_type == regular_timestep_comp)) {
        amrex::Abort("Maestro::NodalProj - invalid proj_type");
    }

    // build a multifab sig with 1 ghost cell
    Vector<MultiFab> sig(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        sig[lev].define(grids[lev], dmap[lev], 1, 1);
        // needed to avoid NaNs in filling corner ghost cells with 2 physical boundaries
        sig[lev].setVal(0.);
    }

    // fill in sig
    // initial_projection_comp: 1
    // divu_iters_comp:         1
    // pressure_iters_comp:     rho^1/2   -- (rhoold+rhonew)/2
    // regular_timestep_comp:   rho^n+1/2 -- (rhoold+rhonew)/2
    if (proj_type == initial_projection_comp || proj_type == divu_iters_comp) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            sig[lev].setVal(1.);
        }
    } else if (proj_type == pressure_iters_comp ||
               proj_type == regular_timestep_comp) {
        FillPatch(0.5 * (t_old + t_new), sig, sold, snew, Rho, 0, 1, Rho,
                  bcs_s);
    }

    // build a multifab Vproj (the quantity that will be projected) with 1 ghost cell
    Vector<MultiFab> Vproj(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        Vproj[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        // needed to avoid NaNs in filling corner ghost cells with 2 physical boundaries
        Vproj[lev].setVal(0.);
    }

    // fill in Vproj
    // initial_projection_comp: Utilde^0                        -- uold
    // divu_iters_comp:         Utilde^0                        -- uold
    // pressure_iters_comp:     (Utilde^n+1,* - Utilde^n)/dt    -- (unew-uold)/dt
    // regular_timestep_comp:   (Utilde^n+1,* + dt*gpi/rhohalf) -- unew + dt*gpi/rhohalf
    // note sig is only used for regular_timestep_comp, and currently holds rhohalf
    CreateUvecForProj(proj_type, Vproj, sig);

    bool using_alt_energy_fix = false;
    if (use_alt_energy_fix && proj_type != initial_projection_comp &&
        proj_type != divu_iters_comp) {
        using_alt_energy_fix = true;
    }

    // build a multifab to store a Cartesian version of beta0 with 1 ghost cell
    Vector<MultiFab> beta0_cart(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        beta0_cart[lev].define(grids[lev], dmap[lev], 1, 1);
    }

    // convert beta0 to multi-D MultiFab
    // initial_projection_comp: beta0_old
    // divu_iters_comp:         beta0_old
    // pressure_iters_comp:     (beta0_old+beta0_new)/2
    // regular_timestep_comp:   (beta0_old+beta0_new)/2
    if (proj_type == initial_projection_comp || proj_type == divu_iters_comp) {
        Put1dArrayOnCart(beta0_old, beta0_cart, false, false, bcs_f, 0);
    } else {
        auto beta0_nph = 0.5 * (beta0_old + beta0_new);
        Put1dArrayOnCart(beta0_nph, beta0_cart, false, false, bcs_f, 0);
    }

    // convert Vproj to beta0*Vproj
    for (int lev = 0; lev <= finest_level; ++lev) {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            MultiFab::Multiply(Vproj[lev], beta0_cart[lev], 0, dir, 1, 1);
        }
    }

    // invert sig then multiply by beta0 (or beta0^2 if alt_energy_fix)
    // so sig now holds:
    // initial_projection_comp: beta0
    // divu_iters_comp:         beta0
    // pressure_iters_comp:     beta0/rho
    // regular_timestep_comp:   beta0/rho

    for (int lev = 0; lev <= finest_level; ++lev) {
        sig[lev].invert(1.0, 1);
        MultiFab::Multiply(sig[lev], beta0_cart[lev], 0, 0, 1, 1);
        // if using_alt_energy_fix, multiply sig by beta0 again
        if (using_alt_energy_fix) {
            MultiFab::Multiply(sig[lev], beta0_cart[lev], 0, 0, 1, 1);
        }
    }

    /*
   set_outflow_bcs() was used in IAMR/Projection.cpp
   This routine modifies boundary conditions on phi at outflow.
   Typically these are homogeneous Dirichlet but in some cases you need
   to modify this.  First, if div(u) is "large enough" at outflow you need
   to modify the boundary conditions.  This often happens in combustion
   flame systems (PeleLM),  Second if you have hydrostatic effects
   due to gravity you need to modify the boundary conditions.  This happens
   in IAMR_type runs with strong gravity.  MAESTRO
   requires neither of these since we do not run problems with strong
   div(u) at outflow, and our pi (pressure) does not need to capture
   stratification since we have put this pressure in perturbational form
   by adding (rho-rho0)*g to the RHS of the velocity equation.
 */
    /*
    set_outflow_bcs()
 */

    /*
   SetBoundaryVelocity() is a simplified version of the function
   set_boundary_velocity() in IAMR/Projection.cpp
   basically this sets velocity in ghost cells to zero except at inflow
   right now, the nodal solver expects ghost cells to be filled this way
   in slightly more detail
   1) At non-inflow faces, the normal component of velocity will be completely zero'd
   2) If a face is an inflow face, then the normal velocity at corners just outside inflow faces
     will be zero'd outside of Neumann boundaries (slipWall, noSlipWall, Symmetry)
     BUT will retain non-zero values at periodic corners
 */
    SetBoundaryVelocity(Vproj);

    std::array<LinOpBCType, AMREX_SPACEDIM> mlmg_lobc;
    std::array<LinOpBCType, AMREX_SPACEDIM> mlmg_hibc;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (Geom(0).isPeriodic(idim)) {
            mlmg_lobc[idim] = mlmg_hibc[idim] = LinOpBCType::Periodic;
        } else {
            if (phys_bc[idim] == Outflow) {
                mlmg_lobc[idim] = LinOpBCType::Dirichlet;
            } else {
                mlmg_lobc[idim] = LinOpBCType::Neumann;
            }

            if (phys_bc[AMREX_SPACEDIM + idim] == Outflow) {
                mlmg_hibc[idim] = LinOpBCType::Dirichlet;
            } else {
                mlmg_hibc[idim] = LinOpBCType::Neumann;
            }
        }
    }

    LPInfo info;
    info.setMetricTerm(false);

    if (hg_bottom_solver == 4) {
        info.setAgglomeration(true);
        info.setConsolidation(true);
    } else {
        info.setAgglomeration(false);
        info.setConsolidation(false);
    }

    // Only pass up to defined level to prevent looping over undefined grids.
    MLNodeLaplacian mlndlap(Geom(0, finest_level), grids, dmap, info);
    mlndlap.setGaussSeidel(true);
    mlndlap.setHarmonicAverage(false);

    mlndlap.setDomainBC(mlmg_lobc, mlmg_hibc);

    // set sig in the MLNodeLaplacian object
    for (int ilev = 0; ilev <= finest_level; ++ilev) {
        mlndlap.setSigma(ilev, sig[ilev]);
    }

    // rhstotal is nodal and will contain the complete right-hand-side
    // rhtotal = div(Vproj) + rhcc (rhcc is averaged to nodes)
    Vector<MultiFab> rhstotal(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        rhstotal[lev].define(convert(grids[lev], nodal_flag), dmap[lev], 1, 0);
    }

    // phi is nodal and will contain the solution (need 1 ghost cell)
    Vector<MultiFab> phi(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        phi[lev].define(convert(grids[lev], nodal_flag), dmap[lev], 1, 1);
        phi[lev].setVal(0.);
    }

    // multiply rhcc = beta0*(S-Sbar) by -1 since we want
    // rhstotal to contain div(beta*Vproj) - beta0*(S-Sbar)
    for (int lev = 0; lev <= finest_level; ++lev) {
        rhcc[lev].mult(-1.0, 0, 1, 1);
    }

#ifdef AMREX_USE_CUDA
    // compRHS is non-deterministic on the GPU. If this flag is true, run
    // on the CPU instead
    bool launched;
    if (deterministic_nodal_solve) {
        launched = !Gpu::notInLaunchRegion();
        // turn off GPU
        if (launched) Gpu::setLaunchRegion(false);
    }
#endif
    // Assemble the nodal RHS as the sum of the cell-centered RHS averaged to nodes
    // plus div (beta0*Vproj) on nodes
    // so rhstotal = div(beta*Vproj) - beta0*(S-Sbar)
    mlndlap.compRHS(amrex::GetVecOfPtrs(rhstotal), amrex::GetVecOfPtrs(Vproj),
                    {},  // pass in null rhnd
                    amrex::GetVecOfPtrs(rhcc));
#ifdef AMREX_USE_CUDA
    if (deterministic_nodal_solve) {
        // turn GPU back on
        if (launched) Gpu::setLaunchRegion(true);
    }
#endif

    // restore rhcc
    for (int lev = 0; lev <= finest_level; ++lev) {
        rhcc[lev].mult(-1.0, 0, 1, 1);
    }

    MLMG mlmg(mlndlap);
    mlmg.setVerbose(mg_verbose);
    mlmg.setBottomVerbose(cg_verbose);

    Real abs_tol = -1.;  // disable absolute tolerance
    Real rel_tol = 1.e-3;

    // logic for choosing multigrid tolerance
    // parameters are defined in Maestro.cpp
    // we change the tolerance depending on what part of the algorithm we
    // are in, and other factors including planar vs. spherical, the number
    // of AMR levels, etc.
    if (proj_type == initial_projection_comp) {
        rel_tol = (spherical) ? eps_init_proj_sph : eps_init_proj_cart;
    } else if (proj_type == divu_iters_comp) {
        if (spherical) {
            if (istep_divu_iter == init_divu_iter) {
                rel_tol = eps_divu_sph;
            } else if (istep_divu_iter == init_divu_iter - 1) {
                rel_tol = eps_divu_sph * divu_iter_factor;
            } else if (istep_divu_iter <= init_divu_iter - 2) {
                rel_tol = eps_divu_sph * pow(divu_iter_factor, 2);
            }
        } else {
            if (istep_divu_iter == init_divu_iter) {
                rel_tol = amrex::min(
                    eps_divu_cart * pow(divu_level_factor, finest_level),
                    eps_divu_cart * pow(divu_level_factor, 2));
            } else if (istep_divu_iter == init_divu_iter - 1) {
                rel_tol = amrex::min(eps_divu_cart * divu_iter_factor *
                                         pow(divu_level_factor, finest_level),
                                     eps_divu_cart * divu_iter_factor *
                                         pow(divu_level_factor, 2));
            } else if (istep_divu_iter <= init_divu_iter - 2) {
                rel_tol = amrex::min(eps_divu_cart * pow(divu_iter_factor, 2) *
                                         pow(divu_level_factor, finest_level),
                                     eps_divu_cart * pow(divu_iter_factor, 2) *
                                         pow(divu_level_factor, 2));
            }
        }
    } else if (proj_type == pressure_iters_comp ||
               proj_type == regular_timestep_comp) {
        rel_tol =
            amrex::min(eps_hg_max, eps_hg * pow(hg_level_factor, finest_level));
    }

    // solve for phi
    Print() << "Calling nodal solver" << std::endl;
#ifdef AMREX_USE_CUDA
    if (deterministic_nodal_solve) {
        launched = !Gpu::notInLaunchRegion();
        // turn off GPU
        if (launched) Gpu::setLaunchRegion(false);
    }
#endif
    mlmg.solve(amrex::GetVecOfPtrs(phi), amrex::GetVecOfConstPtrs(rhstotal),
               rel_tol, abs_tol);
#ifdef AMREX_USE_CUDA
    if (deterministic_nodal_solve) {
        // turn GPU back on
        if (launched) Gpu::setLaunchRegion(true);
    }
#endif
    Print() << "Done calling nodal solver" << std::endl;

    // convert beta0*Vproj back to Vproj
    for (int lev = 0; lev <= finest_level; ++lev) {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            MultiFab::Divide(Vproj[lev], beta0_cart[lev], 0, dir, 1, 1);
        }
    }

    // divide sig by beta0 (or beta0^2 if alt_energy_fix)
    // so it now holds
    // initial_projection_comp: 1
    // divu_iters_comp:         1
    // pressure_iters_comp:     1/rho
    // regular_timestep_comp:   1/rho
    for (int lev = 0; lev <= finest_level; ++lev) {
        MultiFab::Divide(sig[lev], beta0_cart[lev], 0, 0, 1, 1);
        // if using_alt_energy_fix, divide sig by beta0 again
        if (using_alt_energy_fix) {
            MultiFab::Divide(sig[lev], beta0_cart[lev], 0, 0, 1, 1);
        }
    }

    // reset sig in the MLNodeLaplacian object so the call to updateVelocity
    // works properly
    for (int ilev = 0; ilev <= finest_level; ++ilev) {
        mlndlap.setSigma(ilev, sig[ilev]);
    }

    // When using the lambdabar term with a closed box we require the average of
    // lambda to be zero to remove the degeneracy in phi. Rescale phi here in
    // order to satisfy that constraint. The steps are as follows:
    // (1) Construct beta0 / (gamma1bar * p0) and take its radial sum, Bsum
    // (2) Construct the horizontal average of phi * beta0 / (gamma1bar * p0),
    //     then take its radial sum, Bphisum
    // (3) Add C = - Bphisum / Bsum onto phi
    // N.B. Pi and lambda are not updated in this step
    if ((proj_type == regular_timestep_comp || proj_type == pressure_iters_comp ) &&
        use_lambdabar_term && zero_average_lambda) {

        // set lambdabar = 1 / (gamma1bar p0)
        const auto pg = 0.25 * (p0_new + p0_old) * (gamma1bar_new + gamma1bar_old) ;
        lambdabar.setVal(1.0);
        lambdabar /= pg;

        // set time-centered beta0_nph_int = beta0 / (gamma1bar p0)
        // then construct the radial sum, Bsum
        auto beta0_nph_int = 0.5 * (beta0_old + beta0_new);
        beta0_nph_int *= lambdabar;
        Real Bsum = 0.;
        auto Bgp = beta0_nph_int.array();
        for (auto r = 0; r < base_geom.nr_fine; ++r) {
            Bsum += Bgp(0, r);
        }

        // copy phi to pi
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(pi[lev], phi[lev], 0, 0, 1, 0);
        }
        // snew(Pi) = beta0 * phi
        MakePiCC(beta0_cart);
        // beta0_nph_int = Average(beta0 * phi)
        Average(snew, beta0_nph_int, Pi);
        // beta0_nph_int = phi * beta0 / (gamma1bar * p0)
        beta0_nph_int *= lambdabar;

        // sum of beta0 phibar / (gamma1bar * p0)
        Real Bphisum = 0.;
        //auto Bphi = beta0_nph.array();
        for (auto r = 0; r < base_geom.nr_fine; ++r) {
            Bphisum += Bgp(0, r);
        }

        // correct phi
        Real C = - Bphisum / Bsum;
        for (int lev = 0; lev <= finest_level; ++lev) {
            phi[lev].plus(C,0);
        }
    }

    // compute a cell-centered grad(phi) from nodal phi
    Vector<MultiFab> gphi(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        gphi[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 0);
    }
    ComputeGradPhi(phi, gphi);

    // if using_alt_energy_fix, convert grad(phi) to beta0*grad(phi)
    if (using_alt_energy_fix) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                MultiFab::Multiply(gphi[lev], beta0_cart[lev], 0, dir, 1, 0);
            }
        }
    }

    // update pi and grad(pi)
    // initial_projection_comp: pi = 0         grad(pi) = 0
    // divu_iters_comp:         pi = 0         grad(pi) = 0
    // pressure_iters_comp:     pi = pi + phi  grad(pi) = grad(pi) + grad(phi)
    // regular_timestep_comp:   pi = phi/dt    grad(pi) = grad(phi)/dt
    if (is_predictor && proj_type == regular_timestep_comp) {
        // just update Pi
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(pi[lev], phi[lev], 0, 0, 1, 0);
            pi[lev].mult(1. / dt);
        }
    } else {
        if (proj_type == initial_projection_comp || proj_type == divu_iters_comp) {
            for (int lev = 0; lev <= finest_level; ++lev) {
                pi[lev].setVal(0.);
                gpi[lev].setVal(0.);
            }
        } else if (proj_type == pressure_iters_comp) {
            for (int lev = 0; lev <= finest_level; ++lev) {
                MultiFab::Add(pi[lev], phi[lev], 0, 0, 1, 0);
                MultiFab::Add(gpi[lev], gphi[lev], 0, 0, AMREX_SPACEDIM, 0);
            }
        } else if (proj_type == regular_timestep_comp) {
            for (int lev = 0; lev <= finest_level; ++lev) {
                MultiFab::Copy(pi[lev], phi[lev], 0, 0, 1, 0);
                MultiFab::Copy(gpi[lev], gphi[lev], 0, 0, AMREX_SPACEDIM, 0);
                pi[lev].mult(1. / dt);
                gpi[lev].mult(1. / dt);
            }
        }
    }

    // update pi, lambdabar
    if (proj_type == pressure_iters_comp ||
        proj_type == regular_timestep_comp) {

        // only update the dt*grav*lambdabar term in unew when performing an
        // interative update
        if (is_predictor && use_lambdabar_term && !spherical &&
            lambda_update_method >= 2) {
            Vector<MultiFab> grav_cart(finest_level + 1);
            Vector<MultiFab> lambdabar_cart(finest_level + 1);
            for (int lev = 0; lev <= finest_level; ++lev) {
                grav_cart[lev].define(grids[lev], dmap[lev],
                                      AMREX_SPACEDIM, 1);
                grav_cart[lev].setVal(0.);

                lambdabar_cart[lev].define(grids[lev], dmap[lev], 1, 1);
                lambdabar_cart[lev].setVal(0.);
            }
            // use the time-centered grav
            auto grav_cell_nph = 0.5 * dt * (grav_cell_old + grav_cell_new);
            Put1dArrayOnCart(grav_cell_nph, grav_cart, false, true, bcs_f, 0);
            Put1dArrayOnCart(lambdabar, lambdabar_cart, false, false, bcs_f, 0);

            // multiply grav by lambdabar, and add the result to unew
            for (int lev = 0; lev <= finest_level; ++lev) {
                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                    MultiFab::Multiply(grav_cart[lev], lambdabar_cart[lev], 0,
                                       dir, 1, 0);
                }
                MultiFab::Subtract(unew[lev], grav_cart[lev], 0, 0,
                              AMREX_SPACEDIM, 0);
                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                    MultiFab::Divide(grav_cart[lev], lambdabar_cart[lev], 0,
                                       dir, 1, 0);
                }
            }

            // average pi from nodes to cell-centers and store in the Pi component of snew
            MakePiCC(beta0_cart);

            // update lambdabar
            auto pg = 0.25 * (p0_new + p0_old) * (gamma1bar_new + gamma1bar_old);
            Average(snew, lambdabar, Pi);
            lambdabar /= pg;

            // multiply grav by lambdabar, and subtract the result from unew
            for (int lev = 0; lev <= finest_level; ++lev) {
                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                    MultiFab::Multiply(grav_cart[lev], lambdabar_cart[lev], 0,
                                       dir, 1, 0);
                }
                MultiFab::Add(unew[lev], grav_cart[lev], 0, 0,
                              AMREX_SPACEDIM, 0);
            }
        } else {

            // average pi from nodes to cell-centers and store in the Pi component of snew
            MakePiCC(beta0_cart);

            // update lambdabar
            if (use_lambdabar_term) {
                auto pg = 0.25 * (p0_new + p0_old) * (gamma1bar_new + gamma1bar_old);
                Average(snew, lambdabar, Pi);
                lambdabar /= pg;
            }
        }
    }

    // quit now if iterating
    if (is_predictor) return;

    // update velocity
    // initial_projection_comp: Utilde^0   = Vproj - sig*grad(phi)
    // divu_iters_comp:         Utilde^0   = Vproj - sig*grad(phi)
    // pressure_iters_comp:     Utilde^n+1 = Utilde^n + dt(Vproj-sig*grad(phi))
    // regular_timestep_comp:   Utilde^n+1 = Vproj - sig*grad(phi)
    if (proj_type == initial_projection_comp || proj_type == divu_iters_comp) {
        // Vproj = Vproj - sig*grad(phi)
        mlndlap.updateVelocity(amrex::GetVecOfPtrs(Vproj),
                               amrex::GetVecOfConstPtrs(phi));
        // Utilde^0   = Vproj - sig*grad(phi)
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(uold[lev], Vproj[lev], 0, 0, AMREX_SPACEDIM, 0);
        }
    } else if (proj_type == pressure_iters_comp) {
        // Vproj = Vproj - sig*grad(phi)
        // we do this manually instead of using mlndlap.updateVelocity() because
        // for alt_energy_fix we neet beta0*grad(phi)
        for (int lev = 0; lev <= finest_level; ++lev) {
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                MultiFab::Multiply(gphi[lev], sig[lev], 0, dir, 1, 0);
            }
            MultiFab::Subtract(Vproj[lev], gphi[lev], 0, 0, AMREX_SPACEDIM, 0);
        }
        // Utilde^n+1 = Utilde^n + dt(Vproj-sig*grad(phi))
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(unew[lev], Vproj[lev], 0, 0, AMREX_SPACEDIM, 0);
            unew[lev].mult(dt);
            MultiFab::Add(unew[lev], uold[lev], 0, 0, AMREX_SPACEDIM, 0);
        }
    } else if (proj_type == regular_timestep_comp) {
        // Vproj = Vproj - sig*grad(phi)
        // we do this manually instead of using mlndlap.updateVelocity() because
        // for alt_energy_fix we neet beta0*grad(phi)
        for (int lev = 0; lev <= finest_level; ++lev) {
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                MultiFab::Multiply(gphi[lev], sig[lev], 0, dir, 1, 0);
            }
            MultiFab::Subtract(Vproj[lev], gphi[lev], 0, 0, AMREX_SPACEDIM, 0);
        }
        // subtract dt*grav*lambdabar term in the one-step update
        if (use_lambdabar_term && !spherical && lambda_update_method == 1) {
            Vector<MultiFab> grav_cart(finest_level + 1);
            Vector<MultiFab> lambdabar_cart(finest_level + 1);
            for (int lev = 0; lev <= finest_level; ++lev) {
                grav_cart[lev].define(grids[lev], dmap[lev],
                                      AMREX_SPACEDIM, 1);
                grav_cart[lev].setVal(0.);

                lambdabar_cart[lev].define(grids[lev], dmap[lev], 1, 1);
                lambdabar_cart[lev].setVal(0.);
            }
            // use the time-centered grav
            auto grav_cell_nph = 0.5 * dt * (grav_cell_old + grav_cell_new);
            Put1dArrayOnCart(grav_cell_nph, grav_cart, false, true, bcs_f, 0);
            Put1dArrayOnCart(lambdabar, lambdabar_cart, false, false, bcs_f, 0);

            // multiply last grav component by dt*lambdabar, and add the result
            // to Vproj
            for (int lev = 0; lev <= finest_level; ++lev) {
                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                    MultiFab::Multiply(grav_cart[lev], lambdabar_cart[lev], 0,
                                       dir, 1, 0);
                }
                MultiFab::Subtract(Vproj[lev], grav_cart[lev], 0, 0,
                                   AMREX_SPACEDIM, 0);
            }
        }
        // Utilde^n+1 = Vproj - sig*grad(phi)
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(unew[lev], Vproj[lev], 0, 0, AMREX_SPACEDIM, 0);
        }
    }

    // average fine data onto coarser cells
    AverageDown(uold, 0, AMREX_SPACEDIM);
    AverageDown(unew, 0, AMREX_SPACEDIM);
    AverageDown(gpi, 0, AMREX_SPACEDIM);

    // fill ghost cells
    FillPatch(t_new, unew, unew, unew, 0, 0, AMREX_SPACEDIM, 0, bcs_u, 1);
    FillPatch(t_new, uold, uold, uold, 0, 0, AMREX_SPACEDIM, 0, bcs_u, 1);
}

// fill in Vproj
// initial_projection_comp: Utilde^0                        -- uold
// divu_iters_comp:         Utilde^0                        -- uold
// pressure_iters_comp:     (Utilde^n+1,* - Utilde^n)/dt    -- (unew-uold)/dt
// regular_timestep_comp:   (Utilde^n+1,* + dt*gpi/rhohalf) -- unew + dt*gpi/rhohalf
// note sig is only used for regular_timestep_comp, and currently holds rhohalf
void Maestro::CreateUvecForProj(int proj_type, Vector<MultiFab>& Vproj,
                                const Vector<MultiFab>& sig) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::CreateUvecForProj()", CreateUvecForProj);

    if (proj_type == initial_projection_comp || proj_type == divu_iters_comp) {
        // Vproj = uold
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(Vproj[lev], uold[lev], 0, 0, AMREX_SPACEDIM, 0);
        }
    } else if (proj_type == pressure_iters_comp) {
        // Vproj = (unew-uold)/dt
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(Vproj[lev], unew[lev], 0, 0, AMREX_SPACEDIM, 0);
            MultiFab::Saxpy(Vproj[lev], -1.0, uold[lev], 0, 0, AMREX_SPACEDIM,
                            0);
            Vproj[lev].mult(1 / dt);
        }

    } else if (proj_type == regular_timestep_comp) {
        // Vproj = unew + dt*gpi/rhohalf
        for (int lev = 0; lev <= finest_level; ++lev) {
            // Vproj = unew
            MultiFab::Copy(Vproj[lev], unew[lev], 0, 0, AMREX_SPACEDIM, 0);
            // convert gpi to gpi/rhohalf
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                MultiFab::Divide(gpi[lev], sig[lev], 0, dir, 1, 0);
            }
            // Vproj = unew + dt*gpi/rhohalf
            MultiFab::Saxpy(Vproj[lev], dt, gpi[lev], 0, 0, AMREX_SPACEDIM, 0);
            // revert gpi/rhohalf back to gpi
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                MultiFab::Multiply(gpi[lev], sig[lev], 0, dir, 1, 0);
            }
        }

        // subtract out dt*grav*lambdabar term in the one-step update
        if (use_lambdabar_term && !spherical && lambda_update_method == 1) {
            Vector<MultiFab> grav_cart(finest_level + 1);
            Vector<MultiFab> lambdabar_cart(finest_level + 1);
            for (int lev = 0; lev <= finest_level; ++lev) {
                grav_cart[lev].define(grids[lev], dmap[lev],
                                      AMREX_SPACEDIM, 1);
                grav_cart[lev].setVal(0.);

                lambdabar_cart[lev].define(grids[lev], dmap[lev], 1, 1);
                lambdabar_cart[lev].setVal(0.);
            }
            // use the time-centered grav
            auto grav_cell_nph = 0.5 * dt * (grav_cell_old + grav_cell_new);
            Put1dArrayOnCart(grav_cell_nph, grav_cart, false, true, bcs_f, 0);
            Put1dArrayOnCart(lambdabar, lambdabar_cart, false, false, bcs_f, 0);

            // multiply grav by lambdabar, and add the result to Vproj
            for (int lev = 0; lev <= finest_level; ++lev) {
                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                    MultiFab::Multiply(grav_cart[lev], lambdabar_cart[lev], 0,
                                       dir, 1, 0);
                }
                MultiFab::Add(Vproj[lev], grav_cart[lev], 0, 0,
                              AMREX_SPACEDIM, 0);
            }
        }
    } else {
        amrex::Abort("MaestroNodalProj: invalid proj_type");
    }

    Real time;
    if (proj_type == initial_projection_comp || proj_type == divu_iters_comp) {
        time = t_old;
    } else {
        time = t_new;
    }

    // fill ghost cells
    AverageDown(Vproj, 0, AMREX_SPACEDIM);
    FillPatch(time, Vproj, Vproj, Vproj, 0, 0, AMREX_SPACEDIM, 0, bcs_u, 1);
}

void Maestro::SetBoundaryVelocity(Vector<MultiFab>& vel) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::SetBoundaryVelocity()", SetBoundaryVelocity);

    // 1) At non-inflow faces, the normal component of velocity will be completely zero'd
    // 2) If a face is an inflow face, then the normal velocity at corners just outside inflow faces
    //                                will be zero'd outside of Neumann boundaries
    //                                (slipWall, noSlipWall, Symmetry)
    //                                BUT will retain non-zero values at periodic corners

    for (int lev = 0; lev <= finest_level; lev++) {
        const BoxArray& grids_lev = grids[lev];
        const Box& domainBox = geom[lev].Domain();

        for (int idir = 0; idir < BL_SPACEDIM; idir++) {
            if (phys_bc[idir] != Inflow &&
                phys_bc[AMREX_SPACEDIM + idir] != Inflow) {
                vel[lev].setBndry(0.0, idir, 1);
            } else {
#ifdef _OPENMP
#pragma omp parallel
#endif
                for (MFIter mfi(vel[lev]); mfi.isValid(); ++mfi) {
                    int i = mfi.index();
                    FArrayBox& v_fab = (vel[lev])[mfi];

                    const Box& reg = grids_lev[i];
                    const Box& bxg1 = amrex::grow(reg, 1);

                    BoxList bxlist(reg);

                    if (phys_bc[idir] == Inflow &&
                        reg.smallEnd(idir) == domainBox.smallEnd(idir)) {
                        Box bx;  // bx is the region we *protect* from zero'ing

                        bx = amrex::adjCellLo(reg, idir);

                        for (int odir = 0; odir < BL_SPACEDIM; odir++) {
                            if (odir != idir) {
                                if (geom[lev].isPeriodic(odir)) {
                                    bx.grow(odir, 1);
                                }
                                if (reg.bigEnd(odir) !=
                                    domainBox.bigEnd(odir)) {
                                    bx.growHi(odir, 1);
                                }
                                if (reg.smallEnd(odir) !=
                                    domainBox.smallEnd(odir)) {
                                    bx.growLo(odir, 1);
                                }
                            }
                        }
                        bxlist.push_back(bx);
                    }

                    if (phys_bc[AMREX_SPACEDIM + idir] == Inflow &&
                        reg.bigEnd(idir) == domainBox.bigEnd(idir)) {
                        Box bx;  // bx is the region we *protect* from zero'ing

                        bx = amrex::adjCellHi(reg, idir);

                        for (int odir = 0; odir < BL_SPACEDIM; odir++) {
                            if (odir != idir) {
                                if (geom[lev].isPeriodic(odir)) {
                                    bx.grow(odir, 1);
                                }
                                if (reg.bigEnd(odir) !=
                                    domainBox.bigEnd(odir)) {
                                    bx.growHi(odir, 1);
                                }
                                if (reg.smallEnd(odir) !=
                                    domainBox.smallEnd(odir)) {
                                    bx.growLo(odir, 1);
                                }
                            }
                        }
                        bxlist.push_back(bx);
                    }

                    BoxList bxlist2 = amrex::complementIn(bxg1, bxlist);

                    for (auto it : bxlist2) {
                        Box ovlp = it & v_fab.box();
                        if (ovlp.ok()) {
                            v_fab.setVal<RunOn::Device>(0.0, ovlp, idir, 1);
                        }
                    }
                }  // end loop over grids
            }      // end if/else logic for inflow
        }          // end loop over direction
    }              // end loop over levels
}

// given a nodal phi, compute grad(phi) at cell centers
void Maestro::ComputeGradPhi(Vector<MultiFab>& phi, Vector<MultiFab>& gphi) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::ComputeGradPhi()", ComputeGradPhi);

    for (int lev = 0; lev <= finest_level; ++lev) {
        MultiFab& gphi_mf = gphi[lev];

        const auto dx = geom[lev].CellSizeArray();

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(gphi_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the tile's valid region
            const Box& tilebox = mfi.tilebox();

            const Array4<const Real> phi_arr = phi[lev].array(mfi);
            const Array4<Real> gphi_arr = gphi[lev].array(mfi);

#if (AMREX_SPACEDIM == 2)
            ParallelFor(tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                gphi_arr(i, j, k, 0) =
                    0.5 *
                    (phi_arr(i + 1, j, k) + phi_arr(i + 1, j + 1, k) -
                     phi_arr(i, j, k) - phi_arr(i, j + 1, k)) /
                    dx[0];
                gphi_arr(i, j, k, 1) =
                    0.5 *
                    (phi_arr(i, j + 1, k) + phi_arr(i + 1, j + 1, k) -
                     phi_arr(i, j, k) - phi_arr(i + 1, j, k)) /
                    dx[1];
            });
#else
            ParallelFor(tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                gphi_arr(i, j, k, 0) =
                    0.25 *
                    (phi_arr(i + 1, j, k) + phi_arr(i + 1, j + 1, k) +
                     phi_arr(i + 1, j, k + 1) + phi_arr(i + 1, j + 1, k + 1) -
                     phi_arr(i, j, k) - phi_arr(i, j + 1, k) -
                     phi_arr(i, j, k + 1) - phi_arr(i, j + 1, k + 1)) /
                    dx[0];
                gphi_arr(i, j, k, 1) =
                    0.25 *
                    (phi_arr(i, j + 1, k) + phi_arr(i + 1, j + 1, k) +
                     phi_arr(i, j + 1, k + 1) + phi_arr(i + 1, j + 1, k + 1) -
                     phi_arr(i, j, k) - phi_arr(i + 1, j, k) -
                     phi_arr(i, j, k + 1) - phi_arr(i + 1, j, k + 1)) /
                    dx[1];
                gphi_arr(i, j, k, 2) =
                    0.25 *
                    (phi_arr(i, j, k + 1) + phi_arr(i + 1, j, k + 1) +
                     phi_arr(i, j + 1, k + 1) + phi_arr(i + 1, j + 1, k + 1) -
                     phi_arr(i, j, k) - phi_arr(i + 1, j, k) -
                     phi_arr(i, j + 1, k) - phi_arr(i + 1, j + 1, k)) /
                    dx[2];
            });
#endif
        }
    }
}

// average nodal pi to cell-centers and put in the Pi component of snew
void Maestro::MakePiCC(const Vector<MultiFab>& beta0_cart) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakePiCC()", MakePiCC);

    const bool use_alt_energy_fix_loc = use_alt_energy_fix;

    for (int lev = 0; lev <= finest_level; ++lev) {
        MultiFab& snew_mf = snew[lev];

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(snew_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the tile's valid region
            const Box& tilebox = mfi.tilebox();

            const Array4<const Real> pi_arr = pi[lev].array(mfi);
            const Array4<Real> pi_cc = snew[lev].array(mfi, Pi);
            const Array4<const Real> beta0_cart_arr =
                beta0_cart[lev].array(mfi);

            ParallelFor(tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
#if (AMREX_SPACEDIM == 2)
                pi_cc(i, j, k) =
                    0.25 * (pi_arr(i, j, k) + pi_arr(i + 1, j, k) +
                            pi_arr(i, j + 1, k) + pi_arr(i + 1, j + 1, k));
#else
                pi_cc(i,j,k) = 0.125 * (pi_arr(i,j,k) + pi_arr(i+1,j,k)
                    + pi_arr(i,j+1,k) + pi_arr(i,j,k+1)
                    + pi_arr(i+1,j+1,k) + pi_arr(i+1,j,k+1)
                    + pi_arr(i,j+1,k+1) + pi_arr(i+1,j+1,k+1));
#endif
                if (use_alt_energy_fix_loc) {
                    pi_cc(i, j, k) = pi_cc(i, j, k) * beta0_cart_arr(i, j, k);
                }
            });
        }
    }
}

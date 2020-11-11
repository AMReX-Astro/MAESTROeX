
#include <AMReX_VisMF.H>
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

// umac enters with face-centered, time-centered Utilde^* and should leave with Utilde
// macphi is the solution to the elliptic solve and
//   enters as either zero, or the solution to the predictor MAC projection
// macrhs enters as beta0*(S-Sbar)
// beta0 is a 1d cell-centered array
void Maestro::MacProj(Vector<std::array<MultiFab, AMREX_SPACEDIM> >& umac,
                      Vector<MultiFab>& macphi, const Vector<MultiFab>& macrhs,
                      const BaseState<Real>& beta0, const bool is_predictor) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MacProj()", MacProj);

    // this will hold solver RHS = macrhs - div(beta0*umac)
    Vector<MultiFab> solverrhs(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        solverrhs[lev].define(grids[lev], dmap[lev], 1, 0);
    }

    // we also need beta0 at edges
    // allocate AND compute it here
    BaseState<Real> beta0_edge(base_geom.max_radial_level + 1,
                               base_geom.nr_fine + 1);

    Vector<std::array<MultiFab, AMREX_SPACEDIM> > beta0_cart_edge(finest_level +
                                                                  1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        AMREX_D_TERM(
            beta0_cart_edge[lev][0].define(convert(grids[lev], nodal_flag_x),
                                           dmap[lev], 1, 1);
            , beta0_cart_edge[lev][1].define(convert(grids[lev], nodal_flag_y),
                                             dmap[lev], 1, 1);
            , beta0_cart_edge[lev][2].define(convert(grids[lev], nodal_flag_z),
                                             dmap[lev], 1, 1););
    }

    if (spherical) {
#if (AMREX_SPACEDIM == 3)
        MakeS0mac(beta0, beta0_cart_edge);
#endif
    } else {
        CelltoEdge(beta0, beta0_edge);
    }

    // convert Utilde^* to beta0*Utilde^*
    int mult_or_div;
    if (!spherical) {
        mult_or_div = 1;
        MultFacesByBeta0(umac, beta0, beta0_edge, mult_or_div);
    } else {  // spherical
        for (int lev = 0; lev <= finest_level; ++lev) {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                MultiFab::Multiply(umac[lev][idim], beta0_cart_edge[lev][idim],
                                   0, 0, 1, 0);
            }
        }
    }

    // compute the RHS for the solve, RHS = macrhs - div(beta0*umac)
    AverageDownFaces(umac);
    ComputeMACSolverRHS(solverrhs, macrhs, umac);

    // create a MultiFab filled with rho and 1 ghost cell.
    // if this is the predictor mac projection, use rho^n
    // if this is the corrector mac projection, use (1/2)(rho^n + rho^{n+1,*})
    Vector<MultiFab> rho(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        rho[lev].define(grids[lev], dmap[lev], 1, 1);
        // needed to avoid NaNs in filling corner ghost cells with 2 physical boundaries
        rho[lev].setVal(0.);
    }
    Real rho_time = (is_predictor) ? t_old : 0.5 * (t_old + t_new);
    FillPatch(rho_time, rho, sold, snew, Rho, 0, 1, Rho, bcs_s);

    // coefficients for solver
    Vector<MultiFab> acoef(finest_level + 1);
    Vector<std::array<MultiFab, AMREX_SPACEDIM> > face_bcoef(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        acoef[lev].define(grids[lev], dmap[lev], 1, 0);
        AMREX_D_TERM(
            face_bcoef[lev][0].define(convert(grids[lev], nodal_flag_x),
                                      dmap[lev], 1, 0);
            , face_bcoef[lev][1].define(convert(grids[lev], nodal_flag_y),
                                        dmap[lev], 1, 0);
            , face_bcoef[lev][2].define(convert(grids[lev], nodal_flag_z),
                                        dmap[lev], 1, 0););
    }

    // set cell-centered A coefficient to zero
    for (int lev = 0; lev <= finest_level; ++lev) {
        acoef[lev].setVal(0.);
    }

    // OR 1) average face-centered B coefficients to rho
    for (int lev = 0; lev <= finest_level; ++lev) {
        amrex::average_cellcenter_to_face(GetArrOfPtrs(face_bcoef[lev]),
                                          rho[lev], geom[lev]);
    }

    // AND 2) invert B coefficients to 1/rho
    for (int lev = 0; lev <= finest_level; ++lev) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            face_bcoef[lev][idim].invert(1.0, 0, 1);
        }
    }

    // Make sure that the fine edges average down onto the coarse edges (edge_restriction)
    AverageDownFaces(face_bcoef);

    // multiply face-centered B coefficients by beta0 so they contain beta0/rho
    if (!spherical) {
        mult_or_div = 1;
        MultFacesByBeta0(face_bcoef, beta0, beta0_edge, mult_or_div);
        if (use_alt_energy_fix) {
            MultFacesByBeta0(face_bcoef, beta0, beta0_edge, mult_or_div);
        }
    } else {  //spherical
        for (int lev = 0; lev <= finest_level; ++lev) {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                MultiFab::Multiply(face_bcoef[lev][idim],
                                   beta0_cart_edge[lev][idim], 0, 0, 1, 0);
                if (use_alt_energy_fix) {
                    MultiFab::Multiply(face_bcoef[lev][idim],
                                       beta0_cart_edge[lev][idim], 0, 0, 1, 0);
                }
            }
        }
    }

    // Set up implicit solve using MLABecLaplacian class
    //
    LPInfo info;
    info.setMetricTerm(false);

    if (mg_bottom_solver == 4) {
        info.setAgglomeration(true);
        info.setConsolidation(true);
    } else {
        info.setAgglomeration(false);
        info.setConsolidation(false);
    }

    // Only pass up to defined level to prevent looping over undefined grids.
    MLABecLaplacian mlabec(Geom(0, finest_level), grids, dmap, info);

    // order of stencil
    int linop_maxorder = 2;
    mlabec.setMaxOrder(linop_maxorder);

    // set boundaries for mlabec using velocity bc's
    SetMacSolverBCs(mlabec);

    for (int lev = 0; lev <= finest_level; ++lev) {
        mlabec.setLevelBC(lev, &macphi[lev]);
    }

    mlabec.setScalars(0.0, 1.0);

    for (int lev = 0; lev <= finest_level; ++lev) {
        mlabec.setACoeffs(lev, acoef[lev]);
        mlabec.setBCoeffs(lev, amrex::GetArrOfConstPtrs(face_bcoef[lev]));
    }

    // solve -div B grad phi = RHS

    // build an MLMG solver
    MLMG mac_mlmg(mlabec);

    // set solver parameters
    mac_mlmg.setVerbose(mg_verbose);
    mac_mlmg.setBottomVerbose(cg_verbose);

    // tolerance parameters taken from original MAESTRO fortran code
    const Real mac_tol_abs = -1.e0;
    const Real mac_tol_rel =
        amrex::min(eps_mac * pow(mac_level_factor, finest_level), eps_mac_max);

    // solve for phi
    mac_mlmg.solve(GetVecOfPtrs(macphi), GetVecOfConstPtrs(solverrhs),
                   mac_tol_rel, mac_tol_abs);

    // update velocity, beta0 * Utilde = beta0 * Utilde^* - B grad phi

    // storage for "-B grad_phi"
    Vector<std::array<MultiFab, AMREX_SPACEDIM> > mac_fluxes(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        AMREX_D_TERM(
            mac_fluxes[lev][0].define(convert(grids[lev], nodal_flag_x),
                                      dmap[lev], 1, 0);
            , mac_fluxes[lev][1].define(convert(grids[lev], nodal_flag_y),
                                        dmap[lev], 1, 0);
            , mac_fluxes[lev][2].define(convert(grids[lev], nodal_flag_z),
                                        dmap[lev], 1, 0););
    }

    Vector<std::array<MultiFab*, AMREX_SPACEDIM> > mac_fluxptr(finest_level +
                                                               1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        // fluxes computed are "-B grad phi"
        mac_fluxptr[lev] = GetArrOfPtrs(mac_fluxes[lev]);
    }
    mac_mlmg.getFluxes(mac_fluxptr);

    // Make sure that the fine edges average down onto the coarse edges (edge_restriction)
    AverageDownFaces(mac_fluxes);

    for (int lev = 0; lev <= finest_level; ++lev) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            // add -B grad phi to beta0*Utilde
            MultiFab::Add(umac[lev][idim], mac_fluxes[lev][idim], 0, 0, 1, 0);
        }
    }

    // convert beta0*Utilde to Utilde
    if (!spherical) {
        mult_or_div = 0;
        MultFacesByBeta0(umac, beta0, beta0_edge, mult_or_div);
    } else {
        for (int lev = 0; lev <= finest_level; ++lev) {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                MultiFab::Divide(umac[lev][idim], beta0_cart_edge[lev][idim], 0,
                                 0, 1, 0);
            }
        }
    }

    if (finest_level == 0) {
        // fill periodic ghost cells
        for (int lev = 0; lev <= finest_level; ++lev) {
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                umac[lev][d].FillBoundary(geom[lev].periodicity());
            }
        }

        // fill ghost cells behind physical boundaries
        FillUmacGhost(umac);
    } else {
        // edge_restriction for velocities
        AverageDownFaces(umac);

        // fill level n ghost cells using interpolation from level n-1 data
        FillPatchUedge(umac);
    }
}

// multiply (or divide) face-data by beta0
void Maestro::MultFacesByBeta0(
    Vector<std::array<MultiFab, AMREX_SPACEDIM> >& edge,
    const BaseState<Real>& beta0_s, const BaseState<Real>& beta0_edge_s,
    const int& mult_or_div) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MultFacesByBeta0()", MultFacesByBeta0);

    const auto beta0 = beta0_s.const_array();
    const auto beta0_edge = beta0_edge_s.const_array();

    // write an MFIter loop to convert edge -> beta0*edge OR beta0*edge -> edge
    for (int lev = 0; lev <= finest_level; ++lev) {
        // loop over boxes
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(sold[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of valid region
            const Box& xbx = mfi.nodaltilebox(0);
            const Box& ybx = mfi.nodaltilebox(1);
#if (AMREX_SPACEDIM == 3)
            const Box& zbx = mfi.nodaltilebox(2);
#endif

            const Array4<Real> uedge = edge[lev][0].array(mfi);
            const Array4<Real> vedge = edge[lev][1].array(mfi);
#if (AMREX_SPACEDIM == 3)
            const Array4<Real> wedge = edge[lev][2].array(mfi);
#endif

            if (mult_or_div == 1) {
                ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
#if (AMREX_SPACEDIM == 2)
                    int r = j;
#else 
                    int r = k;
#endif
                    uedge(i, j, k) *= beta0(lev, r);
                });

                ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
#if (AMREX_SPACEDIM == 2)
                    vedge(i, j, k) *= beta0_edge(lev, j);
#else 
                    vedge(i,j,k) *= beta0(lev,k);
#endif
                });

#if (AMREX_SPACEDIM == 3)
                ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    wedge(i, j, k) *= beta0_edge(lev, k);
                });
#endif
            } else {
                ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
#if (AMREX_SPACEDIM == 2)
                    int r = j;
#else 
                    int r = k;
#endif
                    uedge(i, j, k) /= beta0(lev, r);
                });

                ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
#if (AMREX_SPACEDIM == 2)
                    vedge(i, j, k) /= beta0_edge(lev, j);
#else 
                    vedge(i,j,k) /= beta0(lev,k);
#endif
                });

#if (AMREX_SPACEDIM == 3)
                ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    wedge(i, j, k) /= beta0_edge(lev, k);
                });
#endif
            }
        }
    }
}

// compute the RHS for the solve, RHS = macrhs - div(beta0*umac)
void Maestro::ComputeMACSolverRHS(
    Vector<MultiFab>& solverrhs, const Vector<MultiFab>& macrhs,
    const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& umac) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::ComputeMACSolverRHS()", ComputeMACSolverRHS);

    // Note that umac = beta0*mac
    for (int lev = 0; lev <= finest_level; ++lev) {
        // loop over boxes
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(solverrhs[lev], TilingIfNotGPU()); mfi.isValid();
             ++mfi) {
            // Get the index space of valid region
            const Box& tileBox = mfi.tilebox();

            const auto dx = geom[lev].CellSizeArray();

            const Array4<Real> solverrhs_arr = solverrhs[lev].array(mfi);
            const Array4<const Real> macrhs_arr = macrhs[lev].array(mfi);
            const Array4<const Real> uedge = umac[lev][0].array(mfi);
            const Array4<const Real> vedge = umac[lev][1].array(mfi);
#if (AMREX_SPACEDIM == 3)
            const Array4<const Real> wedge = umac[lev][2].array(mfi);
#endif

            ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                // Compute newrhs = oldrhs - div(Uedge)
                solverrhs_arr(i, j, k) =
                    macrhs_arr(i, j, k) -
                    ((uedge(i + 1, j, k) - uedge(i, j, k)) / dx[0] +
                     (vedge(i, j + 1, k) - vedge(i, j, k)) / dx[1]
#if (AMREX_SPACEDIM == 3)
                     + (wedge(i, j, k + 1) - wedge(i, j, k)) / dx[2]
#endif
                    );
            });
        }
    }
}

// Average bcoefs at faces using inverse of rho
void Maestro::AvgFaceBcoeffsInv(
    Vector<std::array<MultiFab, AMREX_SPACEDIM> >& facebcoef,
    const Vector<MultiFab>& rhocc) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::AvgFaceBcoeffsInv()", AvgFaceBcoeffsInv);

    // write an MFIter loop
    for (int lev = 0; lev <= finest_level; ++lev) {
        // loop over boxes
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(rhocc[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of valid region
            const Box& tileBox = mfi.tilebox();
            const Box& xbx = amrex::growHi(tileBox, 0, 1);
            const Box& ybx = amrex::growHi(tileBox, 1, 1);
#if (AMREX_SPACEDIM == 3)
            const Box& zbx = amrex::growHi(tileBox, 2, 1);
#endif

            const Array4<Real> xbcoef = facebcoef[lev][0].array(mfi);
            const Array4<Real> ybcoef = facebcoef[lev][1].array(mfi);
#if (AMREX_SPACEDIM == 3)
            const Array4<Real> zbcoef = facebcoef[lev][2].array(mfi);
#endif
            const Array4<const Real> rhocc_arr = rhocc[lev].array(mfi);

            ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                xbcoef(i, j, k) =
                    2.0 / (rhocc_arr(i, j, k) + rhocc_arr(i - 1, j, k));
            });

            ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                ybcoef(i, j, k) =
                    2.0 / (rhocc_arr(i, j, k) + rhocc_arr(i, j - 1, k));
            });
#if (AMREX_SPACEDIM == 3)
            ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                zbcoef(i, j, k) =
                    2.0 / (rhocc_arr(i, j, k) + rhocc_arr(i, j, k - 1));
            });
#endif
        }
    }
}

// Set boundaries for MAC velocities
void Maestro::SetMacSolverBCs(MLABecLaplacian& mlabec) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::SetMacSolverBCs()", SetMacSolverBCs);

    // build array of boundary conditions needed by MLABecLaplacian
    std::array<LinOpBCType, AMREX_SPACEDIM> mlmg_lobc;
    std::array<LinOpBCType, AMREX_SPACEDIM> mlmg_hibc;

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (Geom(0).isPeriodic(idim)) {
            mlmg_lobc[idim] = mlmg_hibc[idim] = LinOpBCType::Periodic;
        } else {
            // lo-side BCs
            if (phys_bc[idim] == Outflow) {
                mlmg_lobc[idim] = LinOpBCType::Dirichlet;
            } else {
                mlmg_lobc[idim] = LinOpBCType::Neumann;
            }

            // hi-side BCs
            if (phys_bc[AMREX_SPACEDIM + idim] == Outflow) {
                mlmg_hibc[idim] = LinOpBCType::Dirichlet;
            } else {
                mlmg_hibc[idim] = LinOpBCType::Neumann;
            }
        }
    }

    mlabec.setDomainBC(mlmg_lobc, mlmg_hibc);
}

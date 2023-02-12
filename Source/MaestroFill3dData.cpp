/*
   A collection of routines for mapping 1D arrays onto multi-D cartesian MultiFabs
 */

#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void Maestro::Put1dArrayOnCart(const BaseState<Real>& s0,
                               Vector<MultiFab>& s0_cart,
                               const bool is_input_edge_centered,
                               const bool is_output_a_vector,
                               const Vector<BCRec>& bcs, const int sbccomp,
                               const int variable_type) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Put1dArrayOnCart()", Put1dArrayOnCart);

    int ng = s0_cart[0].nGrow();
    if (ng > 0 && bcs.empty()) {
        Abort("Put1dArrayOnCart with ghost cells requires bcs input");
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
        Put1dArrayOnCart(lev, s0, s0_cart, is_input_edge_centered,
                         is_output_a_vector, bcs, sbccomp);
    }

    int ncomp = is_output_a_vector ? AMREX_SPACEDIM : 1;

    // set covered coarse cells to be the average of overlying fine cells
    AverageDown(s0_cart, 0, ncomp);

    // fill ghost cells using first-order extrapolation
    if (ng > 0) {
        FillPatch(t_old, s0_cart, s0_cart, s0_cart, 0, 0, ncomp, sbccomp, bcs,
                  variable_type);
    }
}

void Maestro::Put1dArrayOnCart(const int lev, const BaseState<Real>& s0,
                               Vector<MultiFab>& s0_cart,
                               const bool is_input_edge_centered,
                               const bool is_output_a_vector,
                               const Vector<BCRec>& bcs, const int sbccomp) {
    Put1dArrayOnCart(lev, s0, s0_cart[lev], is_input_edge_centered,
                     is_output_a_vector, bcs, sbccomp);
}

void Maestro::Put1dArrayOnCart(const int lev, const BaseState<Real>& s0,
                               MultiFab& s0_cart,
                               const bool is_input_edge_centered,
                               const bool is_output_a_vector,
                               [[maybe_unused]] const Vector<BCRec>& bcs,
                               [[maybe_unused]] const int sbccomp) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Put1dArrayOnCart_lev()", Put1dArrayOnCart);

    const auto dx = geom[lev].CellSizeArray();
    const auto prob_lo = geom[lev].ProbLoArray();
    const auto& center_p = center;

    const auto& r_edge_loc = base_geom.r_edge_loc;
    const auto& r_cc_loc = base_geom.r_cc_loc;
    const auto s0_arr = s0.const_array();

    const int nr_fine = base_geom.nr_fine;
    const int w0_interp_type_loc = w0_interp_type;

    // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(s0_cart, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        // Get the index space of the valid region
        const Box& tileBox = mfi.tilebox();

        const Array4<Real> s0_cart_arr = s0_cart.array(mfi);

        if (!spherical) {
            const int outcomp = is_output_a_vector ? AMREX_SPACEDIM - 1 : 0;

            ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                const int r = AMREX_SPACEDIM == 2 ? j : k;

                s0_cart_arr(i, j, k, outcomp) =
                    is_input_edge_centered
                        ? 0.5 * (s0_arr(lev, r) + s0_arr(lev, r + 1))
                        : s0_arr(lev, r);
            });

        } else {
            const Array4<const int> cc_to_r = cell_cc_to_r[lev].array(mfi);

            if (use_exact_base_state) {
                if (is_input_edge_centered) {
                    // we implemented three different ideas for computing s0_cart,
                    // where s0 is edge-centered.
                    // 1.  Piecewise constant
                    // 2.  Piecewise linear
                    // 3.  Quadratic

                    ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j,
                                                              int k) {
                        Real x =
                            prob_lo[0] + (Real(i) + 0.5) * dx[0] - center_p[0];
                        Real y =
                            prob_lo[1] + (Real(j) + 0.5) * dx[1] - center_p[1];
                        Real z =
                            prob_lo[2] + (Real(k) + 0.5) * dx[2] - center_p[2];

                        Real radius = sqrt(x * x + y * y + z * z);
                        int index = cc_to_r(i, j, k);

                        Real rfac;
                        if (index < nr_fine) {
                            rfac =
                                (radius - r_edge_loc(0, index + 1)) /
                                (r_cc_loc(0, index + 1) - r_cc_loc(0, index));
                        } else {
                            rfac =
                                (radius - r_edge_loc(0, index + 1)) /
                                (r_cc_loc(0, index) - r_cc_loc(0, index - 1));
                        }

                        Real s0_cart_val = 0.0;

                        if (w0_interp_type_loc == 1) {
                            s0_cart_val = rfac > 0.5 ? s0_arr(0, index + 1)
                                                     : s0_arr(0, index);

                        } else if (w0_interp_type_loc == 2) {
                            if (index < nr_fine) {
                                s0_cart_val = rfac * s0_arr(0, index + 1) +
                                              (1.0 - rfac) * s0_arr(0, index);
                            } else {
                                s0_cart_val = s0_arr(0, nr_fine);
                            }

                        } else if (w0_interp_type_loc == 3) {
                            if (index <= 0) {
                                index = 0;
                            } else if (index >= nr_fine - 1) {
                                index = nr_fine - 2;
                            } else if (radius - r_edge_loc(0, index) <
                                       r_edge_loc(0, index + 1)) {
                                index--;
                            }

                            s0_cart_val = QuadInterp(
                                radius, r_edge_loc(0, index),
                                r_edge_loc(0, index + 1),
                                r_edge_loc(0, index + 2), s0_arr(0, index),
                                s0_arr(0, index + 1), s0_arr(0, index + 2));
                        }

                        if (is_output_a_vector) {
                            s0_cart_arr(i, j, k, 0) = s0_cart_val * x / radius;
                            s0_cart_arr(i, j, k, 1) = s0_cart_val * y / radius;
                            s0_cart_arr(i, j, k, 2) = s0_cart_val * z / radius;
                        } else {
                            s0_cart_arr(i, j, k, 0) = s0_cart_val;
                        }
                    });

                } else {  // is_input_edge_centered = 0
                    // we directly inject the spherical values into each cell center
                    // because s0 is also bin-centered.

                    ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j,
                                                              int k) {
                        Real x =
                            prob_lo[0] + (Real(i) + 0.5) * dx[0] - center_p[0];
                        Real y =
                            prob_lo[1] + (Real(j) + 0.5) * dx[1] - center_p[1];
                        Real z =
                            prob_lo[2] + (Real(k) + 0.5) * dx[2] - center_p[2];

                        Real radius = sqrt(x * x + y * y + z * z);
                        int index = cc_to_r(i, j, k);

                        Real s0_cart_val = s0_arr(0, index);

                        if (is_output_a_vector) {
                            s0_cart_arr(i, j, k, 0) = s0_cart_val * x / radius;
                            s0_cart_arr(i, j, k, 1) = s0_cart_val * y / radius;
                            s0_cart_arr(i, j, k, 2) = s0_cart_val * z / radius;
                        } else {
                            s0_cart_arr(i, j, k, 0) = s0_cart_val;
                        }
                    });
                }  // is_input_edge_centered

            } else {  // use_exact_base_state = 0

                const Real drf = base_geom.dr_fine;

                if (is_input_edge_centered) {
                    // we implemented three different ideas for computing s0_cart,
                    // where s0 is edge-centered.
                    // 1.  Piecewise constant
                    // 2.  Piecewise linear
                    // 3.  Quadratic

                    ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j,
                                                              int k) {
                        Real x =
                            prob_lo[0] + (Real(i) + 0.5) * dx[0] - center_p[0];
                        Real y =
                            prob_lo[1] + (Real(j) + 0.5) * dx[1] - center_p[1];
                        Real z =
                            prob_lo[2] + (Real(k) + 0.5) * dx[2] - center_p[2];

                        Real radius = sqrt(x * x + y * y + z * z);

                        auto index = int(radius / drf);
                        Real rfac = (radius - Real(index) * drf) / drf;
                        Real s0_cart_val = 0.0;

                        if (w0_interp_type_loc == 1) {
                            s0_cart_val = rfac > 0.5 ? s0_arr(0, index + 1)
                                                     : s0_arr(0, index);

                        } else if (w0_interp_type_loc == 2) {
                            if (index < nr_fine) {
                                s0_cart_val = rfac * s0_arr(0, index + 1) +
                                              (1.0 - rfac) * s0_arr(0, index);
                            } else {
                                s0_cart_val = s0_arr(0, nr_fine);
                            }

                        } else if (w0_interp_type_loc == 3) {
                            if (index <= 0) {
                                index = 0;
                            } else if (index >= nr_fine - 1) {
                                index = nr_fine - 2;
                            } else if (radius - r_edge_loc(0, index) <
                                       r_edge_loc(0, index + 1)) {
                                index--;
                            }

                            s0_cart_val = QuadInterp(
                                radius, r_edge_loc(0, index),
                                r_edge_loc(0, index + 1),
                                r_edge_loc(0, index + 2), s0_arr(0, index),
                                s0_arr(0, index + 1), s0_arr(0, index + 2));
                        }

                        if (is_output_a_vector) {
                            s0_cart_arr(i, j, k, 0) = s0_cart_val * x / radius;
                            s0_cart_arr(i, j, k, 1) = s0_cart_val * y / radius;
                            s0_cart_arr(i, j, k, 2) = s0_cart_val * z / radius;
                        } else {
                            s0_cart_arr(i, j, k, 0) = s0_cart_val;
                        }
                    });

                } else {  // is_input_edge_centered = 0

                    const int s0_interp_type_loc = s0_interp_type;

                    // we currently have three different ideas for computing s0_cart,
                    // where s0 is bin-centered.
                    // 1.  Piecewise constant
                    // 2.  Piecewise linear
                    // 3.  Quadratic
                    ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j,
                                                              int k) {
                        Real x =
                            prob_lo[0] + (Real(i) + 0.5) * dx[0] - center_p[0];
                        Real y =
                            prob_lo[1] + (Real(j) + 0.5) * dx[1] - center_p[1];
                        Real z =
                            prob_lo[2] + (Real(k) + 0.5) * dx[2] - center_p[2];

                        Real radius = sqrt(x * x + y * y + z * z);
                        auto index = int(radius / drf);

                        Real s0_cart_val = 0.0;

                        if (s0_interp_type_loc == 1) {
                            s0_cart_val = s0_arr(0, index);

                        } else if (s0_interp_type_loc == 2) {
                            if (radius >= r_cc_loc(0, index)) {
                                if (index >= nr_fine - 1) {
                                    s0_cart_val = s0_arr(0, nr_fine - 1);
                                } else {
                                    s0_cart_val =
                                        s0_arr(0, index + 1) *
                                            (radius - r_cc_loc(0, index)) /
                                            drf +
                                        s0_arr(0, index) *
                                            (r_cc_loc(0, index + 1) - radius) /
                                            drf;
                                }
                            } else {
                                if (index == 0) {
                                    s0_cart_val = s0_arr(0, index);
                                } else if (index > nr_fine - 1) {
                                    s0_cart_val = s0_arr(0, nr_fine - 1);
                                } else {
                                    s0_cart_val =
                                        s0_arr(0, index) *
                                            (radius - r_cc_loc(0, index - 1)) /
                                            drf +
                                        s0_arr(0, index - 1) *
                                            (r_cc_loc(0, index) - radius) / drf;
                                }
                            }
                        } else if (s0_interp_type_loc == 3) {
                            if (index == 0) {
                                index = 1;
                            } else if (index >= nr_fine - 1) {
                                index = nr_fine - 2;
                            }

                            s0_cart_val = QuadInterp(
                                radius, r_cc_loc(0, index - 1),
                                r_cc_loc(0, index), r_cc_loc(0, index + 1),
                                s0_arr(0, index - 1), s0_arr(0, index),
                                s0_arr(0, index + 1));
                        }

                        if (is_output_a_vector) {
                            s0_cart_arr(i, j, k, 0) = s0_cart_val * x / radius;
                            s0_cart_arr(i, j, k, 1) = s0_cart_val * y / radius;
                            s0_cart_arr(i, j, k, 2) = s0_cart_val * z / radius;
                        } else {
                            s0_cart_arr(i, j, k, 0) = s0_cart_val;
                        }
                    });
                }  // is_input_edge_centered
            }      // use_exact_base_state
        }
    }
}

AMREX_GPU_DEVICE
Real QuadInterp(const Real x, const Real x0, const Real x1, const Real x2,
                const Real y0, const Real y1, const Real y2, bool limit) {
    Real y = y0 + (y1 - y0) / (x1 - x0) * (x - x0) +
             ((y2 - y1) / (x2 - x1) - (y1 - y0) / (x1 - x0)) / (x2 - x0) *
                 (x - x0) * (x - x1);

    if (limit) {
        if (y > amrex::max(y0, amrex::max(y1, y2))) {
            y = amrex::max(y0, amrex::max(y1, y2));
        }
        if (y < amrex::min(y0, amrex::min(y1, y2))) {
            y = amrex::min(y0, amrex::min(y1, y2));
        }
    }

    return y;
}

void Maestro::Addw0(Vector<std::array<MultiFab, AMREX_SPACEDIM> >& u_edge,
                    [[maybe_unused]] const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& w0mac,
                    const Real& mult) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Addw0()", Addw0);

    const auto w0_arr = w0.const_array();

    for (int lev = 0; lev <= finest_level; ++lev) {
        // need one cell-centered MF for the MFIter
        MultiFab& sold_mf = sold[lev];

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(sold_mf); mfi.isValid(); ++mfi) {
            const Array4<Real> vedge = u_edge[lev][1].array(mfi);
#if (AMREX_SPACEDIM == 3)
            const Array4<Real> uedge = u_edge[lev][0].array(mfi);
            const Array4<Real> wedge = u_edge[lev][2].array(mfi);
#endif

            if (!spherical) {
#if (AMREX_SPACEDIM == 2)
                const Box& ybx =
                    amrex::grow(mfi.nodaltilebox(1), amrex::IntVect(1, 0));

                ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    vedge(i, j, k) += mult * w0_arr(lev, j);
                });
#else
                const Box& zbx =
                    amrex::grow(mfi.nodaltilebox(2), amrex::IntVect(1, 1, 0));

                ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    wedge(i, j, k) += mult * w0_arr(lev, k);
                });
#endif
            } else {
#if (AMREX_SPACEDIM == 3)

                const Box& xbx = mfi.nodaltilebox(0);
                const Box& ybx = mfi.nodaltilebox(1);
                const Box& zbx = mfi.nodaltilebox(2);

                const Array4<const Real> w0macx = w0mac[lev][0].array(mfi);
                const Array4<const Real> w0macy = w0mac[lev][1].array(mfi);
                const Array4<const Real> w0macz = w0mac[lev][2].array(mfi);

                ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    uedge(i, j, k) += mult * w0macx(i, j, k);
                });

                ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    vedge(i, j, k) += mult * w0macy(i, j, k);
                });

                ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    wedge(i, j, k) += mult * w0macz(i, j, k);
                });
#else
                Abort("Addw0: Spherical is not valid for DIM < 3");
#endif
            }  //end spherical
        }
    }

    if (finest_level == 0) {
#ifdef AMREX_USE_CUDA
        // FillBoundary can be non-deterministic on the GPU for non-cell
        // centered data (like u_edge here). If this flag is true, run
        // on the CPU instead
        bool launched;
        if (deterministic_nodal_solve) {
            launched = !Gpu::notInLaunchRegion();
            // turn off GPU
            if (launched) Gpu::setLaunchRegion(false);
        }
#endif
        // fill periodic ghost cells
        for (int lev = 0; lev <= finest_level; ++lev) {
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                u_edge[lev][d].FillBoundary(geom[lev].periodicity());
            }
        }

#ifdef AMREX_USE_CUDA
        if (deterministic_nodal_solve) {
            // turn GPU back on
            if (launched) Gpu::setLaunchRegion(true);
        }
#endif

        // fill ghost cells behind physical boundaries
        FillUmacGhost(u_edge);

    } else {
        // edge_restriction
        AverageDownFaces(u_edge);

        // fill all ghost cells for edge-based velocity field
        FillPatchUedge(u_edge);
    }
}

#if (AMREX_SPACEDIM == 3)
void Maestro::MakeW0mac(Vector<std::array<MultiFab, AMREX_SPACEDIM> >& w0mac) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeW0mac()", MakeW0mac);

    if (!spherical) {
        Abort("Error: only call MakeW0mac for spherical");
    }

    if (w0mac[0][0].nGrow() != 1) {
        Abort("Error: MakeW0mac assumes one ghost cell");
    }

    // make nodal w0
    Vector<MultiFab> w0_nodal(finest_level + 1);

    for (int lev = 0; lev <= finest_level; ++lev) {
        w0_nodal[lev].define(convert(grids[lev], nodal_flag), dmap[lev],
                             AMREX_SPACEDIM, 1);
        w0_nodal[lev].setVal(0.);
    }

    const int nr_fine = base_geom.nr_fine;
    const int w0mac_interp_type_loc = w0mac_interp_type;
    const Real drf = base_geom.dr_fine;
    const auto w0_arr = w0.array();
    const auto r_edge_loc = base_geom.r_edge_loc;
    const auto& center_p = center;

    for (int lev = 0; lev <= finest_level; ++lev) {
        const auto dx = geom[lev].CellSizeArray();
        const auto prob_lo = geom[lev].ProbLoArray();

        // get references to the MultiFabs at level lev
        MultiFab& w0cart_mf = w0_cart[lev];

        if (w0mac_interp_type == 4) {
            // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(w0cart_mf, TilingIfNotGPU()); mfi.isValid();
                 ++mfi) {
                // Get the index space of the valid region
                const Box& gntbx = mfi.grownnodaltilebox(-1, 1);

                const Array4<Real> w0_nodal_arr = w0_nodal[lev].array(mfi);

                ParallelFor(gntbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    Real x = prob_lo[0] + Real(i) * dx[0] - center_p[0];
                    Real y = prob_lo[1] + Real(j) * dx[1] - center_p[1];
                    Real z = prob_lo[2] + Real(k) * dx[2] - center_p[2];

                    Real radius = sqrt(x * x + y * y + z * z);
                    auto index = int(radius / drf);
                    Real rfac = (radius - Real(index) * drf) / drf;

                    Real w0_cart_val;
                    if (index < nr_fine) {
                        w0_cart_val = rfac * w0_arr(0, index + 1) +
                                      (1.0 - rfac) * w0_arr(0, index);
                    } else {
                        w0_cart_val = w0_arr(0, nr_fine);
                    }

                    if (radius == 0.0) {
                        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                            w0_nodal_arr(i, j, k, n) = w0_cart_val;
                        }
                    } else {
                        w0_nodal_arr(i, j, k, 0) = w0_cart_val * x / radius;
                        w0_nodal_arr(i, j, k, 1) = w0_cart_val * y / radius;
                        w0_nodal_arr(i, j, k, 2) = w0_cart_val * z / radius;
                    }
                });
            }
        }

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(w0cart_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& xbx = mfi.grownnodaltilebox(0, 1);
            const Box& ybx = mfi.grownnodaltilebox(1, 1);
            const Box& zbx = mfi.grownnodaltilebox(2, 1);

            const Array4<const Real> w0_nodal_arr = w0_nodal[lev].array(mfi);
            const Array4<Real> w0macx = w0mac[lev][0].array(mfi);
            const Array4<Real> w0macy = w0mac[lev][1].array(mfi);
            const Array4<Real> w0macz = w0mac[lev][2].array(mfi);
            const Array4<const Real> w0_cart_arr = w0_cart[lev].array(mfi);

            if (w0mac_interp_type == 1) {
                ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    w0macx(i, j, k) = 0.5 * (w0_cart_arr(i - 1, j, k, 0) +
                                             w0_cart_arr(i, j, k, 0));
                });

                ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    w0macy(i, j, k) = 0.5 * (w0_cart_arr(i, j - 1, k, 1) +
                                             w0_cart_arr(i, j, k, 1));
                });

                ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    w0macz(i, j, k) = 0.5 * (w0_cart_arr(i, j, k - 1, 2) +
                                             w0_cart_arr(i, j, k, 2));
                });

            } else if (w0mac_interp_type == 2 || w0mac_interp_type == 3) {
                ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    Real x = prob_lo[0] + Real(i) * dx[0] - center_p[0];
                    Real y = prob_lo[1] + (Real(j) + 0.5) * dx[1] - center_p[1];
                    Real z = prob_lo[2] + (Real(k) + 0.5) * dx[2] - center_p[2];

                    Real radius = sqrt(x * x + y * y + z * z);
                    auto index = int(radius / drf);
                    Real w0_cart_val;

                    if (w0mac_interp_type_loc == 2) {
                        Real rfac = (radius - Real(index) * drf) / drf;

                        if (index < nr_fine) {
                            w0_cart_val = rfac * w0_arr(0, index + 1) +
                                          (1.0 - rfac) * w0_arr(0, index);
                        } else {
                            w0_cart_val = w0_arr(0, nr_fine);
                        }
                    } else {
                        if (index <= 0) {
                            index = 0;
                        } else if (index >= nr_fine - 1) {
                            index = nr_fine - 2;
                        } else if (radius - r_edge_loc(0, index) <
                                   r_edge_loc(0, index + 1)) {
                            index--;
                        }

                        w0_cart_val = QuadInterp(
                            radius, r_edge_loc(0, index),
                            r_edge_loc(0, index + 1), r_edge_loc(0, index + 2),
                            w0_arr(0, index), w0_arr(0, index + 1),
                            w0_arr(0, index + 2));
                    }

                    w0macx(i, j, k) = w0_cart_val * x / radius;
                });

                ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    Real x = prob_lo[0] + (Real(i) + 0.5) * dx[0] - center_p[0];
                    Real y = prob_lo[1] + Real(j) * dx[1] - center_p[1];
                    Real z = prob_lo[2] + (Real(k) + 0.5) * dx[2] - center_p[2];

                    Real radius = sqrt(x * x + y * y + z * z);
                    auto index = int(radius / drf);
                    Real w0_cart_val;

                    if (w0mac_interp_type_loc == 2) {
                        Real rfac = (radius - Real(index) * drf) / drf;

                        if (index < nr_fine) {
                            w0_cart_val = rfac * w0_arr(0, index + 1) +
                                          (1.0 - rfac) * w0_arr(0, index);
                        } else {
                            w0_cart_val = w0_arr(0, nr_fine);
                        }
                    } else {
                        if (index <= 0) {
                            index = 0;
                        } else if (index >= nr_fine - 1) {
                            index = nr_fine - 2;
                        } else if (radius - r_edge_loc(0, index) <
                                   r_edge_loc(0, index + 1)) {
                            index--;
                        }

                        w0_cart_val = QuadInterp(
                            radius, r_edge_loc(0, index),
                            r_edge_loc(0, index + 1), r_edge_loc(0, index + 2),
                            w0_arr(0, index), w0_arr(0, index + 1),
                            w0_arr(0, index + 2));
                    }

                    w0macy(i, j, k) = w0_cart_val * y / radius;
                });

                ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    Real x = prob_lo[0] + (Real(i) + 0.5) * dx[0] - center_p[0];
                    Real y = prob_lo[1] + (Real(j) + 0.5) * dx[1] - center_p[1];
                    Real z = prob_lo[2] + Real(k) * dx[2] - center_p[2];

                    Real radius = sqrt(x * x + y * y + z * z);
                    auto index = int(radius / drf);
                    Real w0_cart_val;

                    if (w0mac_interp_type_loc == 2) {
                        Real rfac = (radius - Real(index) * drf) / drf;

                        if (index < nr_fine) {
                            w0_cart_val = rfac * w0_arr(0, index + 1) +
                                          (1.0 - rfac) * w0_arr(0, index);
                        } else {
                            w0_cart_val = w0_arr(0, nr_fine);
                        }
                    } else {
                        if (index <= 0) {
                            index = 0;
                        } else if (index >= nr_fine - 1) {
                            index = nr_fine - 2;
                        } else if (radius - r_edge_loc(0, index) <
                                   r_edge_loc(0, index + 1)) {
                            index--;
                        }

                        w0_cart_val = QuadInterp(
                            radius, r_edge_loc(0, index),
                            r_edge_loc(0, index + 1), r_edge_loc(0, index + 2),
                            w0_arr(0, index), w0_arr(0, index + 1),
                            w0_arr(0, index + 2));
                    }

                    w0macz(i, j, k) = w0_cart_val * z / radius;
                });

            } else if (w0mac_interp_type == 4) {
                ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    w0macx(i, j, k) = 0.25 * (w0_nodal_arr(i, j, k, 0) +
                                              w0_nodal_arr(i, j + 1, k, 0) +
                                              w0_nodal_arr(i, j, k + 1, 0) +
                                              w0_nodal_arr(i, j + 1, k + 1, 0));
                });

                ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    w0macy(i, j, k) = 0.25 * (w0_nodal_arr(i, j, k, 1) +
                                              w0_nodal_arr(i + 1, j, k, 1) +
                                              w0_nodal_arr(i, j, k + 1, 1) +
                                              w0_nodal_arr(i + 1, j, k + 1, 1));
                });

                ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    w0macz(i, j, k) = 0.25 * (w0_nodal_arr(i, j, k, 2) +
                                              w0_nodal_arr(i + 1, j, k, 2) +
                                              w0_nodal_arr(i, j + 1, k, 2) +
                                              w0_nodal_arr(i + 1, j + 1, k, 2));
                });
            }
        }
    }
}

void Maestro::MakeS0mac(const BaseState<Real>& s0,
                        Vector<std::array<MultiFab, AMREX_SPACEDIM> >& s0mac) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeS0mac()", MakeS0mac);

    if (!spherical) {
        Abort("Error: only call MakeS0mac for spherical");
    }

    // Construct a cartesian version of s0
    Vector<MultiFab> s0_cart(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        s0_cart[lev].define(grids[lev], dmap[lev], 1, 2);
        s0_cart[lev].setVal(0.);
    }

    if (s0mac_interp_type == 1) {
        Put1dArrayOnCart(s0, s0_cart, false, false, bcs_f, 0);
    }

    if (s0mac[0][0].nGrow() != 1) {
        Abort("Error: MakeS0mac assumes one ghost cell");
    }

    const int nr_fine = base_geom.nr_fine;
    const Real drf = base_geom.dr_fine;
    const auto& r_cc_loc = base_geom.r_cc_loc;
    const auto center_p = center;
    const auto s0_arr = s0.const_array();

    for (int lev = 0; lev <= finest_level; ++lev) {
        const auto dx = geom[lev].CellSizeArray();
        const auto prob_lo = geom[lev].ProbLoArray();

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(s0_cart[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region
            const Box& xbx = mfi.grownnodaltilebox(0, 1);
            const Box& ybx = mfi.grownnodaltilebox(1, 1);
            const Box& zbx = mfi.grownnodaltilebox(2, 1);

            const Array4<Real> s0macx = s0mac[lev][0].array(mfi);
            const Array4<Real> s0macy = s0mac[lev][1].array(mfi);
            const Array4<Real> s0macz = s0mac[lev][2].array(mfi);
            const Array4<const Real> s0_cart_arr = s0_cart[lev].array(mfi);

            if (use_exact_base_state) {
                // we currently have three different ideas for computing s0mac
                // 1.  Interpolate s0 to cell centers, then average to edges
                // 2.  Interpolate s0 to edges directly using linear interpolation
                // 3.  Interpolate s0 to edges directly using quadratic interpolation
                // 4.  Interpolate s0 to nodes, then average to edges

                if (s0mac_interp_type == 1) {
                    ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        s0macx(i, j, k) = 0.5 * (s0_cart_arr(i - 1, j, k) +
                                                 s0_cart_arr(i, j, k));
                    });

                    ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        s0macy(i, j, k) = 0.5 * (s0_cart_arr(i, j - 1, k) +
                                                 s0_cart_arr(i, j, k));
                    });

                    ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        s0macz(i, j, k) = 0.5 * (s0_cart_arr(i, j, k - 1) +
                                                 s0_cart_arr(i, j, k));
                    });

                } else if (s0mac_interp_type == 2) {
                    ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        Real x = prob_lo[0] + Real(i) * dx[0] - center_p[0];
                        Real y =
                            prob_lo[1] + (Real(j) + 0.5) * dx[1] - center_p[1];
                        Real z =
                            prob_lo[2] + (Real(k) + 0.5) * dx[2] - center_p[2];

                        Real radius = sqrt(x * x + y * y + z * z);
                        auto index = (int)amrex::Math::round(
                            radius * radius / (dx[0] * dx[0]) - 0.375);
                        // closest radial index to edge-centered point

                        if (radius >= r_cc_loc(0, index)) {
                            Real dri =
                                r_cc_loc(0, index + 1) - r_cc_loc(0, index);
                            if (index >= nr_fine - 1) {
                                s0macx(i, j, k) = s0_arr(0, nr_fine - 1);
                            } else {
                                s0macx(i, j, k) =
                                    s0_arr(0, index + 1) *
                                        (radius - r_cc_loc(0, index)) / dri +
                                    s0_arr(0, index) *
                                        (r_cc_loc(0, index + 1) - radius) / dri;
                            }
                        } else {
                            Real dri =
                                r_cc_loc(0, index) - r_cc_loc(0, index - 1);
                            if (index == 0) {
                                s0macx(i, j, k) = s0_arr(0, index);
                            } else if (index > nr_fine - 1) {
                                s0macx(i, j, k) = s0_arr(0, nr_fine - 1);
                            } else {
                                s0macx(i, j, k) =
                                    s0_arr(0, index) *
                                        (radius - r_cc_loc(0, index - 1)) /
                                        dri +
                                    s0_arr(0, index - 1) *
                                        (r_cc_loc(0, index) - radius) / dri;
                            }
                        }
                    });

                    ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        Real x =
                            prob_lo[0] + (Real(i) + 0.5) * dx[0] - center_p[0];
                        Real y = prob_lo[1] + Real(j) * dx[1] - center_p[1];
                        Real z =
                            prob_lo[2] + (Real(k) + 0.5) * dx[2] - center_p[2];

                        Real radius = sqrt(x * x + y * y + z * z);
                        auto index = (int)amrex::Math::round(
                            radius * radius / (dx[1] * dx[1]) - 0.375);
                        // closest radial index to edge-centered point

                        if (radius >= r_cc_loc(0, index)) {
                            Real dri =
                                r_cc_loc(0, index + 1) - r_cc_loc(0, index);
                            if (index >= nr_fine - 1) {
                                s0macy(i, j, k) = s0_arr(0, nr_fine - 1);
                            } else {
                                s0macy(i, j, k) =
                                    s0_arr(0, index + 1) *
                                        (radius - r_cc_loc(0, index)) / dri +
                                    s0_arr(0, index) *
                                        (r_cc_loc(0, index + 1) - radius) / dri;
                            }
                        } else {
                            Real dri =
                                r_cc_loc(0, index) - r_cc_loc(0, index - 1);
                            if (index == 0) {
                                s0macy(i, j, k) = s0_arr(0, index);
                            } else if (index > nr_fine - 1) {
                                s0macy(i, j, k) = s0_arr(0, nr_fine - 1);
                            } else {
                                s0macy(i, j, k) =
                                    s0_arr(0, index) *
                                        (radius - r_cc_loc(0, index - 1)) /
                                        dri +
                                    s0_arr(0, index - 1) *
                                        (r_cc_loc(0, index) - radius) / dri;
                            }
                        }
                    });

                    ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        Real x =
                            prob_lo[0] + (Real(i) + 0.5) * dx[0] - center_p[0];
                        Real y =
                            prob_lo[1] + (Real(j) + 0.5) * dx[1] - center_p[1];
                        Real z = prob_lo[2] + Real(k) * dx[2] - center_p[2];

                        Real radius = sqrt(x * x + y * y + z * z);
                        auto index = (int)amrex::Math::round(
                            radius * radius / (dx[2] * dx[2]) - 0.375);
                        // closest radial index to edge-centered point

                        if (radius >= r_cc_loc(0, index)) {
                            Real dri =
                                r_cc_loc(0, index + 1) - r_cc_loc(0, index);
                            if (index >= nr_fine - 1) {
                                s0macz(i, j, k) = s0_arr(0, nr_fine - 1);
                            } else {
                                s0macz(i, j, k) =
                                    s0_arr(0, index + 1) *
                                        (radius - r_cc_loc(0, index)) / dri +
                                    s0_arr(0, index) *
                                        (r_cc_loc(0, index + 1) - radius) / dri;
                            }
                        } else {
                            Real dri =
                                r_cc_loc(0, index) - r_cc_loc(0, index - 1);
                            if (index == 0) {
                                s0macz(i, j, k) = s0_arr(0, index);
                            } else if (index > nr_fine - 1) {
                                s0macz(i, j, k) = s0_arr(0, nr_fine - 1);
                            } else {
                                s0macz(i, j, k) =
                                    s0_arr(0, index) *
                                        (radius - r_cc_loc(0, index - 1)) /
                                        dri +
                                    s0_arr(0, index - 1) *
                                        (r_cc_loc(0, index) - radius) / dri;
                            }
                        }
                    });

                } else if (s0mac_interp_type == 3) {
                    ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        Real x = prob_lo[0] + Real(i) * dx[0] - center_p[0];
                        Real y =
                            prob_lo[1] + (Real(j) + 0.5) * dx[1] - center_p[1];
                        Real z =
                            prob_lo[2] + (Real(k) + 0.5) * dx[2] - center_p[2];

                        Real radius = sqrt(x * x + y * y + z * z);
                        auto index = (int)amrex::Math::round(
                            radius * radius / (dx[0] * dx[0]) - 0.375);
                        // closest radial index to edge-centered point

                        // index refers to the center point in the quadratic stencil.
                        // we need to modify this if we're too close to the edge
                        if (index == 0) {
                            index = 1;
                        } else if (index >= nr_fine - 1) {
                            index = nr_fine - 2;
                        }

                        s0macx(i, j, k) = QuadInterp(
                            radius, r_cc_loc(0, index - 1), r_cc_loc(0, index),
                            r_cc_loc(0, index + 1), s0_arr(0, index - 1),
                            s0_arr(0, index), s0_arr(0, index + 1));
                    });

                    ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        Real x =
                            prob_lo[0] + (Real(i) + 0.5) * dx[0] - center_p[0];
                        Real y = prob_lo[1] + Real(j) * dx[1] - center_p[1];
                        Real z =
                            prob_lo[2] + (Real(k) + 0.5) * dx[2] - center_p[2];

                        Real radius = sqrt(x * x + y * y + z * z);
                        auto index = (int)amrex::Math::round(
                            radius * radius / (dx[1] * dx[1]) - 0.375);
                        // closest radial index to edge-centered point

                        // index refers to the center point in the quadratic stencil.
                        // we need to modify this if we're too close to the edge
                        if (index == 0) {
                            index = 1;
                        } else if (index >= nr_fine - 1) {
                            index = nr_fine - 2;
                        }

                        s0macy(i, j, k) = QuadInterp(
                            radius, r_cc_loc(0, index - 1), r_cc_loc(0, index),
                            r_cc_loc(0, index + 1), s0_arr(0, index - 1),
                            s0_arr(0, index), s0_arr(0, index + 1));
                    });

                    ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        Real x =
                            prob_lo[0] + (Real(i) + 0.5) * dx[0] - center_p[0];
                        Real y =
                            prob_lo[1] + (Real(j) + 0.5) * dx[1] - center_p[1];
                        Real z = prob_lo[2] + Real(k) * dx[2] - center_p[2];

                        Real radius = sqrt(x * x + y * y + z * z);
                        auto index = (int)amrex::Math::round(
                            radius * radius / (dx[2] * dx[2]) - 0.375);
                        // closest radial index to edge-centered point

                        // index refers to the center point in the quadratic stencil.
                        // we need to modify this if we're too close to the edge
                        if (index == 0) {
                            index = 1;
                        } else if (index >= nr_fine - 1) {
                            index = nr_fine - 2;
                        }

                        s0macz(i, j, k) = QuadInterp(
                            radius, r_cc_loc(0, index - 1), r_cc_loc(0, index),
                            r_cc_loc(0, index + 1), s0_arr(0, index - 1),
                            s0_arr(0, index), s0_arr(0, index + 1));
                    });
                }

            } else {  // use_exact_base_state = 0

                if (s0mac_interp_type == 1) {
                    ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        s0macx(i, j, k) = 0.5 * (s0_cart_arr(i - 1, j, k) +
                                                 s0_cart_arr(i, j, k));
                    });

                    ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        s0macy(i, j, k) = 0.5 * (s0_cart_arr(i, j - 1, k) +
                                                 s0_cart_arr(i, j, k));
                    });

                    ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        s0macz(i, j, k) = 0.5 * (s0_cart_arr(i, j, k - 1) +
                                                 s0_cart_arr(i, j, k));
                    });

                } else if (s0mac_interp_type == 2) {
                    ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        Real x = prob_lo[0] + Real(i) * dx[0] - center_p[0];
                        Real y =
                            prob_lo[1] + (Real(j) + 0.5) * dx[1] - center_p[1];
                        Real z =
                            prob_lo[2] + (Real(k) + 0.5) * dx[2] - center_p[2];

                        Real radius = sqrt(x * x + y * y + z * z);
                        auto index = int(radius / drf);

                        if (radius >= r_cc_loc(0, index)) {
                            if (index >= nr_fine - 1) {
                                s0macx(i, j, k) = s0_arr(0, nr_fine - 1);
                            } else {
                                s0macx(i, j, k) =
                                    s0_arr(0, index + 1) *
                                        (radius - r_cc_loc(0, index)) / drf +
                                    s0_arr(0, index) *
                                        (r_cc_loc(0, index + 1) - radius) / drf;
                            }
                        } else {
                            if (index == 0) {
                                s0macx(i, j, k) = s0_arr(0, index);
                            } else if (index > nr_fine - 1) {
                                s0macx(i, j, k) = s0_arr(0, nr_fine - 1);
                            } else {
                                s0macx(i, j, k) =
                                    s0_arr(0, index) *
                                        (radius - r_cc_loc(0, index - 1)) /
                                        drf +
                                    s0_arr(0, index - 1) *
                                        (r_cc_loc(0, index) - radius) / drf;
                            }
                        }
                    });

                    ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        Real x =
                            prob_lo[0] + (Real(i) + 0.5) * dx[0] - center_p[0];
                        Real y = prob_lo[1] + Real(j) * dx[1] - center_p[1];
                        Real z =
                            prob_lo[2] + (Real(k) + 0.5) * dx[2] - center_p[2];

                        Real radius = sqrt(x * x + y * y + z * z);
                        auto index = int(radius / drf);

                        if (radius >= r_cc_loc(0, index)) {
                            if (index >= nr_fine - 1) {
                                s0macy(i, j, k) = s0_arr(0, nr_fine - 1);
                            } else {
                                s0macy(i, j, k) =
                                    s0_arr(0, index + 1) *
                                        (radius - r_cc_loc(0, index)) / drf +
                                    s0_arr(0, index) *
                                        (r_cc_loc(0, index + 1) - radius) / drf;
                            }
                        } else {
                            if (index == 0) {
                                s0macy(i, j, k) = s0_arr(0, index);
                            } else if (index > nr_fine - 1) {
                                s0macy(i, j, k) = s0_arr(0, nr_fine - 1);
                            } else {
                                s0macy(i, j, k) =
                                    s0_arr(0, index) *
                                        (radius - r_cc_loc(0, index - 1)) /
                                        drf +
                                    s0_arr(0, index - 1) *
                                        (r_cc_loc(0, index) - radius) / drf;
                            }
                        }
                    });

                    ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        Real x =
                            prob_lo[0] + (Real(i) + 0.5) * dx[0] - center_p[0];
                        Real y =
                            prob_lo[1] + (Real(j) + 0.5) * dx[1] - center_p[1];
                        Real z = prob_lo[2] + Real(k) * dx[2] - center_p[2];

                        Real radius = sqrt(x * x + y * y + z * z);
                        auto index = int(radius / drf);

                        if (radius >= r_cc_loc(0, index)) {
                            if (index >= nr_fine - 1) {
                                s0macz(i, j, k) = s0_arr(0, nr_fine - 1);
                            } else {
                                s0macz(i, j, k) =
                                    s0_arr(0, index + 1) *
                                        (radius - r_cc_loc(0, index)) / drf +
                                    s0_arr(0, index) *
                                        (r_cc_loc(0, index + 1) - radius) / drf;
                            }
                        } else {
                            if (index == 0) {
                                s0macz(i, j, k) = s0_arr(0, index);
                            } else if (index > nr_fine - 1) {
                                s0macz(i, j, k) = s0_arr(0, nr_fine - 1);
                            } else {
                                s0macz(i, j, k) =
                                    s0_arr(0, index) *
                                        (radius - r_cc_loc(0, index - 1)) /
                                        drf +
                                    s0_arr(0, index - 1) *
                                        (r_cc_loc(0, index) - radius) / drf;
                            }
                        }
                    });

                } else if (s0mac_interp_type == 3) {
                    ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        Real x = prob_lo[0] + Real(i) * dx[0] - center_p[0];
                        Real y =
                            prob_lo[1] + (Real(j) + 0.5) * dx[1] - center_p[1];
                        Real z =
                            prob_lo[2] + (Real(k) + 0.5) * dx[2] - center_p[2];

                        Real radius = sqrt(x * x + y * y + z * z);
                        auto index = int(radius / drf);

                        if (index == 0) {
                            index = 1;
                        } else if (index >= nr_fine - 1) {
                            index = nr_fine - 2;
                        }

                        s0macx(i, j, k) = QuadInterp(
                            radius, r_cc_loc(0, index - 1), r_cc_loc(0, index),
                            r_cc_loc(0, index + 1), s0_arr(0, index - 1),
                            s0_arr(0, index), s0_arr(0, index + 1));
                    });

                    ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        Real x =
                            prob_lo[0] + (Real(i) + 0.5) * dx[0] - center_p[0];
                        Real y = prob_lo[1] + Real(j) * dx[1] - center_p[1];
                        Real z =
                            prob_lo[2] + (Real(k) + 0.5) * dx[2] - center_p[2];

                        Real radius = sqrt(x * x + y * y + z * z);
                        auto index = int(radius / drf);

                        if (index == 0) {
                            index = 1;
                        } else if (index >= nr_fine - 1) {
                            index = nr_fine - 2;
                        }

                        s0macy(i, j, k) = QuadInterp(
                            radius, r_cc_loc(0, index - 1), r_cc_loc(0, index),
                            r_cc_loc(0, index + 1), s0_arr(0, index - 1),
                            s0_arr(0, index), s0_arr(0, index + 1));
                    });

                    ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        Real x =
                            prob_lo[0] + (Real(i) + 0.5) * dx[0] - center_p[0];
                        Real y =
                            prob_lo[1] + (Real(j) + 0.5) * dx[1] - center_p[1];
                        Real z = prob_lo[2] + Real(k) * dx[2] - center_p[2];

                        Real radius = sqrt(x * x + y * y + z * z);
                        auto index = int(radius / drf);

                        if (index == 0) {
                            index = 1;
                        } else if (index >= nr_fine - 1) {
                            index = nr_fine - 2;
                        }

                        s0macz(i, j, k) = QuadInterp(
                            radius, r_cc_loc(0, index - 1), r_cc_loc(0, index),
                            r_cc_loc(0, index + 1), s0_arr(0, index - 1),
                            s0_arr(0, index), s0_arr(0, index + 1));
                    });
                }
            }
        }
    }
}
#endif
void Maestro::MakeNormal() {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeNormal()", MakeNormal);

    // normal is the unit vector in the radial direction (e_r) in spherical
    // coordinates.
    //
    // in terms of Cartesian coordinates, with unit vectors e_x, e_y, e_z,
    //    e_r = sin(theta)cos(phi) e_x + sin(theta)sin(phi) e_y + cos(theta) e_z
    // or
    //    e_r = (x/R) e_x + (y/R) e_y + (z/R) e_z

    if (spherical) {
        const auto& center_p = center;

        for (int lev = 0; lev <= finest_level; ++lev) {
            const auto dx = geom[lev].CellSizeArray();
            const auto prob_lo = geom[lev].ProbLoArray();

            // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(normal[lev], TilingIfNotGPU()); mfi.isValid();
                 ++mfi) {
                const Box& tileBox = mfi.tilebox();
                const Array4<Real> normal_arr = normal[lev].array(mfi);

                ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    Real x = prob_lo[0] + (Real(i) + 0.5) * dx[0] - center_p[0];
                    Real y = prob_lo[1] + (Real(j) + 0.5) * dx[1] - center_p[1];
                    Real z = prob_lo[2] + (Real(k) + 0.5) * dx[2] - center_p[2];

                    Real inv_radius = 1.0 / sqrt(x * x + y * y + z * z);

                    normal_arr(i, j, k, 0) = x * inv_radius;
                    normal_arr(i, j, k, 1) = y * inv_radius;
                    normal_arr(i, j, k, 2) = z * inv_radius;
                });
            }
        }
    }
}

void Maestro::PutDataOnFaces(
    const Vector<MultiFab>& s_cc,
    Vector<std::array<MultiFab, AMREX_SPACEDIM> >& face,
    const bool harmonic_avg) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::PutDataOnFaces()", PutDataOnFaces);

    for (int lev = 0; lev <= finest_level; ++lev) {
        // need one cell-centered MF for the MFIter
        const MultiFab& scc_mf = s_cc[lev];

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(scc_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region
            const Box& xbx = mfi.nodaltilebox(0);
            const Box& ybx = mfi.nodaltilebox(1);
#if (AMREX_SPACEDIM == 3)
            const Box& zbx = mfi.nodaltilebox(2);
#endif

            const Array4<const Real> scc = s_cc[lev].array(mfi);
            const Array4<Real> facex = face[lev][0].array(mfi);
            const Array4<Real> facey = face[lev][1].array(mfi);
#if (AMREX_SPACEDIM == 3)
            const Array4<Real> facez = face[lev][2].array(mfi);
#endif

            if (harmonic_avg) {
                ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    Real denom = scc(i, j, k) + scc(i - 1, j, k);
                    Real prod = scc(i, j, k) * scc(i - 1, j, k);

                    if (denom != 0.0) {
                        facex(i, j, k) = 2.0 * prod / denom;
                    } else {
                        facex(i, j, k) = 0.5 * denom;
                    }
                });

                ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    Real denom = scc(i, j, k) + scc(i, j - 1, k);
                    Real prod = scc(i, j, k) * scc(i, j - 1, k);

                    if (denom != 0.0) {
                        facey(i, j, k) = 2.0 * prod / denom;
                    } else {
                        facey(i, j, k) = 0.5 * denom;
                    }
                });
#if (AMREX_SPACEDIM == 3)
                ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    Real denom = scc(i, j, k) + scc(i, j, k - 1);
                    Real prod = scc(i, j, k) * scc(i, j, k - 1);

                    if (denom != 0.0) {
                        facez(i, j, k) = 2.0 * prod / denom;
                    } else {
                        facez(i, j, k) = 0.5 * denom;
                    }
                });
#endif
            } else {
                ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    facex(i, j, k) = 0.5 * (scc(i, j, k) + scc(i - 1, j, k));
                });
                ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    facey(i, j, k) = 0.5 * (scc(i, j, k) + scc(i, j - 1, k));
                });
#if (AMREX_SPACEDIM == 3)
                ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    facez(i, j, k) = 0.5 * (scc(i, j, k) + scc(i, j, k - 1));
                });
#endif
            }
        }
    }

    // Make sure that the fine edges average down onto the coarse edges (edge_restriction)
    AverageDownFaces(face);
}

#if (AMREX_SPACEDIM == 3)
void Maestro::MakeCCtoRadii() {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeCCtoRadius()", MakeCCtoRadii);

    for (int lev = 0; lev <= finest_level; ++lev) {
        // get references to the MultiFabs at level lev
        const auto dx = geom[lev].CellSizeArray();
        const auto dx_fine = geom[max_level].CellSizeArray();

        iMultiFab& cc_to_r = cell_cc_to_r[lev];

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(cc_to_r, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            InitBaseStateMapSphr(lev, mfi, dx_fine, dx);
        }
    }
}
#endif

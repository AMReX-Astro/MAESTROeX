/*
   A collection of routines for mapping 1D arrays onto multi-D cartesian MultiFabs
 */

#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void
Maestro::Put1dArrayOnCart (const RealVector& s0,
                           Vector<MultiFab>& s0_cart,
                           int is_input_edge_centered,
                           int is_output_a_vector,
                           const Vector<BCRec>& bcs,
                           int sbccomp, int variable_type)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Put1dArrayOnCart()", Put1dArrayOnCart);

    int ng = s0_cart[0].nGrow();
    if (ng > 0 && bcs.size() == 0) {
        Abort("Put1dArrayOnCart with ghost cells requires bcs input");
    }

    for (int lev=0; lev<=finest_level; ++lev) {
        Put1dArrayOnCart(lev,s0,s0_cart,is_input_edge_centered,
                         is_output_a_vector,bcs,sbccomp);
    }

    int ncomp = is_output_a_vector ? AMREX_SPACEDIM : 1;

    // set covered coarse cells to be the average of overlying fine cells
    AverageDown(s0_cart,0,ncomp);

    // fill ghost cells using first-order extrapolation
    if (ng > 0) {
        FillPatch(t_old, s0_cart, s0_cart, s0_cart, 0, 0, ncomp, sbccomp, bcs,
                  variable_type);
    }
}

void
Maestro::Put1dArrayOnCart (int lev,
                           const RealVector& s0,
                           Vector<MultiFab>& s0_cart,
                           int is_input_edge_centered,
                           int is_output_a_vector,
                           const Vector<BCRec>& bcs,
                           int sbccomp)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Put1dArrayOnCart_lev()",Put1dArrayOnCart);

    // get references to the MultiFabs at level lev
    MultiFab& s0_cart_mf = s0_cart[lev];

    const auto dx = geom[lev].CellSizeArray();
    const auto prob_lo = geom[lev].ProbLoArray();
    GpuArray<Real,AMREX_SPACEDIM> center;

    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        center[n] = 0.5 * (geom[lev].ProbLo(n) + geom[lev].ProbHi(n));
    }

    Real * AMREX_RESTRICT r_edge_loc_p = r_edge_loc.dataPtr();
    Real * AMREX_RESTRICT r_cc_loc_p = r_cc_loc.dataPtr();
    const Real * AMREX_RESTRICT s0_p = s0.dataPtr();

    const int max_lev = max_radial_level+1;
    const int nr_fine_loc = nr_fine;
    const int w0_interp_type_loc = w0_interp_type;

    // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
    for ( MFIter mfi(s0_cart_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        // Get the index space of the valid region
        const Box& tileBox = mfi.tilebox();

        const Array4<Real> s0_cart_arr = s0_cart[lev].array(mfi);

        if (spherical == 0) {

            const int outcomp = is_output_a_vector == 1 ? AMREX_SPACEDIM-1 : 0;

            AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                const int r = AMREX_SPACEDIM == 2 ? j : k;

                s0_cart_arr(i,j,k,outcomp) = is_input_edge_centered == 1 ? 
                    0.5 * (s0_p[lev+r*max_lev] + s0_p[lev+(r+1)*max_lev]) : 
                    s0_p[lev+r*max_lev];
            });

        } else {

            const Array4<const Real> cc_to_r = cell_cc_to_r[lev].array(mfi);

            if (use_exact_base_state) {
                if (is_input_edge_centered) {
                    // we implemented three different ideas for computing s0_cart,
                    // where s0 is edge-centered.
                    // 1.  Piecewise constant
                    // 2.  Piecewise linear
                    // 3.  Quadratic

                    AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {

                        Real x = prob_lo[0] + (Real(i)+0.5) * dx[0] - center[0];
                        Real y = prob_lo[1] + (Real(j)+0.5) * dx[1] - center[1];
                        Real z = prob_lo[2] + (Real(k)+0.5) * dx[2] - center[2];

                        Real radius = sqrt(x*x + y*y + z*z);
                        int index = cc_to_r(i,j,k);

                        Real rfac;
                        if (index < nr_fine) {
                            rfac = (radius - r_edge_loc_p[(index+1)*max_lev]) 
                            / (r_cc_loc_p[(index+1)*max_lev] 
                                - r_cc_loc_p[index*max_lev]);
                        } else {
                            rfac = (radius - r_edge_loc_p[(index+1)*max_lev]) 
                            / (r_cc_loc_p[index*max_lev] 
                                - r_cc_loc_p[(index-1)*max_lev]);
                        }

                        Real s0_cart_val;

                        if (w0_interp_type_loc == 1) {

                            s0_cart_val = rfac > 0.5 ? 
                                s0_p[(index+1)*max_lev] : s0_p[index*max_lev];

                        } else if (w0_interp_type_loc == 2) {
                            
                            if (index < nr_fine_loc) {
                                s0_cart_val = rfac * s0_p[(index+1)*max_lev] 
                                    + (1.0-rfac) * s0_p[index*max_lev];
                            } else {
                                s0_cart_val = s0_p[nr_fine_loc*max_lev];
                            }

                        } else if (w0_interp_type_loc == 3) {
                            if (index <= 0) {
                                index = 0;
                            } else if (index >= nr_fine_loc-1) {
                                index = nr_fine_loc - 2;
                            } else if (radius-r_edge_loc_p[index*max_lev] 
                                    < r_edge_loc_p[(index+1)*max_lev]) {
                                index--;
                            }

                            s0_cart_val = QuadInterp(radius, 
                                r_edge_loc_p[index*max_lev],
                                r_edge_loc_p[(index+1)*max_lev], 
                                r_edge_loc_p[(index+2)*max_lev], 
                                s0_p[index*max_lev],
                                s0_p[(index+1)*max_lev],
                                s0_p[(index+2)*max_lev]);
                        }

                        if (is_output_a_vector) {
                            s0_cart_arr(i,j,k,0) = s0_cart_val * x / radius;
                            s0_cart_arr(i,j,k,1) = s0_cart_val * y / radius;
                            s0_cart_arr(i,j,k,2) = s0_cart_val * z / radius;
                        } else {
                            s0_cart_arr(i,j,k,0) = s0_cart_val;
                        }
                    });

                } else { // is_input_edge_centered = 0
                    // we directly inject the spherical values into each cell center
                    // because s0 is also bin-centered.

                    AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {

                        Real x = prob_lo[0] + (Real(i)+0.5) * dx[0] - center[0];
                        Real y = prob_lo[1] + (Real(j)+0.5) * dx[1] - center[1];
                        Real z = prob_lo[2] + (Real(k)+0.5) * dx[2] - center[2];

                        Real radius = sqrt(x*x + y*y + z*z);
                        int index = cc_to_r(i,j,k);

                        Real s0_cart_val = s0_p[index*max_lev];
                        
                        if (is_output_a_vector) {
                            s0_cart_arr(i,j,k,0) = s0_cart_val * x / radius;
                            s0_cart_arr(i,j,k,1) = s0_cart_val * y / radius;
                            s0_cart_arr(i,j,k,2) = s0_cart_val * z / radius;
                        } else {
                            s0_cart_arr(i,j,k,0) = s0_cart_val;
                        }
                    });
                } // is_input_edge_centered

            } else { // use_exact_base_state = 0

                const Real dr = dr_fine;

                if (is_input_edge_centered) {
                    // we implemented three different ideas for computing s0_cart,
                    // where s0 is edge-centered.
                    // 1.  Piecewise constant
                    // 2.  Piecewise linear
                    // 3.  Quadratic

                    AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {

                        Real x = prob_lo[0] + (Real(i)+0.5) * dx[0] - center[0];
                        Real y = prob_lo[1] + (Real(j)+0.5) * dx[1] - center[1];
                        Real z = prob_lo[2] + (Real(k)+0.5) * dx[2] - center[2];

                        Real radius = sqrt(x*x + y*y + z*z);
                        int index = int(radius / dr);

                        Real rfac = (radius - Real(index) * dr) / dr;
                        Real s0_cart_val = 0.0;

                        if (w0_interp_type_loc == 1) {

                            s0_cart_val = rfac > 0.5 ? 
                                s0_p[(index+1)*max_lev] : s0_p[index*max_lev];

                        } else if (w0_interp_type_loc == 2) {
                            
                            if (index < nr_fine_loc) {
                                s0_cart_val = rfac * s0_p[(index+1)*max_lev] + (1.0-rfac) * s0_p[index*max_lev];
                            } else {
                                s0_cart_val = s0_p[nr_fine_loc*max_lev];
                            }

                        } else if (w0_interp_type_loc == 3) {
                            
                            if (index <= 0) {
                                index = 0;
                            } else if (index >= nr_fine_loc-1) {
                                index = nr_fine_loc - 2;
                            } else if (radius-r_edge_loc_p[index*max_lev] 
                                    < r_edge_loc_p[(index+1)*max_lev]) {
                                index--;
                            }

                            s0_cart_val = QuadInterp(radius, 
                                r_edge_loc_p[index*max_lev],
                                r_edge_loc_p[(index+1)*max_lev], 
                                r_edge_loc_p[(index+2)*max_lev], 
                                s0_p[index*max_lev],
                                s0_p[(index+1)*max_lev],
                                s0_p[(index+2)*max_lev]);
                        }
                        
                        if (is_output_a_vector) {
                            s0_cart_arr(i,j,k,0) = s0_cart_val * x / radius;
                            s0_cart_arr(i,j,k,1) = s0_cart_val * y / radius;
                            s0_cart_arr(i,j,k,2) = s0_cart_val * z / radius;
                        } else {
                            s0_cart_arr(i,j,k,0) = s0_cart_val;
                        }
                    });

                } else { // is_input_edge_centered = 0

                    const int s0_interp_type_loc = s0_interp_type;

                    // we currently have three different ideas for computing s0_cart,
                    // where s0 is bin-centered.
                    // 1.  Piecewise constant
                    // 2.  Piecewise linear
                    // 3.  Quadratic
                    AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {

                        Real x = prob_lo[0] + (Real(i)+0.5) * dx[0] - center[0];
                        Real y = prob_lo[1] + (Real(j)+0.5) * dx[1] - center[1];
                        Real z = prob_lo[2] + (Real(k)+0.5) * dx[2] - center[2];

                        Real radius = sqrt(x*x + y*y + z*z);
                        int index = int(radius / dr);

                        Real s0_cart_val = 0.0;

                        if (s0_interp_type_loc == 1) {

                            s0_cart_val = s0_p[index*max_lev];

                        } else if (s0_interp_type_loc == 2) {

                            if (radius >= r_cc_loc_p[index*max_lev]) {
                                if (index >= nr_fine_loc-1) {
                                    s0_cart_val = s0_p[(nr_fine_loc-1)*max_lev];
                                } else {
                                    s0_cart_val = s0_p[(index+1)*max_lev] 
                                        * (radius-r_cc_loc_p[index*max_lev])/dr 
                                        + s0_p[index*max_lev] 
                                        * (r_cc_loc_p[(index+1)*max_lev]-radius)/dr;
                                }
                            } else {
                                if (index == 0) {
                                    s0_cart_val = s0_p[index*max_lev];
                                } else if (index > nr_fine_loc-1) {
                                    s0_cart_val = s0_p[(nr_fine_loc-1)*max_lev];
                                } else {
                                    s0_cart_val = s0_p[index*max_lev] 
                                        * (radius-r_cc_loc_p[(index-1)*max_lev])/dr 
                                        + s0_p[(index-1)*max_lev] 
                                        * (r_cc_loc_p[index*max_lev]-radius)/dr;
                                }
                            }
                        } else if (s0_interp_type_loc == 3) {
                            if (index == 0) {
                                index = 1;
                            } else if (index >= nr_fine_loc-1) {
                                index = nr_fine_loc-2;
                            }

                            s0_cart_val = QuadInterp(radius, 
                                r_cc_loc_p[(index-1)*max_lev],
                                r_cc_loc_p[index*max_lev],
                                r_cc_loc_p[(index+1)*max_lev], 
                                s0_p[(index-1)*max_lev],
                                s0_p[index*max_lev],
                                s0_p[(index+1)*max_lev]);
                        }
                        
                        if (is_output_a_vector) {
                            s0_cart_arr(i,j,k,0) = s0_cart_val * x / radius;
                            s0_cart_arr(i,j,k,1) = s0_cart_val * y / radius;
                            s0_cart_arr(i,j,k,2) = s0_cart_val * z / radius;
                        } else {
                            s0_cart_arr(i,j,k,0) = s0_cart_val;
                        }
                    });
                } // is_input_edge_centered
            } // use_exact_base_state
        }
    }
}

AMREX_GPU_DEVICE 
Real
Maestro::QuadInterp(const Real x, const Real x0, const Real x1, const Real x2,
                    const Real y0, const Real y1, const Real y2, bool limit) 
{
    Real y = y0 + (y1-y0)/(x1-x0)*(x-x0) 
         + ((y2-y1)/(x2-x1)-(y1-y0)/(x1-x0))/(x2-x0)*(x-x0)*(x-x1);

    if (limit) {
        if (y > max(y0, max(y1, y2))) y = max(y0, max(y1, y2));
        if (y < min(y0, min(y1, y2))) y = min(y0, min(y1, y2));
    }

    return y;
}


void
Maestro::Addw0 (Vector<std::array< MultiFab, AMREX_SPACEDIM > >& u_edge,
                const Vector<std::array< MultiFab,AMREX_SPACEDIM > >& w0mac,
                const Real& mult)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Addw0()",Addw0);

    const int max_lev = max_radial_level+1;
    const Real * AMREX_RESTRICT w0_p = w0.dataPtr();

    for (int lev=0; lev<=finest_level; ++lev) {

        // need one cell-centered MF for the MFIter
        MultiFab& sold_mf = sold[lev];

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(sold_mf); mfi.isValid(); ++mfi ) {
            
            const Array4<Real> vedge = u_edge[lev][1].array(mfi);
#if (AMREX_SPACEDIM == 3)
            const Array4<Real> uedge = u_edge[lev][0].array(mfi);
            const Array4<Real> wedge = u_edge[lev][2].array(mfi);
#endif
            
            if (spherical == 0) {
#if (AMREX_SPACEDIM == 2)
                const Box& ybx = amrex::grow(mfi.nodaltilebox(1), amrex::IntVect(1,0));

                AMREX_PARALLEL_FOR_3D(ybx, i, j, k, {
                    vedge(i,j,k) += mult * w0_p[lev+j*max_lev];
                });
#else
                const Box& zbx = amrex::grow(mfi.nodaltilebox(2), amrex::IntVect(1,1,0));

                AMREX_PARALLEL_FOR_3D(zbx, i, j, k, {
                    wedge(i,j,k) += mult * w0_p[lev+k*max_lev];
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

                AMREX_PARALLEL_FOR_3D(xbx, i, j, k, {
                    uedge(i,j,k) += mult * w0macx(i,j,k);
                });

                AMREX_PARALLEL_FOR_3D(ybx, i, j, k, {
                    vedge(i,j,k) += mult * w0macy(i,j,k);
                });

                AMREX_PARALLEL_FOR_3D(zbx, i, j, k, {
                    wedge(i,j,k) += mult * w0macz(i,j,k);
                });
#else
                Abort("Addw0: Spherical is not valid for DIM < 3");
#endif
            }  //end spherical
        }
    }

    if (finest_level == 0) {
        // fill periodic ghost cells
        for (int lev=0; lev<=finest_level; ++lev) {
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                u_edge[lev][d].FillBoundary(geom[lev].periodicity());
            }
        }

        // fill ghost cells behind physical boundaries
        FillUmacGhost(u_edge);
    } else {
        // edge_restriction
        AverageDownFaces(u_edge);

        // fill all ghost cells for edge-based velocity field
        FillPatchUedge(u_edge);
    }
}


void
Maestro::MakeW0mac (Vector<std::array< MultiFab,AMREX_SPACEDIM > >& w0mac)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeW0mac()",MakeW0mac);

    if (spherical == 0) {
        Abort("Error: only call MakeW0mac for spherical");
    }

    if (w0mac[0][0].nGrow() != 1) {
        Abort("Error: MakeW0mac assumes one ghost cell");
    }

    // make nodal w0
    Vector<MultiFab> w0_nodal(finest_level+1);

    for (int lev=0; lev<=finest_level; ++lev) {
        w0_nodal[lev].define(convert(grids[lev],nodal_flag), dmap[lev], AMREX_SPACEDIM, 1);
        w0_nodal[lev].setVal(0.);
    }

    const int nr_fine_loc = nr_fine;
    const int max_lev = max_radial_level+1;
    const int w0mac_interp_type_loc = w0mac_interp_type;
    const Real dr = dr_fine;
    const Real * AMREX_RESTRICT w0_p = w0.dataPtr();
    Real * AMREX_RESTRICT r_edge_loc_p = r_edge_loc.dataPtr();

    for (int lev=0; lev<=finest_level; ++lev) {
    
        const auto dx = geom[lev].CellSizeArray();
        const auto prob_lo = geom[lev].ProbLoArray();
        GpuArray<Real,AMREX_SPACEDIM> center;

        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            center[n] = 0.5 * (geom[lev].ProbLo(n) + geom[lev].ProbHi(n));
        }

        // get references to the MultiFabs at level lev
        MultiFab& w0cart_mf = w0_cart[lev];

        if (w0mac_interp_type == 4) {

            // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
            for ( MFIter mfi(w0cart_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

                // Get the index space of the valid region
                const Box& gntbx = mfi.grownnodaltilebox(-1, 1);

                const Array4<Real> w0_nodal_arr = w0_nodal[lev].array(mfi);

                AMREX_PARALLEL_FOR_3D(gntbx, i, j, k, {
                    Real x = prob_lo[0] + Real(i) * dx[0] - center[0];
                    Real y = prob_lo[1] + Real(j) * dx[1] - center[1];
                    Real z = prob_lo[2] + Real(k) * dx[2] - center[2];

                    Real radius = sqrt(x*x + y*y + z*z);
                    int index = int(radius / dr);
                    Real rfac = (radius - Real(index) * dr) / dr;

                    Real w0_cart_val;
                    if (index < nr_fine_loc) {
                        w0_cart_val = rfac * w0_p[(index+1)*max_lev] + (1.0-rfac) * w0_p[index*max_lev];
                    } else {
                        w0_cart_val = w0_p[nr_fine_loc*max_lev];
                    }

                    if (radius == 0.0) {
                        for (int n=0; n < AMREX_SPACEDIM; ++n) {
                            w0_nodal_arr(i,j,k,n) = w0_cart_val;
                        }
                    } else {
                        w0_nodal_arr(i,j,k,0) = w0_cart_val * x / radius;
                        w0_nodal_arr(i,j,k,1) = w0_cart_val * y / radius;
                        w0_nodal_arr(i,j,k,2) = w0_cart_val * z / radius;
                    }
                });
            }
        }

#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(w0cart_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            const Box& xbx = mfi.grownnodaltilebox(0, 1);
            const Box& ybx = mfi.grownnodaltilebox(1, 1);
            const Box& zbx = mfi.grownnodaltilebox(2, 1);

            const Array4<const Real> w0_nodal_arr = w0_nodal[lev].array(mfi);
            const Array4<Real> w0macx = w0mac[lev][0].array(mfi);
            const Array4<Real> w0macy = w0mac[lev][1].array(mfi);
            const Array4<Real> w0macz = w0mac[lev][2].array(mfi);
            const Array4<const Real> w0_cart_arr = w0_cart[lev].array(mfi);

            if (w0mac_interp_type == 1) {

                AMREX_PARALLEL_FOR_3D(xbx, i, j, k, {
                    w0macx(i,j,k) = 0.5 * (w0_cart_arr(i-1,j,k,0) + w0_cart_arr(i,j,k,0));
                });

                AMREX_PARALLEL_FOR_3D(ybx, i, j, k, {
                    w0macy(i,j,k) = 0.5 * (w0_cart_arr(i,j-1,k,1) + w0_cart_arr(i,j,k,1));
                });

                AMREX_PARALLEL_FOR_3D(zbx, i, j, k, {
                    w0macz(i,j,k) = 0.5 * (w0_cart_arr(i,j,k-1,2) + w0_cart_arr(i,j,k,2));
                });

            } else if (w0mac_interp_type == 2 || w0mac_interp_type == 3) {

                AMREX_PARALLEL_FOR_3D(xbx, i, j, k, {
                    Real x = prob_lo[0] + Real(i) * dx[0] - center[0];
                    Real y = prob_lo[1] + (Real(j)+0.5) * dx[1] - center[1];
                    Real z = prob_lo[2] + (Real(k)+0.5) * dx[2] - center[2];

                    Real radius = sqrt(x*x + y*y + z*z);
                    int index = int(radius / dr);
                    Real w0_cart_val;

                    if (w0mac_interp_type_loc == 2) {

                        Real rfac = (radius - Real(index)*dr) / dr;

                        if (index < nr_fine_loc) {
                            w0_cart_val = rfac * w0_p[(index+1)*max_lev] + (1.0-rfac) * w0_p[index*max_lev];
                        } else {
                            w0_cart_val = w0_p[nr_fine_loc*max_lev];
                        }
                    } else {
                        if (index <= 0) {
                            index = 0;
                        } else if (index >= nr_fine_loc-1) {
                            index = nr_fine_loc-2;
                        } else if (radius - r_edge_loc_p[index*max_lev] < r_edge_loc_p[(index+1)*max_lev]) {
                            index--;
                        }

                        w0_cart_val = QuadInterp(radius, 
                                    r_edge_loc_p[index*max_lev],
                                    r_edge_loc_p[(index+1)*max_lev], 
                                    r_edge_loc_p[(index+2)*max_lev], 
                                    w0_p[index*max_lev],
                                    w0_p[(index+1)*max_lev],
                                    w0_p[(index+2)*max_lev]);
                    }

                    w0macx(i,j,k) = w0_cart_val * x / radius;
                });

                AMREX_PARALLEL_FOR_3D(ybx, i, j, k, {
                    Real x = prob_lo[0] + (Real(i)+0.5) * dx[0] - center[0];
                    Real y = prob_lo[1] + Real(j) * dx[1] - center[1];
                    Real z = prob_lo[2] + (Real(k)+0.5) * dx[2] - center[2];

                    Real radius = sqrt(x*x + y*y + z*z);
                    int index = int(radius / dr);
                    Real w0_cart_val;

                    if (w0mac_interp_type_loc == 2) {

                        Real rfac = (radius - Real(index)*dr) / dr;

                        if (index < nr_fine_loc) {
                            w0_cart_val = rfac * w0_p[(index+1)*max_lev] + (1.0-rfac) * w0_p[index*max_lev];
                        } else {
                            w0_cart_val = w0_p[nr_fine_loc*max_lev];
                        }
                    } else { 
                        if (index <= 0) {
                            index = 0;
                        } else if (index >= nr_fine_loc-1) {
                            index = nr_fine_loc-2;
                        } else if (radius - r_edge_loc_p[index*max_lev] < r_edge_loc_p[(index+1)*max_lev]) {
                            index--;
                        }

                        w0_cart_val = QuadInterp(radius, 
                                    r_edge_loc_p[index*max_lev],
                                    r_edge_loc_p[(index+1)*max_lev], 
                                    r_edge_loc_p[(index+2)*max_lev], 
                                    w0_p[index*max_lev],
                                    w0_p[(index+1)*max_lev],
                                    w0_p[(index+2)*max_lev]);
                    }

                    w0macy(i,j,k) = w0_cart_val * y / radius;
                });

                AMREX_PARALLEL_FOR_3D(zbx, i, j, k, {
                    Real x = prob_lo[0] + (Real(i)+0.5) * dx[0] - center[0];
                    Real y = prob_lo[1] + (Real(j)+0.5) * dx[1] - center[1];
                    Real z = prob_lo[2] + Real(k) * dx[2] - center[2];

                    Real radius = sqrt(x*x + y*y + z*z);
                    int index = int(radius / dr);
                    Real w0_cart_val;

                    if (w0mac_interp_type_loc == 2) {

                        Real rfac = (radius - Real(index)*dr) / dr;
                        
                        if (index < nr_fine_loc) {
                            w0_cart_val = rfac * w0_p[(index+1)*max_lev] + (1.0-rfac) * w0_p[index*max_lev];
                        } else {
                            w0_cart_val = w0_p[nr_fine_loc*max_lev];
                        }
                    } else {

                        if (index <= 0) {
                            index = 0;
                        } else if (index >= nr_fine_loc-1) {
                            index = nr_fine_loc-2;
                        } else if (radius - r_edge_loc_p[index*max_lev] < r_edge_loc_p[(index+1)*max_lev]) {
                            index--;
                        }

                        w0_cart_val = QuadInterp(radius, 
                                    r_edge_loc_p[index*max_lev],
                                    r_edge_loc_p[(index+1)*max_lev], 
                                    r_edge_loc_p[(index+2)*max_lev], 
                                    w0_p[index*max_lev],
                                    w0_p[(index+1)*max_lev],
                                    w0_p[(index+2)*max_lev]);
                    }

                    w0macz(i,j,k) = w0_cart_val * z / radius;
                });
                
            } else if (w0mac_interp_type == 4) {

                AMREX_PARALLEL_FOR_3D(xbx, i, j, k, {
                    w0macx(i,j,k) = 0.25 * (w0_nodal_arr(i,j,k,0) 
                        + w0_nodal_arr(i,j+1,k,0) + w0_nodal_arr(i,j,k+1,0) 
                        + w0_nodal_arr(i,j+1,k+1,0));
                });

                AMREX_PARALLEL_FOR_3D(ybx, i, j, k, {
                    w0macy(i,j,k) = 0.25 * (w0_nodal_arr(i,j,k,1) 
                        + w0_nodal_arr(i+1,j,k,1) + w0_nodal_arr(i,j,k+1,1) 
                        + w0_nodal_arr(i+1,j,k+1,1));
                });

                AMREX_PARALLEL_FOR_3D(zbx, i, j, k, {
                    w0macz(i,j,k) = 0.25 * (w0_nodal_arr(i,j,k,2) 
                        + w0_nodal_arr(i+1,j,k,2) + w0_nodal_arr(i,j+1,k,2) 
                        + w0_nodal_arr(i+1,j+1,k,2));
                });
            }
        }
    }
}


void
Maestro::MakeS0mac (const RealVector& s0,
                    Vector<std::array< MultiFab,AMREX_SPACEDIM > >& s0mac)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeS0mac()",MakeS0mac);

    if (spherical == 0) {
        Abort("Error: only call MakeS0mac for spherical");
    }

    // Construct a cartesian version of s0
    Vector<MultiFab> s0_cart(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        s0_cart[lev].define(grids[lev], dmap[lev], 1, 2);
        s0_cart[lev].setVal(0.);
    }

    if (s0mac_interp_type == 1) {
        Put1dArrayOnCart(s0, s0_cart, 0, 0, bcs_f, 0);
    }

    if (s0mac[0][0].nGrow() != 1) {
        Abort("Error: MakeS0mac assumes one ghost cell");
    }

    const int nr_fine_loc = nr_fine;
    const int max_lev = max_radial_level+1;
    const Real dr = dr_fine;
    const Real * AMREX_RESTRICT s0_p = s0.dataPtr();
    Real * AMREX_RESTRICT r_cc_loc_p = r_cc_loc.dataPtr();

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& s0cart_mf = s0_cart[lev];
    
        const auto dx = geom[lev].CellSizeArray();
        const auto prob_lo = geom[lev].ProbLoArray();
        GpuArray<Real,AMREX_SPACEDIM> center;

        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            center[n] = 0.5 * (geom[lev].ProbLo(n) + geom[lev].ProbHi(n));
        }

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(s0cart_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

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

                    AMREX_PARALLEL_FOR_3D(xbx, i, j, k, {
                        s0macx(i,j,k) = 0.5 * (s0_cart_arr(i-1,j,k) + s0_cart_arr(i,j,k));
                    });

                    AMREX_PARALLEL_FOR_3D(ybx, i, j, k, {
                        s0macy(i,j,k) = 0.5 * (s0_cart_arr(i,j-1,k) + s0_cart_arr(i,j,k));
                    });

                    AMREX_PARALLEL_FOR_3D(zbx, i, j, k, {
                        s0macz(i,j,k) = 0.5 * (s0_cart_arr(i,j,k-1) + s0_cart_arr(i,j,k));
                    });

                } else if (s0mac_interp_type == 2) {

                    AMREX_PARALLEL_FOR_3D(xbx, i, j, k, {
                        Real x = prob_lo[0] + Real(i) * dx[0] - center[0];
                        Real y = prob_lo[1] + (Real(j)+0.5) * dx[1] - center[1];
                        Real z = prob_lo[2] + (Real(k)+0.5) * dx[2] - center[2];

                        Real radius = sqrt(x*x + y*y + z*z);
                        int index = round(radius*radius / (dx[0]*dx[0]) - 0.375);
                        // closest radial index to edge-centered point

                        if (radius >= r_cc_loc_p[index*max_lev]) {
                            Real dri = r_cc_loc_p[(index+1)*max_lev] 
                                - r_cc_loc_p[index*max_lev];
                            if (index >= nr_fine_loc-1) {
                                s0macx(i,j,k) = s0_p[(nr_fine-1)*max_lev];
                            } else {
                                s0macx(i,j,k) = s0_p[(index+1)*max_lev] 
                                    * (radius-r_cc_loc_p[index*max_lev])/dri
                                    + s0_p[index*max_lev]
                                    * (r_cc_loc_p[(index+1)*max_lev]-radius)/dri;
                            }
                        } else {
                            Real dri = r_cc_loc_p[index*max_lev] 
                                - r_cc_loc_p[(index-1)*max_lev];
                            if (index == 0) {
                                s0macx(i,j,k) = s0_p[index*max_lev];
                            } else if (index > nr_fine_loc-1) {
                                s0macx(i,j,k) = s0_p[(nr_fine-1)*max_lev];
                            } else {
                                s0macx(i,j,k) = s0_p[index*max_lev] 
                                    * (radius-r_cc_loc_p[(index-1)*max_lev])/dri
                                    + s0_p[(index-1)*max_lev]
                                    * (r_cc_loc_p[index*max_lev]-radius)/dri;
                            }
                        }
                    });

                    AMREX_PARALLEL_FOR_3D(ybx, i, j, k, {
                        Real x = prob_lo[0] + (Real(i)+0.5) * dx[0] - center[0];
                        Real y = prob_lo[1] + Real(j) * dx[1] - center[1];
                        Real z = prob_lo[2] + (Real(k)+0.5) * dx[2] - center[2];

                        Real radius = sqrt(x*x + y*y + z*z);
                        int index = round(radius*radius / (dx[1]*dx[1]) - 0.375);
                        // closest radial index to edge-centered point

                        if (radius >= r_cc_loc_p[index*max_lev]) {
                            Real dri = r_cc_loc_p[(index+1)*max_lev] 
                                - r_cc_loc_p[index*max_lev];
                            if (index >= nr_fine_loc-1) {
                                s0macy(i,j,k) = s0_p[(nr_fine-1)*max_lev];
                            } else {
                                s0macy(i,j,k) = s0_p[(index+1)*max_lev] 
                                    * (radius-r_cc_loc_p[index*max_lev])/dri
                                    + s0_p[index*max_lev]
                                    * (r_cc_loc_p[(index+1)*max_lev]-radius)/dri;
                            }
                        } else {
                            Real dri = r_cc_loc_p[index*max_lev] 
                                - r_cc_loc_p[(index-1)*max_lev];
                            if (index == 0) {
                                s0macy(i,j,k) = s0_p[index*max_lev];
                            } else if (index > nr_fine_loc-1) {
                                s0macy(i,j,k) = s0_p[(nr_fine-1)*max_lev];
                            } else {
                                s0macy(i,j,k) = s0_p[index*max_lev] 
                                    * (radius-r_cc_loc_p[(index-1)*max_lev])/dri
                                    + s0_p[(index-1)*max_lev]
                                    * (r_cc_loc_p[index*max_lev]-radius)/dri;
                            }
                        }
                    });

                    AMREX_PARALLEL_FOR_3D(zbx, i, j, k, {
                        Real x = prob_lo[0] + (Real(i)+0.5) * dx[0] - center[0];
                        Real y = prob_lo[1] + (Real(j)+0.5) * dx[1] - center[1];
                        Real z = prob_lo[2] + Real(k) * dx[2] - center[2];

                        Real radius = sqrt(x*x + y*y + z*z);
                        int index = round(radius*radius / (dx[2]*dx[2]) - 0.375);
                        // closest radial index to edge-centered point

                        if (radius >= r_cc_loc_p[index*max_lev]) {
                            Real dri = r_cc_loc_p[(index+1)*max_lev] 
                                - r_cc_loc_p[index*max_lev];
                            if (index >= nr_fine_loc-1) {
                                s0macz(i,j,k) = s0_p[(nr_fine-1)*max_lev];
                            } else {
                                s0macz(i,j,k) = s0_p[(index+1)*max_lev] 
                                    * (radius-r_cc_loc_p[index*max_lev])/dri
                                    + s0_p[index*max_lev]
                                    * (r_cc_loc_p[(index+1)*max_lev]-radius)/dri;
                            }
                        } else {
                            Real dri = r_cc_loc_p[index*max_lev] 
                                - r_cc_loc_p[(index-1)*max_lev];
                            if (index == 0) {
                                s0macz(i,j,k) = s0_p[index*max_lev];
                            } else if (index > nr_fine_loc-1) {
                                s0macz(i,j,k) = s0_p[(nr_fine-1)*max_lev];
                            } else {
                                s0macz(i,j,k) = s0_p[index*max_lev] 
                                    * (radius-r_cc_loc_p[(index-1)*max_lev])/dri
                                    + s0_p[(index-1)*max_lev]
                                    * (r_cc_loc_p[index*max_lev]-radius)/dri;
                            }
                        }
                    });

                } else if (s0mac_interp_type == 3) {

                    AMREX_PARALLEL_FOR_3D(xbx, i, j, k, {
                        Real x = prob_lo[0] + Real(i) * dx[0] - center[0];
                        Real y = prob_lo[1] + (Real(j)+0.5) * dx[1] - center[1];
                        Real z = prob_lo[2] + (Real(k)+0.5) * dx[2] - center[2];

                        Real radius = sqrt(x*x + y*y + z*z);
                        int index = round(radius*radius / (dx[0]*dx[0]) - 0.375);
                        // closest radial index to edge-centered point

                        // index refers to the center point in the quadratic stencil.
                        // we need to modify this if we're too close to the edge
                        if (index == 0) {
                            index = 1;
                        } else if (index >= nr_fine_loc-1) {
                            index = nr_fine_loc-2;
                        }

                        s0macx(i,j,k) = QuadInterp(radius, 
                                            r_cc_loc_p[(index-1)*max_lev],
                                            r_cc_loc_p[index*max_lev], 
                                            r_cc_loc_p[(index+1)*max_lev], 
                                            s0_p[(index-1)*max_lev],
                                            s0_p[index*max_lev],
                                            s0_p[(index+1)*max_lev]);
                    });

                    AMREX_PARALLEL_FOR_3D(ybx, i, j, k, {
                        Real x = prob_lo[0] + (Real(i)+0.5) * dx[0] - center[0];
                        Real y = prob_lo[1] + Real(j) * dx[1] - center[1];
                        Real z = prob_lo[2] + (Real(k)+0.5) * dx[2] - center[2];

                        Real radius = sqrt(x*x + y*y + z*z);
                        int index = round(radius*radius / (dx[1]*dx[1]) - 0.375);
                        // closest radial index to edge-centered point

                        // index refers to the center point in the quadratic stencil.
                        // we need to modify this if we're too close to the edge
                        if (index == 0) {
                            index = 1;
                        } else if (index >= nr_fine_loc-1) {
                            index = nr_fine_loc-2;
                        }

                        s0macy(i,j,k) = QuadInterp(radius, 
                                            r_cc_loc_p[(index-1)*max_lev],
                                            r_cc_loc_p[index*max_lev], 
                                            r_cc_loc_p[(index+1)*max_lev], 
                                            s0_p[(index-1)*max_lev],
                                            s0_p[index*max_lev],
                                            s0_p[(index+1)*max_lev]);
                    });

                    AMREX_PARALLEL_FOR_3D(zbx, i, j, k, {
                        Real x = prob_lo[0] + (Real(i)+0.5) * dx[0] - center[0];
                        Real y = prob_lo[1] + (Real(j)+0.5) * dx[1] - center[1];
                        Real z = prob_lo[2] + Real(k) * dx[2] - center[2];

                        Real radius = sqrt(x*x + y*y + z*z);
                        int index = round(radius*radius / (dx[2]*dx[2]) - 0.375);
                        // closest radial index to edge-centered point

                        // index refers to the center point in the quadratic stencil.
                        // we need to modify this if we're too close to the edge
                        if (index == 0) {
                            index = 1;
                        } else if (index >= nr_fine_loc-1) {
                            index = nr_fine_loc-2;
                        }

                        s0macz(i,j,k) = QuadInterp(radius, 
                                            r_cc_loc_p[(index-1)*max_lev],
                                            r_cc_loc_p[index*max_lev], 
                                            r_cc_loc_p[(index+1)*max_lev], 
                                            s0_p[(index-1)*max_lev],
                                            s0_p[index*max_lev],
                                            s0_p[(index+1)*max_lev]);
                    });
                }

            } else { // use_exact_base_state = 0

                if (s0mac_interp_type == 1) {

                    AMREX_PARALLEL_FOR_3D(xbx, i, j, k, {
                        s0macx(i,j,k) = 0.5 * (s0_cart_arr(i-1,j,k) + s0_cart_arr(i,j,k));
                    });

                    AMREX_PARALLEL_FOR_3D(ybx, i, j, k, {
                        s0macy(i,j,k) = 0.5 * (s0_cart_arr(i,j-1,k) + s0_cart_arr(i,j,k));
                    });

                    AMREX_PARALLEL_FOR_3D(zbx, i, j, k, {
                        s0macz(i,j,k) = 0.5 * (s0_cart_arr(i,j,k-1) + s0_cart_arr(i,j,k));
                    });

                } else if (s0mac_interp_type == 2) {

                    AMREX_PARALLEL_FOR_3D(xbx, i, j, k, {
                        Real x = prob_lo[0] + Real(i) * dx[0] - center[0];
                        Real y = prob_lo[1] + (Real(j)+0.5) * dx[1] - center[1];
                        Real z = prob_lo[2] + (Real(k)+0.5) * dx[2] - center[2];

                        Real radius = sqrt(x*x + y*y + z*z);
                        int index = int(radius / dr);

                        if (radius >= r_cc_loc_p[index*max_lev]) {
                            if (index >= nr_fine_loc-1) {
                                s0macx(i,j,k) = s0_p[(nr_fine-1)*max_lev];
                            } else {
                                s0macx(i,j,k) = s0_p[(index+1)*max_lev] 
                                    * (radius-r_cc_loc_p[index*max_lev])/dr
                                    + s0_p[index*max_lev]
                                    * (r_cc_loc_p[(index+1)*max_lev]-radius)/dr;
                            }
                        } else {
                            if (index == 0) {
                                s0macx(i,j,k) = s0_p[index*max_lev];
                            } else if (index > nr_fine_loc-1) {
                                s0macx(i,j,k) = s0_p[(nr_fine-1)*max_lev];
                            } else {
                                s0macx(i,j,k) = s0_p[index*max_lev] 
                                    * (radius-r_cc_loc_p[(index-1)*max_lev])/dr
                                    + s0_p[(index-1)*max_lev]
                                    * (r_cc_loc_p[index*max_lev]-radius)/dr;
                            }
                        }
                    });

                    AMREX_PARALLEL_FOR_3D(ybx, i, j, k, {
                        Real x = prob_lo[0] + (Real(i)+0.5) * dx[0] - center[0];
                        Real y = prob_lo[1] + Real(j) * dx[1] - center[1];
                        Real z = prob_lo[2] + (Real(k)+0.5) * dx[2] - center[2];

                        Real radius = sqrt(x*x + y*y + z*z);
                        int index = int(radius / dr);

                        if (radius >= r_cc_loc_p[index*max_lev]) {
                            if (index >= nr_fine_loc-1) {
                                s0macy(i,j,k) = s0_p[(nr_fine-1)*max_lev];
                            } else {
                                s0macy(i,j,k) = s0_p[(index+1)*max_lev] 
                                    * (radius-r_cc_loc_p[index*max_lev])/dr
                                    + s0_p[index*max_lev]
                                    * (r_cc_loc_p[(index+1)*max_lev]-radius)/dr;
                            }
                        } else {
                            if (index == 0) {
                                s0macy(i,j,k) = s0_p[index*max_lev];
                            } else if (index > nr_fine_loc-1) {
                                s0macy(i,j,k) = s0_p[(nr_fine-1)*max_lev];
                            } else {
                                s0macy(i,j,k) = s0_p[index*max_lev] 
                                    * (radius-r_cc_loc_p[(index-1)*max_lev])/dr
                                    + s0_p[(index-1)*max_lev]
                                    * (r_cc_loc_p[index*max_lev]-radius)/dr;
                            }
                        }
                    });

                    AMREX_PARALLEL_FOR_3D(zbx, i, j, k, {
                        Real x = prob_lo[0] + (Real(i)+0.5) * dx[0] - center[0];
                        Real y = prob_lo[1] + (Real(j)+0.5) * dx[1] - center[1];
                        Real z = prob_lo[2] + Real(k) * dx[2] - center[2];

                        Real radius = sqrt(x*x + y*y + z*z);
                        int index = int(radius / dr);

                        if (radius >= r_cc_loc_p[index*max_lev]) {
                            if (index >= nr_fine_loc-1) {
                                s0macz(i,j,k) = s0_p[(nr_fine-1)*max_lev];
                            } else {
                                s0macz(i,j,k) = s0_p[(index+1)*max_lev] 
                                    * (radius-r_cc_loc_p[index*max_lev])/dr
                                    + s0_p[index*max_lev]
                                    * (r_cc_loc_p[(index+1)*max_lev]-radius)/dr;
                            }
                        } else {
                            if (index == 0) {
                                s0macz(i,j,k) = s0_p[index*max_lev];
                            } else if (index > nr_fine_loc-1) {
                                s0macz(i,j,k) = s0_p[(nr_fine-1)*max_lev];
                            } else {
                                s0macz(i,j,k) = s0_p[index*max_lev] 
                                    * (radius-r_cc_loc_p[(index-1)*max_lev])/dr
                                    + s0_p[(index-1)*max_lev]
                                    * (r_cc_loc_p[index*max_lev]-radius)/dr;
                            }
                        }
                    });

                } else if (s0mac_interp_type == 3) {

                    AMREX_PARALLEL_FOR_3D(xbx, i, j, k, {
                        Real x = prob_lo[0] + Real(i) * dx[0] - center[0];
                        Real y = prob_lo[1] + (Real(j)+0.5) * dx[1] - center[1];
                        Real z = prob_lo[2] + (Real(k)+0.5) * dx[2] - center[2];

                        Real radius = sqrt(x*x + y*y + z*z);
                        int index = int(radius / dr);

                        if (index == 0) {
                            index = 1;
                        } else if (index >= nr_fine_loc-1) {
                            index = nr_fine_loc-2;
                        }

                        s0macx(i,j,k) = QuadInterp(radius, 
                                            r_cc_loc_p[(index-1)*max_lev],
                                            r_cc_loc_p[index*max_lev], 
                                            r_cc_loc_p[(index+1)*max_lev], 
                                            s0_p[(index-1)*max_lev],
                                            s0_p[index*max_lev],
                                            s0_p[(index+1)*max_lev]);
                    });

                    AMREX_PARALLEL_FOR_3D(ybx, i, j, k, {
                        Real x = prob_lo[0] + (Real(i)+0.5) * dx[0] - center[0];
                        Real y = prob_lo[1] + Real(j) * dx[1] - center[1];
                        Real z = prob_lo[2] + (Real(k)+0.5) * dx[2] - center[2];

                        Real radius = sqrt(x*x + y*y + z*z);
                        int index = int(radius / dr);

                        if (index == 0) {
                            index = 1;
                        } else if (index >= nr_fine_loc-1) {
                            index = nr_fine_loc-2;
                        }

                        s0macy(i,j,k) = QuadInterp(radius, 
                                            r_cc_loc_p[(index-1)*max_lev],
                                            r_cc_loc_p[index*max_lev], 
                                            r_cc_loc_p[(index+1)*max_lev], 
                                            s0_p[(index-1)*max_lev],
                                            s0_p[index*max_lev],
                                            s0_p[(index+1)*max_lev]);
                    });

                    AMREX_PARALLEL_FOR_3D(zbx, i, j, k, {
                        Real x = prob_lo[0] + (Real(i)+0.5) * dx[0] - center[0];
                        Real y = prob_lo[1] + (Real(j)+0.5) * dx[1] - center[1];
                        Real z = prob_lo[2] + Real(k) * dx[2] - center[2];

                        Real radius = sqrt(x*x + y*y + z*z);
                        int index = int(radius / dr);

                        if (index == 0) {
                            index = 1;
                        } else if (index >= nr_fine_loc-1) {
                            index = nr_fine_loc-2;
                        }

                        s0macz(i,j,k) = QuadInterp(radius, 
                                            r_cc_loc_p[(index-1)*max_lev],
                                            r_cc_loc_p[index*max_lev], 
                                            r_cc_loc_p[(index+1)*max_lev], 
                                            s0_p[(index-1)*max_lev],
                                            s0_p[index*max_lev],
                                            s0_p[(index+1)*max_lev]);
                    });
                }
            }
        }
    }
}

void
Maestro::MakeNormal ()
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeNormal()", MakeNormal);

    // normal is the unit vector in the radial direction (e_r) in spherical
    // coordinates.
    // 
    // in terms of Cartesian coordinates, with unit vectors e_x, e_y, e_z,
    //    e_r = sin(theta)cos(phi) e_x + sin(theta)sin(phi) e_y + cos(theta) e_z
    // or
    //    e_r = (x/R) e_x + (y/R) e_y + (z/R) e_z

    if (spherical == 1) {

        for (int lev=0; lev<=finest_level; ++lev) {

            const auto dx = geom[lev].CellSizeArray();
            const auto prob_lo = geom[lev].ProbLoArray();
            GpuArray<Real,AMREX_SPACEDIM> center;

            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                center[n] = 0.5 * (geom[lev].ProbLo(n) + geom[lev].ProbHi(n));
            }

            // get references to the MultiFabs at level lev
            MultiFab& normal_mf = normal[lev];

            // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
    #ifdef _OPENMP
    #pragma omp parallel
    #endif
            for ( MFIter mfi(normal_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

                const Box& tileBox = mfi.tilebox();
                const Array4<Real> normal_arr = normal[lev].array(mfi);

                AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                    Real x = prob_lo[0] + (Real(i)+0.5) * dx[0] - center[0];
                    Real y = prob_lo[1] + (Real(j)+0.5) * dx[1] - center[1];
                    Real z = prob_lo[2] + (Real(k)+0.5) * dx[2] - center[2];

                    Real inv_radius = 1.0 / sqrt(x*x + y*y + z*z);

                    normal_arr(i,j,k,0) = x * inv_radius;
                    normal_arr(i,j,k,1) = y * inv_radius;
                    normal_arr(i,j,k,2) = z * inv_radius;
                });
            }
        }
    }
}


void
Maestro::PutDataOnFaces(const Vector<MultiFab>& s_cc,
                        Vector<std::array< MultiFab, AMREX_SPACEDIM > >& face,
                        const int harmonic_avg) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::PutDataOnFaces()",PutDataOnFaces);

    for (int lev=0; lev<=finest_level; ++lev) {
        
        // need one cell-centered MF for the MFIter
        const MultiFab& scc_mf = s_cc[lev];

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(scc_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

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
                AMREX_PARALLEL_FOR_3D(xbx, i, j, k, {
                    Real denom = scc(i,j,k) + scc(i-1,j,k);
                    Real prod  = scc(i,j,k) * scc(i-1,j,k);

                    if (denom != 0.0) {
                        facex(i,j,k) = 2.0 * prod / denom;
                    } else {
                        facex(i,j,k) = 0.5 * denom;
                    }
                });

                AMREX_PARALLEL_FOR_3D(ybx, i, j, k, {
                    Real denom = scc(i,j,k) + scc(i,j-1,k);
                    Real prod  = scc(i,j,k) * scc(i,j-1,k);

                    if (denom != 0.0) {
                        facey(i,j,k) = 2.0 * prod / denom;
                    } else {
                        facey(i,j,k) = 0.5 * denom;
                    }
                });
#if (AMREX_SPACEDIM == 3)
                AMREX_PARALLEL_FOR_3D(zbx, i, j, k, {
                    Real denom = scc(i,j,k) + scc(i,j,k-1);
                    Real prod  = scc(i,j,k) * scc(i,j,k-1);

                    if (denom != 0.0) {
                        facez(i,j,k) = 2.0 * prod / denom;
                    } else {
                        facez(i,j,k) = 0.5 * denom;
                    }
                });
#endif
            } else {
                AMREX_PARALLEL_FOR_3D(xbx, i, j, k, {
                    facex(i,j,k) = 0.5 * (scc(i,j,k) + scc(i-1,j,k));
                });
                AMREX_PARALLEL_FOR_3D(ybx, i, j, k, {
                    facey(i,j,k) = 0.5 * (scc(i,j,k) + scc(i,j-1,k));
                });
#if (AMREX_SPACEDIM == 3)
                AMREX_PARALLEL_FOR_3D(zbx, i, j, k, {
                    facez(i,j,k) = 0.5 * (scc(i,j,k) + scc(i,j,k-1));
                });
#endif
            }
        }
    }

    // Make sure that the fine edges average down onto the coarse edges (edge_restriction)
    AverageDownFaces(face);
}


void
Maestro::MakeCCtoRadii ()
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeCCtoRadius()",MakeCCtoRadii);

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const Real* dx = geom[lev].CellSize();
        const Real* dx_fine = geom[max_level].CellSize();

        MultiFab& cc_to_r = cell_cc_to_r[lev];

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(cc_to_r, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            const Box& tilebox = mfi.tilebox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
#pragma gpu box(tilebox)
            init_base_state_map_sphr(AMREX_INT_ANYD(tilebox.loVect()), 
                     AMREX_INT_ANYD(tilebox.hiVect()), 
                     BL_TO_FORTRAN_ANYD(cc_to_r[mfi]),
                                     AMREX_REAL_ANYD(dx_fine),
                                     AMREX_REAL_ANYD(dx));
        }
    }
}

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
    // MultiFab& cc_to_r = cell_cc_to_r[lev];
    GpuArray<Real,AMREX_SPACEDIM> dx;
    GpuArray<Real,AMREX_SPACEDIM> center;
    GpuArray<Real,AMREX_SPACEDIM> prob_lo;

    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        dx[n] = geom[lev].CellSize()[n];
        center[n] = 0.5 * (geom[lev].ProbLo(n) + geom[lev].ProbHi(n));
        prob_lo[n] = geom[lev].ProbLo(n);
    }

    Real * AMREX_RESTRICT r_edge_loc_p = r_edge_loc.dataPtr();
    Real * AMREX_RESTRICT r_cc_loc_p = r_cc_loc.dataPtr();

    // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
    for ( MFIter mfi(s0_cart_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

    	// Get the index space of the valid region
    	const Box& tileBox = mfi.tilebox();

        int max_lev = max_radial_level;
        int nr_fine_loc = nr_fine;

        const Array4<Real> s0_cart_arr = s0_cart[lev].array(mfi);
        const Real * AMREX_RESTRICT s0_p = s0.dataPtr();

    	// call fortran subroutine
    	// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
    	// lo/hi coordinates (including ghost cells), and/or the # of components
    	// We will also pass "validBox", which specifies the "valid" region.
    	if (spherical == 0) {

            AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {

                const int outcomp = is_output_a_vector == 1 ? AMREX_SPACEDIM-1 : 0;
                int r = AMREX_SPACEDIM == 2 ? j : k;

                s0_cart_arr(i,j,k,outcomp) = is_input_edge_centered == 1 ? 
                    0.5 * (s0_p[lev+r*(max_lev+1)] + s0_p[lev+(r+1)*(max_lev+1)]) : 
                    s0_p[lev+r*(max_lev+1)];
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
                    if (w0_interp_type == 1) {

                        AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {

                            Real x = prob_lo[0] + (Real(i)+0.5) * dx[0] - center[0];
                            Real y = prob_lo[1] + (Real(j)+0.5) * dx[1] - center[1];
                            Real z = prob_lo[2] + (Real(k)+0.5) * dx[2] - center[2];

                            Real radius = sqrt(x*x + y*y + z*z);
                            int index = cc_to_r(i,j,k);

                            Real rfac = (radius - r_edge_loc_p[(index+1)*(max_lev+1)]) 
                                / (r_cc_loc_p[(index+1)*(max_lev+1)] 
                                   - r_cc_loc_p[index*(max_lev+1)]);

                            Real s0_cart_val = rfac > 0.5 ? 
                                s0_p[(index+1)*(max_lev+1)] : s0_p[index*(max_lev+1)];
                            
                            if (is_output_a_vector) {
                                s0_cart_arr(i,j,k,0) = s0_cart_val * x / radius;
                                s0_cart_arr(i,j,k,1) = s0_cart_val * y / radius;
                                s0_cart_arr(i,j,k,2) = s0_cart_val * z / radius;
                            } else {
                                s0_cart_arr(i,j,k,0) = s0_cart_val;
                            }
                        });

                    } else if (w0_interp_type == 2) {

                        AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {

                            Real x = prob_lo[0] + (Real(i)+0.5) * dx[0] - center[0];
                            Real y = prob_lo[1] + (Real(j)+0.5) * dx[1] - center[1];
                            Real z = prob_lo[2] + (Real(k)+0.5) * dx[2] - center[2];

                            Real radius = sqrt(x*x + y*y + z*z);
                            int index = cc_to_r(i,j,k);

                            Real rfac = (radius - r_edge_loc_p[(index+1)*(max_lev+1)]) 
                                / (r_cc_loc_p[(index+1)*(max_lev+1)] 
                                   - r_cc_loc_p[index*(max_lev+1)]);

                            Real s0_cart_val;
                            if (index < nr_fine_loc) {
                                s0_cart_val = rfac * s0_p[(index+1)*(max_lev+1)] 
                                    + (1.0-rfac) * s0_p[index*(max_lev+1)];
                            } else {
                                s0_cart_val = s0_p[index*(max_lev+1)];
                            }

                            if (is_output_a_vector) {
                                s0_cart_arr(i,j,k,0) = s0_cart_val * x / radius;
                                s0_cart_arr(i,j,k,1) = s0_cart_val * y / radius;
                                s0_cart_arr(i,j,k,2) = s0_cart_val * z / radius;
                            } else {
                                s0_cart_arr(i,j,k,0) = s0_cart_val;
                            }
                        });

                    } else if (w0_interp_type == 3) {

                        AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {

                            Real x = prob_lo[0] + (Real(i)+0.5) * dx[0] - center[0];
                            Real y = prob_lo[1] + (Real(j)+0.5) * dx[1] - center[1];
                            Real z = prob_lo[2] + (Real(k)+0.5) * dx[2] - center[2];

                            Real radius = sqrt(x*x + y*y + z*z);
                            int index = cc_to_r(i,j,k) + 1;

                            if (index <= 0) {
                                index = 0;
                            } else if (index >= nr_fine_loc-1) {
                                index = nr_fine_loc - 2;
                            } else if (radius-r_edge_loc_p[index*(max_lev+1)] 
                                       < r_edge_loc_p[(index+1)*(max_lev+1)]) {
                                index--;
                            }

                            Real s0_cart_val = QuadInterp(radius, 
                                r_edge_loc_p[index*(max_lev+1)],
                                r_edge_loc_p[(index+1)*(max_lev+1)], 
                                r_edge_loc_p[(index+2)*(max_lev+1)], 
                                s0_p[index*(max_lev+1)],
                                s0_p[(index+1)*(max_lev+1)],
                                s0_p[(index+2)*(max_lev+1)]);

                            if (is_output_a_vector) {
                                s0_cart_arr(i,j,k,0) = s0_cart_val * x / radius;
                                s0_cart_arr(i,j,k,1) = s0_cart_val * y / radius;
                                s0_cart_arr(i,j,k,2) = s0_cart_val * z / radius;
                            } else {
                                s0_cart_arr(i,j,k,0) = s0_cart_val;
                            }
                        });
                    } else {
                        Abort("Error: w0_interp_type not defined");
                    }
                } else { // is_input_edge_centered = 0
                    // we directly inject the spherical values into each cell center
                    // because s0 is also bin-centered.

                    AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {

                        Real x = prob_lo[0] + (Real(i)+0.5) * dx[0] - center[0];
                        Real y = prob_lo[1] + (Real(j)+0.5) * dx[1] - center[1];
                        Real z = prob_lo[2] + (Real(k)+0.5) * dx[2] - center[2];

                        Real radius = sqrt(x*x + y*y + z*z);
                        int index = cc_to_r(i,j,k);

                        Real s0_cart_val = s0_p[index*(max_lev+1)];
                        
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

                    if (w0_interp_type == 1) {

                        AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {

                            Real x = prob_lo[0] + (Real(i)+0.5) * dx[0] - center[0];
                            Real y = prob_lo[1] + (Real(j)+0.5) * dx[1] - center[1];
                            Real z = prob_lo[2] + (Real(k)+0.5) * dx[2] - center[2];

                            Real radius = sqrt(x*x + y*y + z*z);
                            int index = int(radius / dr);

                            Real rfac = (radius - Real(index) * dr) / dr;

                            Real s0_cart_val = rfac > 0.5 ? 
                                s0_p[(index+1)*(max_lev+1)] : s0_p[index*(max_lev+1)];
                            
                            if (is_output_a_vector) {
                                s0_cart_arr(i,j,k,0) = s0_cart_val * x / radius;
                                s0_cart_arr(i,j,k,1) = s0_cart_val * y / radius;
                                s0_cart_arr(i,j,k,2) = s0_cart_val * z / radius;
                            } else {
                                s0_cart_arr(i,j,k,0) = s0_cart_val;
                            }
                        });

                    } else if (w0_interp_type == 2) {

                        AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {

                            Real x = prob_lo[0] + (Real(i)+0.5) * dx[0] - center[0];
                            Real y = prob_lo[1] + (Real(j)+0.5) * dx[1] - center[1];
                            Real z = prob_lo[2] + (Real(k)+0.5) * dx[2] - center[2];

                            Real radius = sqrt(x*x + y*y + z*z);
                            int index = int(radius / dr);

                            Real rfac = (radius - Real(index) * dr) / dr;

                            Real s0_cart_val;
                            if (index < nr_fine_loc) {
                                s0_cart_val = rfac * s0_p[(index+1)*(max_lev+1)] + (1.0-rfac) * s0_p[index*(max_lev+1)];
                            } else {
                                s0_cart_val = s0_p[nr_fine_loc*(max_lev+1)];
                            }

                            if (is_output_a_vector) {
                                s0_cart_arr(i,j,k,0) = s0_cart_val * x / radius;
                                s0_cart_arr(i,j,k,1) = s0_cart_val * y / radius;
                                s0_cart_arr(i,j,k,2) = s0_cart_val * z / radius;
                            } else {
                                s0_cart_arr(i,j,k,0) = s0_cart_val;
                            }
                        });

                    } else if (w0_interp_type == 3) {

                        AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {

                            Real x = prob_lo[0] + (Real(i)+0.5) * dx[0] - center[0];
                            Real y = prob_lo[1] + (Real(j)+0.5) * dx[1] - center[1];
                            Real z = prob_lo[2] + (Real(k)+0.5) * dx[2] - center[2];

                            Real radius = sqrt(x*x + y*y + z*z);
                            int index = int(radius / dr);

                            if (index <= 0) {
                                index = 0;
                            } else if (index >= nr_fine_loc-1) {
                                index = nr_fine_loc - 2;
                            } else if (radius-r_edge_loc_p[index*(max_lev+1)] 
                                       < r_edge_loc_p[(index+1)*(max_lev+1)]) {
                                index--;
                            }

                            Real s0_cart_val = QuadInterp(radius, 
                                r_edge_loc_p[index*(max_lev+1)],
                                r_edge_loc_p[(index+1)*(max_lev+1)], 
                                r_edge_loc_p[(index+2)*(max_lev+1)], 
                                s0_p[index*(max_lev+1)],
                                s0_p[(index+1)*(max_lev+1)],
                                s0_p[(index+2)*(max_lev+1)]);

                            if (is_output_a_vector) {
                                s0_cart_arr(i,j,k,0) = s0_cart_val * x / radius;
                                s0_cart_arr(i,j,k,1) = s0_cart_val * y / radius;
                                s0_cart_arr(i,j,k,2) = s0_cart_val * z / radius;
                            } else {
                                s0_cart_arr(i,j,k,0) = s0_cart_val;
                            }
                        });
                    } else {
                        Abort("Error: w0_interp_type not defined");
                    }
                } else { // is_input_edge_centered = 0
                    // we currently have three different ideas for computing s0_cart,
                    // where s0 is bin-centered.
                    // 1.  Piecewise constant
                    // 2.  Piecewise linear
                    // 3.  Quadratic

                    if (s0_interp_type == 1) {

                        AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {

                            Real x = prob_lo[0] + (Real(i)+0.5) * dx[0] - center[0];
                            Real y = prob_lo[1] + (Real(j)+0.5) * dx[1] - center[1];
                            Real z = prob_lo[2] + (Real(k)+0.5) * dx[2] - center[2];

                            Real radius = sqrt(x*x + y*y + z*z);
                            int index = int(radius / dr);

                            Real s0_cart_val = s0_p[index*(max_lev+1)];
                            
                            if (is_output_a_vector) {
                                s0_cart_arr(i,j,k,0) = s0_cart_val * x / radius;
                                s0_cart_arr(i,j,k,1) = s0_cart_val * y / radius;
                                s0_cart_arr(i,j,k,2) = s0_cart_val * z / radius;
                            } else {
                                s0_cart_arr(i,j,k,0) = s0_cart_val;
                            }
                        });

                    } else if (s0_interp_type == 2) {

                        AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {

                            Real x = prob_lo[0] + (Real(i)+0.5) * dx[0] - center[0];
                            Real y = prob_lo[1] + (Real(j)+0.5) * dx[1] - center[1];
                            Real z = prob_lo[2] + (Real(k)+0.5) * dx[2] - center[2];

                            Real radius = sqrt(x*x + y*y + z*z);
                            int index = int(radius / dr);

                            Real s0_cart_val;
                            
                            if (radius >= r_cc_loc_p[index*(max_lev+1)]) {
                                if (index >= nr_fine_loc-1) {
                                    s0_cart_val = s0_p[(nr_fine_loc-1)*(max_lev+1)];
                                } else {
                                    s0_cart_val = s0_p[(index+1)*(max_lev+1)] 
                                        * (radius-r_cc_loc_p[index*(max_lev+1)])/dr 
                                        + s0_p[index*(max_lev+1)] 
                                        * (r_cc_loc_p[(index+1)*(max_lev+1)]-radius)/dr;
                                }
                            } else {
                                if (index == 0) {
                                    s0_cart_val = s0_p[index*(max_lev+1)];
                                } else if (index > nr_fine_loc-1) {
                                    s0_cart_val = s0_p[(nr_fine_loc-1)*(max_lev+1)];
                                } else {
                                    s0_cart_val = s0_p[index*(max_lev+1)] 
                                        * (radius-r_cc_loc_p[(index-1)*(max_lev+1)])/dr 
                                        + s0_p[(index-1)*(max_lev+1)] 
                                        * (r_cc_loc_p[index*(max_lev+1)]-radius)/dr;
                                }
                            }
                            
                            if (is_output_a_vector) {
                                s0_cart_arr(i,j,k,0) = s0_cart_val * x / radius;
                                s0_cart_arr(i,j,k,1) = s0_cart_val * y / radius;
                                s0_cart_arr(i,j,k,2) = s0_cart_val * z / radius;
                            } else {
                                s0_cart_arr(i,j,k,0) = s0_cart_val;
                            }
                        });

                    } else if (s0_interp_type == 3) {

                        AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {

                            Real x = prob_lo[0] + (Real(i)+0.5) * dx[0] - center[0];
                            Real y = prob_lo[1] + (Real(j)+0.5) * dx[1] - center[1];
                            Real z = prob_lo[2] + (Real(k)+0.5) * dx[2] - center[2];

                            Real radius = sqrt(x*x + y*y + z*z);
                            int index = int(radius / dr);

                            if (index == 0) {
                                index = 1;
                            } else if (index >= nr_fine_loc-1) {
                                index = nr_fine_loc-2;
                            }

                            Real s0_cart_val = QuadInterp(radius, 
                                r_edge_loc_p[(index-1)*(max_lev+1)],
                                r_edge_loc_p[index*(max_lev+1)],
                                r_edge_loc_p[(index+1)*(max_lev+1)], 
                                s0_p[(index-1)*(max_lev+1)],
                                s0_p[index*(max_lev+1)],
                                s0_p[(index+1)*(max_lev+1)]);
                            
                            if (is_output_a_vector) {
                                s0_cart_arr(i,j,k,0) = s0_cart_val * x / radius;
                                s0_cart_arr(i,j,k,1) = s0_cart_val * y / radius;
                                s0_cart_arr(i,j,k,2) = s0_cart_val * z / radius;
                            } else {
                                s0_cart_arr(i,j,k,0) = s0_cart_val;
                            }
                        });
                    } // s0_interp_type
                } // is_input_edge_centere
            } // use_exact_base_state
    	}
    }
}

AMREX_GPU_DEVICE 
Real
Maestro::QuadInterp(const Real x, const Real x0, const Real x1, const Real x2,
                    const Real y0, const Real y1, const Real y2) 
{
    Real y = y0 + (y1-y0)/(x1-x0)*(x-x0) 
         + ((y2-y1)/(x2-x1)-(y1-y0)/(x1-x0))/(x2-x0)*(x-x0)*(x-x1);

    if (y > max(y0, max(y1, y2))) y = max(y0, max(y1, y2));
    if (y < min(y0, min(y1, y2))) y = min(y0, min(y1, y2));

    return y;
}


void
Maestro::Addw0 (Vector<std::array< MultiFab, AMREX_SPACEDIM > >& uedge,
                const Vector<std::array< MultiFab,AMREX_SPACEDIM > >& w0mac,
                const Real& mult)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Addw0()",Addw0);

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& uedge_mf = uedge[lev][0];
        MultiFab& vedge_mf = uedge[lev][1];
#if (AMREX_SPACEDIM == 3)
        MultiFab& wedge_mf = uedge[lev][2];
        const MultiFab& w0macx_mf = w0mac[lev][0];
        const MultiFab& w0macy_mf = w0mac[lev][1];
        const MultiFab& w0macz_mf = w0mac[lev][2];
#endif

        // need one cell-centered MF for the MFIter
        MultiFab& sold_mf = sold[lev];

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(sold_mf); mfi.isValid(); ++mfi ) {
            
            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            if (spherical == 0) {

#if (AMREX_SPACEDIM == 2)
                const Box& ybx = amrex::grow(mfi.nodaltilebox(1), amrex::IntVect(1,0));
#else
                const Box& zbx = amrex::grow(mfi.nodaltilebox(2), amrex::IntVect(1,1,0));
#endif

#if (AMREX_SPACEDIM == 2)
#pragma gpu box(ybx)
                addw0(AMREX_INT_ANYD(ybx.loVect()), 
                      AMREX_INT_ANYD(ybx.hiVect()), 
                      lev,
                      BL_TO_FORTRAN_ANYD(vedge_mf[mfi]),
                      w0.dataPtr(),mult);
#else 
#pragma gpu box(zbx)
                addw0(AMREX_INT_ANYD(zbx.loVect()), 
                      AMREX_INT_ANYD(zbx.hiVect()), 
                      lev,
                      BL_TO_FORTRAN_ANYD(wedge_mf[mfi]),
                      w0.dataPtr(),mult);
#endif

            } else {

#if (AMREX_SPACEDIM == 3)

                const Box& xbx = mfi.nodaltilebox(0); 
                const Box& ybx = mfi.nodaltilebox(1); 
                const Box& zbx = mfi.nodaltilebox(2); 
                
#pragma gpu box(xbx)
                addw0_sphr(AMREX_INT_ANYD(xbx.loVect()), 
                           AMREX_INT_ANYD(xbx.hiVect()),
                           BL_TO_FORTRAN_ANYD(uedge_mf[mfi]),
                           BL_TO_FORTRAN_ANYD(w0macx_mf[mfi]),
                           mult);
#pragma gpu box(ybx)
                addw0_sphr(AMREX_INT_ANYD(ybx.loVect()), 
                           AMREX_INT_ANYD(ybx.hiVect()),
                           BL_TO_FORTRAN_ANYD(vedge_mf[mfi]),
                           BL_TO_FORTRAN_ANYD(w0macy_mf[mfi]),
                           mult);
#pragma gpu box(zbx)
                addw0_sphr(AMREX_INT_ANYD(zbx.loVect()), 
                           AMREX_INT_ANYD(zbx.hiVect()),
                           BL_TO_FORTRAN_ANYD(wedge_mf[mfi]),
                           BL_TO_FORTRAN_ANYD(w0macz_mf[mfi]),
                           mult);
#else
                Abort("Addw0: Spherical is not valid for DIM < 3");
#endif
            }             //end spherical

        }
    }

    if (finest_level == 0) {
        // fill periodic ghost cells
        for (int lev=0; lev<=finest_level; ++lev) {
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                uedge[lev][d].FillBoundary(geom[lev].periodicity());
            }
        }

        // fill ghost cells behind physical boundaries
        FillUmacGhost(uedge);
    } else {
        // edge_restriction
        AverageDownFaces(uedge);

        // fill all ghost cells for edge-based velocity field
        FillPatchUedge(uedge);
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

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& w0macx_mf = w0mac[lev][0];
        MultiFab& w0macy_mf = w0mac[lev][1];
        MultiFab& w0macz_mf = w0mac[lev][2];
        MultiFab& w0cart_mf = w0_cart[lev];
        MultiFab& w0_nodal_mf = w0_nodal[lev];
        const Real* dx = geom[lev].CellSize();

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(w0cart_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& gntbx = mfi.grownnodaltilebox(-1, 1);

            if (w0mac_interp_type == 4) {
#pragma gpu box(gntbx)
                make_w0mac_nodal(AMREX_INT_ANYD(gntbx.loVect()), 
                            AMREX_INT_ANYD(gntbx.hiVect()),
                            w0.dataPtr(),
                            BL_TO_FORTRAN_ANYD(w0_nodal_mf[mfi]), 
                            AMREX_REAL_ANYD(dx));

            }
        }

#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(w0cart_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            const Box& xbx = mfi.grownnodaltilebox(0, 1);
            const Box& ybx = mfi.grownnodaltilebox(1, 1);
            const Box& zbx = mfi.grownnodaltilebox(2, 1);

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
#pragma gpu box(xbx)
            make_w0mac_sphr(AMREX_INT_ANYD(xbx.loVect()), 
                            AMREX_INT_ANYD(xbx.hiVect()),1,
                            w0.dataPtr(),
                            BL_TO_FORTRAN_ANYD(w0macx_mf[mfi]),
                            BL_TO_FORTRAN_ANYD(w0cart_mf[mfi]), 
                            w0cart_mf.nComp(),
                            BL_TO_FORTRAN_ANYD(w0_nodal_mf[mfi]),
                            AMREX_REAL_ANYD(dx), r_edge_loc.dataPtr());

#pragma gpu box(ybx)
            make_w0mac_sphr(AMREX_INT_ANYD(ybx.loVect()), 
                            AMREX_INT_ANYD(ybx.hiVect()),2,
                            w0.dataPtr(),
                            BL_TO_FORTRAN_ANYD(w0macy_mf[mfi]),
                            BL_TO_FORTRAN_ANYD(w0cart_mf[mfi]), 
                            w0cart_mf.nComp(),
                            BL_TO_FORTRAN_ANYD(w0_nodal_mf[mfi]),
                            AMREX_REAL_ANYD(dx), r_edge_loc.dataPtr());

#pragma gpu box(zbx)
            make_w0mac_sphr(AMREX_INT_ANYD(zbx.loVect()), 
                            AMREX_INT_ANYD(zbx.hiVect()),3,
                            w0.dataPtr(),
                            BL_TO_FORTRAN_ANYD(w0macz_mf[mfi]),
                            BL_TO_FORTRAN_ANYD(w0cart_mf[mfi]), 
                            w0cart_mf.nComp(),
                            BL_TO_FORTRAN_ANYD(w0_nodal_mf[mfi]),
                            AMREX_REAL_ANYD(dx), r_edge_loc.dataPtr());
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

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& s0macx_mf = s0mac[lev][0];
        MultiFab& s0macy_mf = s0mac[lev][1];
        MultiFab& s0macz_mf = s0mac[lev][2];
        MultiFab& s0cart_mf = s0_cart[lev];
        const Real* dx = geom[lev].CellSize();

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(s0cart_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& xbx = mfi.grownnodaltilebox(0, 1);
            const Box& ybx = mfi.grownnodaltilebox(1, 1);
            const Box& zbx = mfi.grownnodaltilebox(2, 1);

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            if (use_exact_base_state) {
#pragma gpu box(xbx)
                make_s0mac_sphr_irreg(AMREX_INT_ANYD(xbx.loVect()), 
                            AMREX_INT_ANYD(xbx.hiVect()),1,
                            s0.dataPtr(),
                            BL_TO_FORTRAN_ANYD(s0macx_mf[mfi]),
                            BL_TO_FORTRAN_ANYD(s0cart_mf[mfi]),
                            AMREX_REAL_ANYD(dx), r_cc_loc.dataPtr());
#pragma gpu box(ybx)
                make_s0mac_sphr_irreg(AMREX_INT_ANYD(ybx.loVect()), 
                            AMREX_INT_ANYD(ybx.hiVect()),2,
                            s0.dataPtr(),
                            BL_TO_FORTRAN_ANYD(s0macy_mf[mfi]),
                            BL_TO_FORTRAN_ANYD(s0cart_mf[mfi]),
                            AMREX_REAL_ANYD(dx), r_cc_loc.dataPtr());
#pragma gpu box(zbx)
                make_s0mac_sphr_irreg(AMREX_INT_ANYD(zbx.loVect()), 
                            AMREX_INT_ANYD(zbx.hiVect()),3,
                            s0.dataPtr(),
                            BL_TO_FORTRAN_ANYD(s0macz_mf[mfi]),
                            BL_TO_FORTRAN_ANYD(s0cart_mf[mfi]),
                            AMREX_REAL_ANYD(dx), r_cc_loc.dataPtr());
            } else {

#pragma gpu box(xbx)
                make_s0mac_sphr(AMREX_INT_ANYD(xbx.loVect()), 
                        AMREX_INT_ANYD(xbx.hiVect()),1,
                        s0.dataPtr(),
                        BL_TO_FORTRAN_ANYD(s0macx_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(s0cart_mf[mfi]),
                        AMREX_REAL_ANYD(dx), r_cc_loc.dataPtr());
#pragma gpu box(ybx)   
                make_s0mac_sphr(AMREX_INT_ANYD(ybx.loVect()), 
                        AMREX_INT_ANYD(ybx.hiVect()),2,
                        s0.dataPtr(),
                        BL_TO_FORTRAN_ANYD(s0macy_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(s0cart_mf[mfi]),
                        AMREX_REAL_ANYD(dx), r_cc_loc.dataPtr());
#pragma gpu box(zbx)
                make_s0mac_sphr(AMREX_INT_ANYD(zbx.loVect()), 
                        AMREX_INT_ANYD(zbx.hiVect()),3,
                        s0.dataPtr(),
                        BL_TO_FORTRAN_ANYD(s0macz_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(s0cart_mf[mfi]),
                        AMREX_REAL_ANYD(dx), r_cc_loc.dataPtr());
            }
        }
    }
}


void
Maestro::MakeNormal ()
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeNormal()",MakeNormal);

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& normal_mf = normal[lev];
        const Real* dx = geom[lev].CellSize();

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(normal_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            const Box& tileBox = mfi.tilebox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
#pragma gpu box(tileBox)
            make_normal(AMREX_INT_ANYD(tileBox.loVect()),
                        AMREX_INT_ANYD(tileBox.hiVect()),
                        BL_TO_FORTRAN_ANYD(normal_mf[mfi]), 
                        AMREX_REAL_ANYD(dx));
        }
    }
}


void
Maestro::PutDataOnFaces(const Vector<MultiFab>& s_cc,
                        Vector<std::array< MultiFab, AMREX_SPACEDIM > >& face,
                        int harmonic_avg) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::PutDataOnFaces()",PutDataOnFaces);

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& facex_mf = face[lev][0];
        MultiFab& facey_mf = face[lev][1];
#if (AMREX_SPACEDIM == 3)
        MultiFab& facez_mf = face[lev][2];
#endif
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

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // x-direction
#pragma gpu box(xbx)
            put_data_on_faces(AMREX_INT_ANYD(xbx.loVect()), 
                              AMREX_INT_ANYD(xbx.hiVect()),1,
                              BL_TO_FORTRAN_ANYD(scc_mf[mfi]),
                              BL_TO_FORTRAN_ANYD(facex_mf[mfi]),
                              harmonic_avg);
            // y-direction
#pragma gpu box(ybx)
            put_data_on_faces(AMREX_INT_ANYD(ybx.loVect()), 
                              AMREX_INT_ANYD(ybx.hiVect()),2,
                              BL_TO_FORTRAN_ANYD(scc_mf[mfi]),
                              BL_TO_FORTRAN_ANYD(facey_mf[mfi]),
                              harmonic_avg);
            // z-direction
#if (AMREX_SPACEDIM == 3)
#pragma gpu box(zbx)
            put_data_on_faces(AMREX_INT_ANYD(zbx.loVect()), 
                              AMREX_INT_ANYD(zbx.hiVect()),3,
                              BL_TO_FORTRAN_ANYD(scc_mf[mfi]),
                              BL_TO_FORTRAN_ANYD(facez_mf[mfi]),
                              harmonic_avg);
#endif
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

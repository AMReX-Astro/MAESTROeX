
#include <Maestro.H>
#include <MaestroBCThreads.H>

using namespace amrex;

void
Maestro::MakeUtrans (Vector<MultiFab>& utilde,
                     const Vector<MultiFab>& ufull,
                     Vector<std::array< MultiFab, AMREX_SPACEDIM > >& utrans,
                     const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& w0mac)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeUtrans()",MakeUtrans);

    for (int lev=0; lev<=finest_level; ++lev) {

        // Get the index space and grid spacing of the domain
        const Box& domainBox = geom[lev].Domain();
        const Real* dx = geom[lev].CellSize();

        // get references to the MultiFabs at level lev
              MultiFab& utilde_mf  = utilde[lev];
        const MultiFab& ufull_mf   = ufull[lev];
              MultiFab& utrans_mf  = utrans[lev][0];
              MultiFab& vtrans_mf  = utrans[lev][1];
              MultiFab Ip, Im;
              Ip.define(grids[lev],dmap[lev],AMREX_SPACEDIM,2);
              Im.define(grids[lev],dmap[lev],AMREX_SPACEDIM,2);
#if (AMREX_SPACEDIM == 3)
              MultiFab& wtrans_mf  = utrans[lev][2];
        const MultiFab& w0macx_mf  = w0mac[lev][0];
        const MultiFab& w0macy_mf  = w0mac[lev][1];
        const MultiFab& w0macz_mf  = w0mac[lev][2];
#endif
        const MultiFab& w0_mf = w0_cart[lev];

#ifdef AMREX_USE_CUDA
        int* bc_f = prepare_bc(bcs_u[0].data(), AMREX_SPACEDIM);
#else
        const int* bc_f = bcs_u[0].data();
#endif
        MultiFab u_mf, v_mf, w_mf;

        if (ppm_type == 0) {
           u_mf.define(grids[lev],dmap[lev],1,utilde[lev].nGrow());
           v_mf.define(grids[lev],dmap[lev],1,utilde[lev].nGrow());

           MultiFab::Copy(u_mf, utilde[lev], 0, 0, 1, utilde[lev].nGrow());
           MultiFab::Copy(v_mf, utilde[lev], 1, 0, 1, utilde[lev].nGrow());

#if (AMREX_SPACEDIM == 3)
           w_mf.define(grids[lev],dmap[lev],1,utilde[lev].nGrow());
           MultiFab::Copy(w_mf, utilde[lev], 2, 0, 1, utilde[lev].nGrow());
#endif

        } else if (ppm_type == 1 || ppm_type == 2) {

           u_mf.define(grids[lev],dmap[lev],1,ufull[lev].nGrow());
           v_mf.define(grids[lev],dmap[lev],1,ufull[lev].nGrow());

           MultiFab::Copy(u_mf, ufull[lev], 0, 0, 1, ufull[lev].nGrow());
           MultiFab::Copy(v_mf, ufull[lev], 1, 0, 1, ufull[lev].nGrow());

#if (AMREX_SPACEDIM == 3)
           w_mf.define(grids[lev],dmap[lev],1,ufull[lev].nGrow());
           MultiFab::Copy(w_mf, ufull[lev], 2, 0, 1, ufull[lev].nGrow());
#endif
        }

#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(utilde_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const Box& obx = amrex::grow(tileBox, 1);
            const Box& xbx = mfi.nodaltilebox(0);
            const Box& ybx = mfi.nodaltilebox(1);
#if (AMREX_SPACEDIM == 3)
            const Box& zbx = mfi.nodaltilebox(2);
#endif

#if (AMREX_SPACEDIM == 2)

            if (ppm_type == 0) {
                // we're going to reuse Ip here as slopex as it has the
                // correct number of ghost zones
                // x-direction
                Slopex(obx, u_mf.array(mfi), 
                       Ip.array(mfi), 
                       domainBox, bcs_u, 
                       1,0);
            } else {

                PPM_2d(obx, utilde_mf.array(mfi), 
                       u_mf.array(mfi), v_mf.array(mfi), 
                       Ip.array(mfi), Im.array(mfi), 
                       domainBox, bcs_u, dx, 
                       false, 0, 0);
           }

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
#pragma gpu box(xbx)
            mkutrans_2d(AMREX_INT_ANYD(xbx.loVect()),
                        AMREX_INT_ANYD(xbx.hiVect()),
                        1,
                        AMREX_INT_ANYD(domainBox.loVect()),
                        AMREX_INT_ANYD(domainBox.hiVect()),
                        BL_TO_FORTRAN_ANYD(utilde_mf[mfi]), utilde_mf.nComp(), utilde_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(ufull_mf[mfi]), ufull_mf.nComp(), ufull_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(utrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(vtrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(Ip[mfi]),
                        BL_TO_FORTRAN_ANYD(Im[mfi]),
                        BL_TO_FORTRAN_ANYD(w0_mf[mfi]), 
                        AMREX_REAL_ANYD(dx), dt, bc_f,
                        phys_bc.dataPtr());

            if (ppm_type == 0) {
                // we're going to reuse Im here as slopey as it has the
                // correct number of ghost zones
                Slopey(obx, v_mf.array(mfi), 
                       Im.array(mfi), 
                       domainBox, bcs_u, 
                       1,1);
            } else {
                PPM_2d(obx, utilde_mf.array(mfi), 
                       u_mf.array(mfi), v_mf.array(mfi), 
                       Ip.array(mfi), Im.array(mfi), 
                       domainBox, bcs_u, dx, 
                       false, 1, 1);
           }

#pragma gpu box(ybx)
           mkutrans_2d(AMREX_INT_ANYD(ybx.loVect()),
                       AMREX_INT_ANYD(ybx.hiVect()),
                       2,
                       AMREX_INT_ANYD(domainBox.loVect()),
                       AMREX_INT_ANYD(domainBox.hiVect()),
                       BL_TO_FORTRAN_ANYD(utilde_mf[mfi]), utilde_mf.nComp(), utilde_mf.nGrow(),
                       BL_TO_FORTRAN_ANYD(ufull_mf[mfi]), ufull_mf.nComp(), ufull_mf.nGrow(),
                       BL_TO_FORTRAN_ANYD(utrans_mf[mfi]),
                       BL_TO_FORTRAN_ANYD(vtrans_mf[mfi]),
                       BL_TO_FORTRAN_ANYD(Ip[mfi]),
                       BL_TO_FORTRAN_ANYD(Im[mfi]),
                       BL_TO_FORTRAN_ANYD(w0_mf[mfi]), 
                       AMREX_REAL_ANYD(dx), dt, bc_f,
                       phys_bc.dataPtr());
#elif (AMREX_SPACEDIM == 3)

            // x-direction
            if (ppm_type == 0) {
                // we're going to reuse Ip here as slopex as it has the
                // correct number of ghost zones
                Slopex(obx, u_mf.array(mfi), 
                       Ip.array(mfi), 
                       domainBox, bcs_u, 
                       1,0);
            } else {
                PPM_3d(obx, utilde_mf.array(mfi), 
                       u_mf.array(mfi), v_mf.array(mfi), w_mf.array(mfi),
                       Ip.array(mfi), Im.array(mfi), 
                       domainBox, bcs_u, dx, 
                       false, 0, 0);
            }

#pragma gpu box(xbx)
            mkutrans_3d(AMREX_INT_ANYD(xbx.loVect()),
                        AMREX_INT_ANYD(xbx.hiVect()),
                        1,
                        AMREX_INT_ANYD(domainBox.loVect()),
                        AMREX_INT_ANYD(domainBox.hiVect()),
                        BL_TO_FORTRAN_ANYD(utilde_mf[mfi]), utilde_mf.nComp(), utilde_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(ufull_mf[mfi]), ufull_mf.nComp(), ufull_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(utrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(vtrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(wtrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(Ip[mfi]),
                        BL_TO_FORTRAN_ANYD(Im[mfi]),
                        BL_TO_FORTRAN_ANYD(w0macx_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(w0macy_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(w0macz_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(w0_mf[mfi]), 
                        AMREX_REAL_ANYD(dx), dt, bc_f, phys_bc.dataPtr());

            // y-direction
            if (ppm_type == 0) {
                // we're going to reuse Im here as slopey as it has the
                // correct number of ghost zones
                Slopey(obx, v_mf.array(mfi), 
                       Im.array(mfi), 
                       domainBox, bcs_u, 
                       1,1);
            } else {
                PPM_3d(obx, utilde_mf.array(mfi), 
                       u_mf.array(mfi), v_mf.array(mfi), w_mf.array(mfi),
                       Ip.array(mfi), Im.array(mfi), 
                       domainBox, bcs_u, dx, 
                       false, 1, 1);
            }

#pragma gpu box(ybx)
            mkutrans_3d(AMREX_INT_ANYD(ybx.loVect()),
                        AMREX_INT_ANYD(ybx.hiVect()),
                        2,
                        AMREX_INT_ANYD(domainBox.loVect()),
                        AMREX_INT_ANYD(domainBox.hiVect()),
                        BL_TO_FORTRAN_ANYD(utilde_mf[mfi]), utilde_mf.nComp(), utilde_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(ufull_mf[mfi]), ufull_mf.nComp(), ufull_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(utrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(vtrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(wtrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(Ip[mfi]),
                        BL_TO_FORTRAN_ANYD(Im[mfi]),
                        BL_TO_FORTRAN_ANYD(w0macx_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(w0macy_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(w0macz_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(w0_mf[mfi]), 
                        AMREX_REAL_ANYD(dx), dt, bc_f, phys_bc.dataPtr());

            // z-direction
            if (ppm_type == 0) {
                // we're going to reuse Im here as slopez as it has the
                // correct number of ghost zones
                Slopez(obx, w_mf.array(mfi), 
                       Im.array(mfi), 
                       domainBox, bcs_u,  
                       1,2);
            } else {
                PPM_3d(obx, utilde_mf.array(mfi), 
                       u_mf.array(mfi), v_mf.array(mfi), w_mf.array(mfi),
                       Ip.array(mfi), Im.array(mfi), 
                       domainBox, bcs_u, dx, 
                       false, 2, 2);
            }

#pragma gpu box(zbx)
            mkutrans_3d(AMREX_INT_ANYD(zbx.loVect()),
                        AMREX_INT_ANYD(zbx.hiVect()),
                        3,
                        AMREX_INT_ANYD(domainBox.loVect()),
                        AMREX_INT_ANYD(domainBox.hiVect()),
                        BL_TO_FORTRAN_ANYD(utilde_mf[mfi]), utilde_mf.nComp(), utilde_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(ufull_mf[mfi]), ufull_mf.nComp(), ufull_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(utrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(vtrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(wtrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(Ip[mfi]),
                        BL_TO_FORTRAN_ANYD(Im[mfi]),
                        BL_TO_FORTRAN_ANYD(w0macx_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(w0macy_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(w0macz_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(w0_mf[mfi]), 
                        AMREX_REAL_ANYD(dx), dt, bc_f, phys_bc.dataPtr());
#endif
        } // end MFIter loop

#ifdef AMREX_USE_CUDA
        clean_bc(bc_f);
#endif
    } // end loop over levels

    if (finest_level == 0) {
        // fill periodic ghost cells
        for (int lev=0; lev<=finest_level; ++lev) {
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                    utrans[lev][d].FillBoundary(geom[lev].periodicity());
            }
        }

        // fill ghost cells behind physical boundaries
        FillUmacGhost(utrans);
    } else {
        // edge_restriction
        AverageDownFaces(utrans);

        // fill ghost cells for all levels
        FillPatchUedge(utrans);
    }
}


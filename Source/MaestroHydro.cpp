
#include <Maestro.H>
#include <MaestroHydro_F.H>
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

// #pragma gpu box(obx)
//                 ppm_2d(AMREX_INT_ANYD(obx.loVect()),
//                        AMREX_INT_ANYD(obx.hiVect()),
//                        BL_TO_FORTRAN_ANYD(utilde_mf[mfi]),
//                        utilde_mf.nComp(),
//                        BL_TO_FORTRAN_ANYD(u_mf[mfi]),
//                        BL_TO_FORTRAN_ANYD(v_mf[mfi]),
//                        BL_TO_FORTRAN_ANYD(Ip[mfi]),
//                        BL_TO_FORTRAN_ANYD(Im[mfi]),
//                        AMREX_INT_ANYD(domainBox.loVect()),
//                        AMREX_INT_ANYD(domainBox.hiVect()),
//                        bc_f, AMREX_REAL_ANYD(dx), dt, false,
//                        1,1,AMREX_SPACEDIM);
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
                // PPM_2d(obx, utilde_mf.array(mfi), 
                //        u_mf.array(mfi), v_mf.array(mfi), 
                //        Ip.array(mfi), Im.array(mfi), 
                //        domainBox, bcs_u, dx, 
                //        false, 1, 1);

#pragma gpu box(obx)
                ppm_2d(AMREX_INT_ANYD(obx.loVect()),
                       AMREX_INT_ANYD(obx.hiVect()),
                       BL_TO_FORTRAN_ANYD(utilde_mf[mfi]),
                       utilde_mf.nComp(),
                       BL_TO_FORTRAN_ANYD(u_mf[mfi]),
                       BL_TO_FORTRAN_ANYD(v_mf[mfi]),
                       BL_TO_FORTRAN_ANYD(Ip[mfi]),
                       BL_TO_FORTRAN_ANYD(Im[mfi]),
                       AMREX_INT_ANYD(domainBox.loVect()),
                       AMREX_INT_ANYD(domainBox.hiVect()),
                       bc_f, AMREX_REAL_ANYD(dx), dt, false,
                       2,2,AMREX_SPACEDIM);
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
// #pragma gpu box(obx)
//                 ppm_3d(AMREX_INT_ANYD(obx.loVect()),
//                        AMREX_INT_ANYD(obx.hiVect()),
//                        BL_TO_FORTRAN_ANYD(utilde_mf[mfi]),
//                        utilde_mf.nComp(),
//                        BL_TO_FORTRAN_ANYD(u_mf[mfi]),
//                        BL_TO_FORTRAN_ANYD(v_mf[mfi]),
//                        BL_TO_FORTRAN_ANYD(w_mf[mfi]),
//                        BL_TO_FORTRAN_ANYD(Ip[mfi]),
//                        BL_TO_FORTRAN_ANYD(Im[mfi]),
//                        AMREX_INT_ANYD(domainBox.loVect()),
//                        AMREX_INT_ANYD(domainBox.hiVect()),
//                        bc_f, AMREX_REAL_ANYD(dx), dt, false,
//                        1,1,AMREX_SPACEDIM,false);
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
                // PPM_3d(obx, utilde_mf.array(mfi), 
                //        u_mf.array(mfi), v_mf.array(mfi), w_mf.array(mfi),
                //        Ip.array(mfi), Im.array(mfi), 
                //        domainBox, bcs_u, dx, 
                //        false, 1, 1);
#pragma gpu box(obx)
                ppm_3d(AMREX_INT_ANYD(obx.loVect()),
                       AMREX_INT_ANYD(obx.hiVect()),
                       BL_TO_FORTRAN_ANYD(utilde_mf[mfi]),
                       utilde_mf.nComp(),
                       BL_TO_FORTRAN_ANYD(u_mf[mfi]),
                       BL_TO_FORTRAN_ANYD(v_mf[mfi]),
                       BL_TO_FORTRAN_ANYD(w_mf[mfi]),
                       BL_TO_FORTRAN_ANYD(Ip[mfi]),
                       BL_TO_FORTRAN_ANYD(Im[mfi]),
                       AMREX_INT_ANYD(domainBox.loVect()),
                       AMREX_INT_ANYD(domainBox.hiVect()),
                       bc_f, AMREX_REAL_ANYD(dx), dt, false,
                       2,2,AMREX_SPACEDIM,false);
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
                // PPM_3d(obx, utilde_mf.array(mfi), 
                //        u_mf.array(mfi), v_mf.array(mfi), w_mf.array(mfi),
                //        Ip.array(mfi), Im.array(mfi), 
                //        domainBox, bcs_u, dx, 
                //        false, 2, 2);
#pragma gpu box(obx)
                ppm_3d(AMREX_INT_ANYD(obx.loVect()),
                        AMREX_INT_ANYD(obx.hiVect()),
                        BL_TO_FORTRAN_ANYD(utilde_mf[mfi]),
                        utilde_mf.nComp(),
                        BL_TO_FORTRAN_ANYD(u_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(v_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(w_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(Ip[mfi]),
                        BL_TO_FORTRAN_ANYD(Im[mfi]),
                        AMREX_INT_ANYD(domainBox.loVect()),
                        AMREX_INT_ANYD(domainBox.hiVect()),
                        bc_f, AMREX_REAL_ANYD(dx), dt, false,
                        3,3,AMREX_SPACEDIM,false);
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

void
Maestro::VelPred (Vector<MultiFab>& utilde,
                  const Vector<MultiFab>& ufull,
                  const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& utrans,
                  Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
		          const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& w0mac,
                  const Vector<MultiFab>& force)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::VelPred()",VelPred);

    for (int lev=0; lev<=finest_level; ++lev) {

        // Get the index space and grid spacing of the domain
        const Box& domainBox = geom[lev].Domain();
        const Real* dx = geom[lev].CellSize();

        // get references to the MultiFabs at level lev
              MultiFab& utilde_mf  = utilde[lev];
        const MultiFab& ufull_mf   = ufull[lev];
              MultiFab& umac_mf    = umac[lev][0];
        const MultiFab& utrans_mf  = utrans[lev][0];
        const MultiFab& vtrans_mf  = utrans[lev][1];
              MultiFab& vmac_mf    = umac[lev][1];

        // Ipu and Imu are always the u velocity component
        // Ipv and Imv are always the v velocity component
        MultiFab Ipu, Imu, Ipv, Imv;
        Ipu.define(grids[lev],dmap[lev],AMREX_SPACEDIM,1);
        Imu.define(grids[lev],dmap[lev],AMREX_SPACEDIM,1);
        Ipv.define(grids[lev],dmap[lev],AMREX_SPACEDIM,1);
        Imv.define(grids[lev],dmap[lev],AMREX_SPACEDIM,1);

        MultiFab Ipfx, Imfx, Ipfy, Imfy;
        Ipfx.define(grids[lev],dmap[lev],AMREX_SPACEDIM,1);
        Imfx.define(grids[lev],dmap[lev],AMREX_SPACEDIM,1);
        Ipfy.define(grids[lev],dmap[lev],AMREX_SPACEDIM,1);
        Imfy.define(grids[lev],dmap[lev],AMREX_SPACEDIM,1);

        MultiFab ulx, urx, uimhx, uly, ury, uimhy;
        ulx.define(grids[lev],dmap[lev],AMREX_SPACEDIM,1);
        urx.define(grids[lev],dmap[lev],AMREX_SPACEDIM,1);
        uimhx.define(grids[lev],dmap[lev],AMREX_SPACEDIM,1);
        uly.define(grids[lev],dmap[lev],AMREX_SPACEDIM,1);
        ury.define(grids[lev],dmap[lev],AMREX_SPACEDIM,1);
        uimhy.define(grids[lev],dmap[lev],AMREX_SPACEDIM,1);
#if (AMREX_SPACEDIM == 3)
        const MultiFab& wtrans_mf  = utrans[lev][2];
        MultiFab& wmac_mf          = umac[lev][2];
        const MultiFab& w0macx_mf  = w0mac[lev][0];
        const MultiFab& w0macy_mf  = w0mac[lev][1];
        const MultiFab& w0macz_mf  = w0mac[lev][2];

        MultiFab Ipw, Imw;
        Ipw.define(grids[lev],dmap[lev],AMREX_SPACEDIM,1);
        Imw.define(grids[lev],dmap[lev],AMREX_SPACEDIM,1);

        MultiFab Ipfz, Imfz;
        Ipfz.define(grids[lev],dmap[lev],AMREX_SPACEDIM,1);
        Imfz.define(grids[lev],dmap[lev],AMREX_SPACEDIM,1);

        MultiFab ulz, urz, uimhz;
        ulz.define(grids[lev],dmap[lev],AMREX_SPACEDIM,1);
        urz.define(grids[lev],dmap[lev],AMREX_SPACEDIM,1);
        uimhz.define(grids[lev],dmap[lev],AMREX_SPACEDIM,1);

        MultiFab uimhyz, uimhzy, vimhxz, vimhzx, wimhxy, wimhyx;
        uimhyz.define(grids[lev],dmap[lev],1,1);
        uimhzy.define(grids[lev],dmap[lev],1,1);
        vimhxz.define(grids[lev],dmap[lev],1,1);
        vimhzx.define(grids[lev],dmap[lev],1,1);
        wimhxy.define(grids[lev],dmap[lev],1,1);
        wimhyx.define(grids[lev],dmap[lev],1,1);
#endif
        const MultiFab& force_mf = force[lev];
        const MultiFab& w0_mf = w0_cart[lev];

#if (AMREX_SPACEDIM == 2)
        MultiFab u_mf, v_mf;

        if (ppm_type == 0) {
           u_mf.define(grids[lev],dmap[lev],1,utilde[lev].nGrow());
           v_mf.define(grids[lev],dmap[lev],1,utilde[lev].nGrow());

           MultiFab::Copy(u_mf, utilde[lev], 0, 0, 1, utilde[lev].nGrow());
           MultiFab::Copy(v_mf, utilde[lev], 1, 0, 1, utilde[lev].nGrow());

        } else if (ppm_type == 1 || ppm_type == 2) {

           u_mf.define(grids[lev],dmap[lev],1,ufull[lev].nGrow());
           v_mf.define(grids[lev],dmap[lev],1,ufull[lev].nGrow());

           MultiFab::Copy(u_mf, ufull[lev], 0, 0, 1, ufull[lev].nGrow());
           MultiFab::Copy(v_mf, ufull[lev], 1, 0, 1, ufull[lev].nGrow());
        }

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(utilde_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const Box& obx = amrex::grow(tileBox, 1);
            const Box& xbx = mfi.nodaltilebox(0);
            const Box& ybx = mfi.nodaltilebox(1);
            const Box& mxbx = amrex::growLo(obx,0, -1);
            const Box& mybx = amrex::growLo(obx,1, -1);

            if (ppm_type == 0) {
                // we're going to reuse Ip here as slopex as it has the
                // correct number of ghost zones
                Slopex(obx, utilde_mf.array(mfi), 
                       Ipu.array(mfi), 
                       domainBox, bcs_u, 
                       AMREX_SPACEDIM,0);
            } else {

                // put u on x and y edges, stored in Ipu and Imu
                PPM_2d(obx, utilde_mf.array(mfi), 
                       u_mf.array(mfi), v_mf.array(mfi), 
                       Ipu.array(mfi), Imu.array(mfi), 
                       domainBox, bcs_u, dx, 
                       false, 0, 0);

                if (ppm_trace_forces == 1) {

                    PPM_2d(obx, force_mf.array(mfi), 
                           u_mf.array(mfi), v_mf.array(mfi), 
                           Ipfx.array(mfi), Imfx.array(mfi), 
                           domainBox, bcs_u, dx, 
                           false, 0, 0);
                }
            }

            if (ppm_type == 0) {
               // we're going to reuse Im here as slopey as it has the
               // correct number of ghost zones
               Slopey(obx, utilde_mf.array(mfi), 
                      Imv.array(mfi), 
                      domainBox, bcs_u, 
                      AMREX_SPACEDIM,0);
            } else {

                // put v on x and y edges, stored in Ipv and Imv
                PPM_2d(obx, utilde_mf.array(mfi), 
                       u_mf.array(mfi), v_mf.array(mfi), 
                       Ipv.array(mfi), Imv.array(mfi), 
                       domainBox, bcs_u, dx, 
                       false, 1, 1);

                if (ppm_trace_forces == 1) {

                    PPM_2d(obx, force_mf.array(mfi), 
                           u_mf.array(mfi), v_mf.array(mfi), 
                           Ipv.array(mfi), Imv.array(mfi), 
                           domainBox, bcs_u, dx, 
                           false, 1, 1);
              }
            }

            VelPredInterface(mfi,
                             utilde_mf.array(mfi),
                             ufull_mf.array(mfi),
                             utrans_mf.array(mfi),
                             vtrans_mf.array(mfi),
                             Imu.array(mfi), Ipu.array(mfi),
                             Imv.array(mfi), Ipv.array(mfi),
                             ulx.array(mfi), urx.array(mfi),
                             uimhx.array(mfi),
                             uly.array(mfi), ury.array(mfi),
                             uimhy.array(mfi),
                             domainBox, dx);

            VelPredVelocities(mfi,
                             utilde_mf.array(mfi),
                             utrans_mf.array(mfi),
                             vtrans_mf.array(mfi),
                             umac_mf.array(mfi), vmac_mf.array(mfi),
                             Imfx.array(mfi), Ipfx.array(mfi),
                             Imfy.array(mfi), Ipfy.array(mfi),
                             ulx.array(mfi), urx.array(mfi),
                             uimhx.array(mfi),
                             uly.array(mfi), ury.array(mfi),
                             uimhy.array(mfi),
                             force_mf.array(mfi),
                             w0_mf.array(mfi),
                             domainBox, dx);
        } // end MFIter loop

#elif (AMREX_SPACEDIM == 3)

        MultiFab u_mf, v_mf, w_mf;

        if (ppm_type == 0) {
           u_mf.define(grids[lev],dmap[lev],1,utilde[lev].nGrow());
           v_mf.define(grids[lev],dmap[lev],1,utilde[lev].nGrow());
           w_mf.define(grids[lev],dmap[lev],1,utilde[lev].nGrow());

           MultiFab::Copy(u_mf, utilde[lev], 0, 0, 1, utilde[lev].nGrow());
           MultiFab::Copy(v_mf, utilde[lev], 1, 0, 1, utilde[lev].nGrow());
           MultiFab::Copy(w_mf, utilde[lev], 2, 0, 1, utilde[lev].nGrow());

        } else if (ppm_type == 1 || ppm_type == 2) {

           u_mf.define(grids[lev],dmap[lev],1,ufull[lev].nGrow());
           v_mf.define(grids[lev],dmap[lev],1,ufull[lev].nGrow());
           w_mf.define(grids[lev],dmap[lev],1,ufull[lev].nGrow());

           MultiFab::Copy(u_mf, ufull[lev], 0, 0, 1, ufull[lev].nGrow());
           MultiFab::Copy(v_mf, ufull[lev], 1, 0, 1, ufull[lev].nGrow());
           MultiFab::Copy(w_mf, ufull[lev], 2, 0, 1, ufull[lev].nGrow());
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
            const Box& zbx = mfi.nodaltilebox(2);
            const Box& mxbx = amrex::growLo(obx,0, -1);
            const Box& mybx = amrex::growLo(obx,1, -1);
            const Box& mzbx = amrex::growLo(obx,2, -1);

            // x-direction
            if (ppm_type == 0) {
                // we're going to reuse Ipu here as slopex as it has the
                // correct number of ghost zones
                Slopex(obx, utilde_mf.array(mfi), 
                       Ipu.array(mfi), 
                       domainBox, bcs_u, 
                       AMREX_SPACEDIM,0);

            } else {

                // put u on x, y, and z edges, stored in Ipu and Imu
                PPM_3d(obx, utilde_mf.array(mfi), 
                       u_mf.array(mfi), v_mf.array(mfi), w_mf.array(mfi),
                       Ipu.array(mfi), Imu.array(mfi), 
                       domainBox, bcs_u, dx, 
                       false, 0, 0);

                if (ppm_trace_forces == 1) {
                    PPM_3d(obx, force_mf.array(mfi), 
                            u_mf.array(mfi), v_mf.array(mfi), w_mf.array(mfi),
                            Ipfx.array(mfi), Imfx.array(mfi), 
                            domainBox, bcs_u, dx, 
                            false, 0, 0);
                }
            }

            // y-direction
            if (ppm_type == 0) {
               // we're going to reuse Imv here as slopey as it has the
               // correct number of ghost zones
               Slopey(obx, utilde_mf.array(mfi), 
                       Imv.array(mfi), 
                       domainBox, bcs_u, 
                       AMREX_SPACEDIM,0);

            } else {
                // put v on x, y, and z edges, stored in Ipv and Imv
                PPM_3d(obx, utilde_mf.array(mfi), 
                       u_mf.array(mfi), v_mf.array(mfi), w_mf.array(mfi),
                       Ipv.array(mfi), Imv.array(mfi), 
                       domainBox, bcs_u, dx, 
                       false, 1, 1);

              if (ppm_trace_forces == 1) {

                PPM_3d(obx, force_mf.array(mfi), 
                        u_mf.array(mfi), v_mf.array(mfi), w_mf.array(mfi),
                        Ipfy.array(mfi), Imfy.array(mfi), 
                        domainBox, bcs_u, dx, 
                        false, 1, 1);
              }
            }

            // z-direction
            if (ppm_type == 0) {
               // we're going to reuse Imw here as slopey as it has the
               // correct number of ghost zones

               Slopez(obx, utilde_mf.array(mfi), 
                      Imw.array(mfi), 
                      domainBox, bcs_u, 
                      AMREX_SPACEDIM,0);

            } else {
                // put w on x, y, and z edges, stored in Ipw and Imw
                PPM_3d(obx, utilde_mf.array(mfi), 
                       u_mf.array(mfi), v_mf.array(mfi), w_mf.array(mfi),
                       Ipw.array(mfi), Imw.array(mfi), 
                       domainBox, bcs_u, dx, 
                       false, 2, 2);

              if (ppm_trace_forces == 1) {

                PPM_3d(obx, force_mf.array(mfi), 
                        u_mf.array(mfi), v_mf.array(mfi), w_mf.array(mfi),
                        Ipfz.array(mfi), Imfz.array(mfi), 
                        domainBox, bcs_u, dx, 
                        false, 2, 2);
              }
            }
            
            VelPredInterface(mfi,
                             utilde_mf.array(mfi),
                             ufull_mf.array(mfi),
                             utrans_mf.array(mfi),
                             vtrans_mf.array(mfi),
                             wtrans_mf.array(mfi),
                             Imu.array(mfi), Ipu.array(mfi),
                             Imv.array(mfi), Ipv.array(mfi),
                             Imw.array(mfi), Ipw.array(mfi),
                             ulx.array(mfi), urx.array(mfi),
                             uimhx.array(mfi),
                             uly.array(mfi), ury.array(mfi),
                             uimhy.array(mfi),
                             ulz.array(mfi), urz.array(mfi),
                             uimhz.array(mfi),
                             domainBox, dx);

            VelPredTransverse(mfi,
                            utilde_mf.array(mfi),
                            utrans_mf.array(mfi),
                            vtrans_mf.array(mfi),
                            wtrans_mf.array(mfi),
                            ulx.array(mfi), urx.array(mfi),
                            uimhx.array(mfi),
                            uly.array(mfi), ury.array(mfi),
                            uimhy.array(mfi),
                            ulz.array(mfi), urz.array(mfi),
                            uimhz.array(mfi),
                            uimhyz.array(mfi), uimhzy.array(mfi), 
                            vimhxz.array(mfi), vimhzx.array(mfi), 
                            wimhxy.array(mfi), wimhyx.array(mfi), 
                            domainBox, dx);

            VelPredVelocities(mfi,
                            utilde_mf.array(mfi),
                            utrans_mf.array(mfi),
                            vtrans_mf.array(mfi),
                            wtrans_mf.array(mfi),
                            umac_mf.array(mfi), vmac_mf.array(mfi),
                            wmac_mf.array(mfi),
                            w0macx_mf.array(mfi), 
                            w0macy_mf.array(mfi),
                            w0macz_mf.array(mfi),
                            Imfx.array(mfi), Ipfx.array(mfi),
                            Imfy.array(mfi), Ipfy.array(mfi),
                            Imfz.array(mfi), Ipfz.array(mfi),
                            ulx.array(mfi), urx.array(mfi),
                            uly.array(mfi), ury.array(mfi),
                            ulz.array(mfi), urz.array(mfi),
                            uimhyz.array(mfi), uimhzy.array(mfi), 
                            vimhxz.array(mfi), vimhzx.array(mfi), 
                            wimhxy.array(mfi), wimhyx.array(mfi), 
                            force_mf.array(mfi),
                            w0_mf.array(mfi),
                            domainBox, dx);
        } // end MFIter loop

#endif // AMREX_SPACEDIM
    } // end loop over levels

    // edge_restriction
    AverageDownFaces(umac);
}

void
Maestro::MakeEdgeScal (Vector<MultiFab>& state,
                       Vector<std::array< MultiFab, AMREX_SPACEDIM > >& sedge,
                       Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
                       Vector<MultiFab>& force,
                       int is_vel, const Vector<BCRec>& bcs, int nbccomp,
                       int start_scomp, int start_bccomp, int num_comp, int is_conservative)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeEdgeScal()", MakeEdgeScal);

    for (int lev=0; lev<=finest_level; ++lev) {

        // Get the index space and grid spacing of the domain
        const Box& domainBox = geom[lev].Domain();
        const Real* dx = geom[lev].CellSize();

        // get references to the MultiFabs at level lev
        const MultiFab& scal_mf   = state[lev];

        MultiFab Ip, Im, Ipf, Imf;
        Ip.define(grids[lev],dmap[lev],AMREX_SPACEDIM,1);
        Im.define(grids[lev],dmap[lev],AMREX_SPACEDIM,1);
        Ipf.define(grids[lev],dmap[lev],AMREX_SPACEDIM,1);
        Imf.define(grids[lev],dmap[lev],AMREX_SPACEDIM,1);

        MultiFab slx, srx, simhx;
        slx.define(grids[lev],dmap[lev],1,1);
        srx.define(grids[lev],dmap[lev],1,1);
        simhx.define(grids[lev],dmap[lev],1,1);

        MultiFab sly, sry, simhy;
        sly.define(grids[lev],dmap[lev],1,1);
        sry.define(grids[lev],dmap[lev],1,1);
        simhy.define(grids[lev],dmap[lev],1,1);

        slx.setVal(0.);
        srx.setVal(0.);
        simhx.setVal(0.);
        sly.setVal(0.);
        sry.setVal(0.);
        simhy.setVal(0.);

#if (AMREX_SPACEDIM == 3)

        MultiFab slopez, divu;
        slopez.define(grids[lev],dmap[lev],1,1);
        divu.define(grids[lev],dmap[lev],1,1);

        MultiFab slz, srz, simhz;
        slz.define(grids[lev],dmap[lev],1,1);
        srz.define(grids[lev],dmap[lev],1,1);
        simhz.define(grids[lev],dmap[lev],1,1);

        MultiFab simhxy, simhxz, simhyx, simhyz, simhzx, simhzy;
        simhxy.define(grids[lev],dmap[lev],1,1);
        simhxz.define(grids[lev],dmap[lev],1,1);
        simhyx.define(grids[lev],dmap[lev],1,1);
        simhyz.define(grids[lev],dmap[lev],1,1);
        simhzx.define(grids[lev],dmap[lev],1,1);
        simhzy.define(grids[lev],dmap[lev],1,1);

        slx.setVal(0.);
        srx.setVal(0.);
        simhx.setVal(0.);
        sly.setVal(0.);
        sry.setVal(0.);
        simhy.setVal(0.);
        slz.setVal(0.);
        srz.setVal(0.);
        simhz.setVal(0.);

        simhxy.setVal(0.);
        simhxz.setVal(0.);
        simhyx.setVal(0.);
        simhyz.setVal(0.);
        simhzx.setVal(0.);
        simhzy.setVal(0.);
#endif
        const MultiFab& force_mf = force[lev];

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#if (AMREX_SPACEDIM == 2)

        Vector<MultiFab> vec_scal_mf(num_comp);
        for (int comp=0; comp < num_comp; ++comp) {
            vec_scal_mf[comp].define(grids[lev],dmap[lev],1,scal_mf.nGrow());

            MultiFab::Copy(vec_scal_mf[comp], scal_mf, start_scomp+comp, 0, 1, scal_mf.nGrow());
        }
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(scal_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const Box& obx = amrex::grow(tileBox, 1);

            // Be careful to pass in comp+1 for fortran indexing
            for (int scomp = start_scomp; scomp < start_scomp + num_comp; ++scomp) {

                int vcomp = scomp - start_scomp;
                int bccomp = start_bccomp + scomp - start_scomp;

                Array4<Real> const scal_arr = state[lev].array(mfi);

                Array4<Real> const umac_arr = umac[lev][0].array(mfi);
                Array4<Real> const vmac_arr = umac[lev][1].array(mfi);

                Array4<Real> const slx_arr = slx.array(mfi);
                Array4<Real> const srx_arr = srx.array(mfi);
                Array4<Real> const sly_arr = sly.array(mfi);
                Array4<Real> const sry_arr = sry.array(mfi);

                Array4<Real> const simhx_arr = simhx.array(mfi);
                Array4<Real> const simhy_arr = simhy.array(mfi);

                if (ppm_type == 0) {
                    // we're going to reuse Ip here as slopex and Im as slopey
                    // as they have the correct number of ghost zones

                    // x-direction
                    Slopex(obx, vec_scal_mf[vcomp].array(mfi), 
                           Ip.array(mfi), 
                           domainBox, bcs, 
                           1,bccomp);

                    // y-direction
                    Slopey(obx, vec_scal_mf[vcomp].array(mfi), 
                           Im.array(mfi), 
                           domainBox, bcs, 
                           1,bccomp);

                } else {

                    PPM_2d(obx, scal_arr, 
                           umac_arr, vmac_arr, 
                           Ip.array(mfi), Im.array(mfi), 
                           domainBox, bcs, dx, 
                           true, scomp, bccomp);

                    if (ppm_trace_forces == 1) {

                        PPM_2d(obx, force[lev].array(mfi), 
                           umac_arr, vmac_arr, 
                           Ipf.array(mfi), Imf.array(mfi), 
                           domainBox, bcs, dx, 
                           true, scomp, bccomp);
                    }
                }

                // Create s_{\i-\half\e_x}^x, etc.

                MakeEdgeScalPredictor(mfi, slx_arr, srx_arr,
                                      sly_arr, sry_arr,
                                      scal_arr, 
                                      Ip.array(mfi), Im.array(mfi), 
                                      umac_arr, vmac_arr, 
                                      simhx_arr, simhy_arr, 
                                      domainBox, bcs, dx,
                                      scomp, bccomp, is_vel);

                Array4<Real> const sedgex_arr = sedge[lev][0].array(mfi);
                Array4<Real> const sedgey_arr = sedge[lev][1].array(mfi);

                // Create sedgelx, etc.

                MakeEdgeScalEdges(mfi, slx_arr, srx_arr,
                                  sly_arr, sry_arr,
                                  scal_arr, 
                                  sedgex_arr, sedgey_arr, 
                                  force[lev].array(mfi),
                                  umac_arr, vmac_arr, 
                                  Ipf.array(mfi), Imf.array(mfi),
                                  simhx_arr, simhy_arr, 
                                  domainBox, bcs, dx,
                                  scomp, bccomp, 
                                  is_vel, is_conservative);
            } // end loop over components
        } // end MFIter loop

#elif (AMREX_SPACEDIM == 3)

        Vector<MultiFab> vec_scal_mf(num_comp);
        for (int comp=0; comp < num_comp; ++comp) {
            vec_scal_mf[comp].define(grids[lev],dmap[lev],1,scal_mf.nGrow());
            vec_scal_mf[comp].setVal(0.);

            MultiFab::Copy(vec_scal_mf[comp], scal_mf, start_scomp+comp, 0, 1, scal_mf.nGrow());
        }

        // Be careful to pass in comp+1 for fortran indexing
        for (int scomp = start_scomp; scomp < start_scomp + num_comp; ++scomp) {

            int vcomp = scomp - start_scomp;
            int bccomp = start_bccomp + scomp - start_scomp;

#ifdef _OPENMP
#pragma omp parallel
#endif
            for ( MFIter mfi(scal_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

                // Get the index space of the valid region
                const Box& tileBox = mfi.tilebox();
                const Box& obx = amrex::grow(tileBox, 1);
                const Box& xbx = mfi.nodaltilebox(0);
                const Box& ybx = mfi.nodaltilebox(1);
                const Box& zbx = mfi.nodaltilebox(2);
                const Box& mxbx = amrex::growLo(obx, 0, -1);
                const Box& mybx = amrex::growLo(obx, 1, -1);
                const Box& mzbx = amrex::growLo(obx, 2, -1);

                Array4<Real> const umac_arr = umac[lev][0].array(mfi);
                Array4<Real> const vmac_arr = umac[lev][1].array(mfi);
                Array4<Real> const wmac_arr = umac[lev][2].array(mfi);

                // make divu 
                if (is_conservative) {
                    MakeDivU(obx, divu.array(mfi), 
                             umac_arr, vmac_arr, wmac_arr, dx);
                }
                          
                if (ppm_type == 0) {
                    // we're going to reuse Ip here as slopex and Im as slopey
                    // as they have the correct number of ghost zones

                    // x-direction
                    Slopex(obx, vec_scal_mf[vcomp].array(mfi), 
                           Ip.array(mfi), 
                           domainBox, bcs, 
                           1,bccomp);

                    // y-direction
                    Slopey(obx, vec_scal_mf[vcomp].array(mfi), 
                           Im.array(mfi), 
                           domainBox, bcs, 
                           1,bccomp);

                    // z-direction
                    Slopez(obx, vec_scal_mf[vcomp].array(mfi), 
                           slopez.array(mfi), 
                           domainBox, bcs, 
                           1,bccomp);

                } else {

                    PPM_3d(obx, state[lev].array(mfi), 
                           umac_arr, vmac_arr, wmac_arr,
                           Ip.array(mfi), Im.array(mfi), 
                           domainBox, bcs, dx, 
                           true, scomp, bccomp);

                    if (ppm_trace_forces == 1) {

                        PPM_3d(obx, force[lev].array(mfi), 
                           umac_arr, vmac_arr, wmac_arr,
                           Ipf.array(mfi), Imf.array(mfi), 
                           domainBox, bcs, dx, 
                           true, scomp, bccomp);
                    }
                }
            }

#ifdef _OPENMP
#pragma omp parallel
#endif
            for ( MFIter mfi(scal_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

                Array4<Real> const scal_arr = state[lev].array(mfi);

                Array4<Real> const umac_arr = umac[lev][0].array(mfi);
                Array4<Real> const vmac_arr = umac[lev][1].array(mfi);
                Array4<Real> const wmac_arr = umac[lev][2].array(mfi);

                Array4<Real> const slx_arr = slx.array(mfi);
                Array4<Real> const srx_arr = srx.array(mfi);
                Array4<Real> const sly_arr = sly.array(mfi);
                Array4<Real> const sry_arr = sry.array(mfi);
                Array4<Real> const slz_arr = slz.array(mfi);
                Array4<Real> const srz_arr = srz.array(mfi);

                Array4<Real> const simhx_arr = simhx.array(mfi);
                Array4<Real> const simhy_arr = simhy.array(mfi);
                Array4<Real> const simhz_arr = simhz.array(mfi);

                // Create s_{\i-\half\e_x}^x, etc.

                MakeEdgeScalPredictor(mfi, slx_arr, srx_arr,
                                      sly_arr, sry_arr,
                                      slz_arr, srz_arr,
                                      scal_arr, 
                                      Ip.array(mfi), Im.array(mfi), 
                                      slopez.array(mfi),
                                      umac_arr, vmac_arr, wmac_arr,
                                      simhx_arr, simhy_arr, simhz_arr,
                                      domainBox, bcs, dx,
                                      scomp, bccomp, is_vel);

                Array4<Real> const simhxy_arr = simhxy.array(mfi);
                Array4<Real> const simhxz_arr = simhxz.array(mfi);
                Array4<Real> const simhyx_arr = simhyx.array(mfi);
                Array4<Real> const simhyz_arr = simhyz.array(mfi);
                Array4<Real> const simhzx_arr = simhzx.array(mfi);
                Array4<Real> const simhzy_arr = simhzy.array(mfi);

                // Create transverse terms, s_{\i-\half\e_x}^{x|y}, etc.

                MakeEdgeScalTransverse(mfi, slx_arr, srx_arr,
                                       sly_arr, sry_arr,
                                       slz_arr, srz_arr,
                                       scal_arr, divu.array(mfi),
                                       umac_arr, vmac_arr, wmac_arr,
                                       simhx_arr, simhy_arr, simhz_arr,
                                       simhxy_arr, simhxz_arr, simhyx_arr,
                                       simhyz_arr, simhzx_arr, simhzy_arr,
                                       domainBox, bcs, dx,
                                       scomp, bccomp, 
                                       is_vel, is_conservative);

                Array4<Real> const sedgex_arr = sedge[lev][0].array(mfi);
                Array4<Real> const sedgey_arr = sedge[lev][1].array(mfi);
                Array4<Real> const sedgez_arr = sedge[lev][2].array(mfi);

                // Create sedgelx, etc.

                MakeEdgeScalEdges(mfi, slx_arr, srx_arr,
                                  sly_arr, sry_arr,
                                  slz_arr, srz_arr, scal_arr, 
                                  sedgex_arr, sedgey_arr, sedgez_arr,
                                  force[lev].array(mfi),
                                  umac_arr, vmac_arr, wmac_arr,
                                  Ipf.array(mfi), Imf.array(mfi),
                                  simhxy_arr, simhxz_arr, simhyx_arr,
                                  simhyz_arr, simhzx_arr, simhzy_arr,
                                  domainBox, bcs, dx,
                                  scomp, bccomp, 
                                  is_vel, is_conservative);
            } // end MFIter loop
        } // end loop over components
#endif
    } // end loop over levels

    // We use edge_restriction for the output velocity if is_vel == 1
    // we do not use edge_restriction for scalars because instead we will use
    // reflux on the fluxes in make_flux.
    if (is_vel == 1) {
        if (reflux_type == 1 || reflux_type == 2) {
            AverageDownFaces(sedge);
        }
    }
}

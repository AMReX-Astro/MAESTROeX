
#include <Maestro.H>
#include <MaestroHydro_F.H>
#include <MaestroBCThreads.H>

using namespace amrex;

void
Maestro::MakeUtrans (const Vector<MultiFab>& utilde,
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
        const MultiFab& utilde_mf  = utilde[lev];
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
#pragma gpu box(obx)
                slopex_2d(AMREX_INT_ANYD(obx.loVect()),
                       AMREX_INT_ANYD(obx.hiVect()),
                       BL_TO_FORTRAN_ANYD(u_mf[mfi]),
                       u_mf.nComp(),
                       BL_TO_FORTRAN_ANYD(Ip[mfi]),Ip.nComp(),
                       AMREX_INT_ANYD(domainBox.loVect()),
                       AMREX_INT_ANYD(domainBox.hiVect()),
                       1,bc_f,AMREX_SPACEDIM,1);

            } else {

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
                       1,1,AMREX_SPACEDIM);
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
#pragma gpu box(obx)
                slopey_2d(AMREX_INT_ANYD(obx.loVect()),
                       AMREX_INT_ANYD(obx.hiVect()),
                       BL_TO_FORTRAN_ANYD(v_mf[mfi]),
                       v_mf.nComp(),
                       BL_TO_FORTRAN_ANYD(Im[mfi]),Im.nComp(),
                       AMREX_INT_ANYD(domainBox.loVect()),
                       AMREX_INT_ANYD(domainBox.hiVect()),
                       1,bc_f,AMREX_SPACEDIM,2);

            } else {

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

#pragma gpu box(obx)
                slopex_2d(AMREX_INT_ANYD(obx.loVect()),
                       AMREX_INT_ANYD(obx.hiVect()),
                       BL_TO_FORTRAN_ANYD(u_mf[mfi]),
                       u_mf.nComp(),
                       BL_TO_FORTRAN_ANYD(Ip[mfi]),Ip.nComp(),
                       AMREX_INT_ANYD(domainBox.loVect()),
                       AMREX_INT_ANYD(domainBox.hiVect()),
                       1,bc_f,AMREX_SPACEDIM,1);
            } else {
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
                       BL_TO_FORTRAN_ANYD(u_mf[mfi]),
                       BL_TO_FORTRAN_ANYD(v_mf[mfi]),
                       AMREX_INT_ANYD(domainBox.loVect()),
                       AMREX_INT_ANYD(domainBox.hiVect()),
                       bc_f, AMREX_REAL_ANYD(dx), dt, false,
                       1,1,AMREX_SPACEDIM,false);
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

#pragma gpu box(obx)
                slopey_2d(AMREX_INT_ANYD(obx.loVect()),
                       AMREX_INT_ANYD(obx.hiVect()),
                       BL_TO_FORTRAN_ANYD(v_mf[mfi]),
                       v_mf.nComp(),
                       BL_TO_FORTRAN_ANYD(Im[mfi]),Im.nComp(),
                       AMREX_INT_ANYD(domainBox.loVect()),
                       AMREX_INT_ANYD(domainBox.hiVect()),
                       1,bc_f,AMREX_SPACEDIM,2);
            } else {
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
                       BL_TO_FORTRAN_ANYD(u_mf[mfi]),
                       BL_TO_FORTRAN_ANYD(v_mf[mfi]),
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

#pragma gpu box(obx)
                slopez_3d(AMREX_INT_ANYD(obx.loVect()),
                       AMREX_INT_ANYD(obx.hiVect()),
                       BL_TO_FORTRAN_ANYD(w_mf[mfi]),
                       w_mf.nComp(),
                       BL_TO_FORTRAN_ANYD(Im[mfi]),Im.nComp(),
                       AMREX_INT_ANYD(domainBox.loVect()),
                       AMREX_INT_ANYD(domainBox.hiVect()),
                       1,bc_f,AMREX_SPACEDIM,3);
            } else {
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
                        BL_TO_FORTRAN_ANYD(u_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(v_mf[mfi]),
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
Maestro::VelPred (const Vector<MultiFab>& utilde,
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
        const MultiFab& utilde_mf  = utilde[lev];
        const MultiFab& ufull_mf   = ufull[lev];
              MultiFab& umac_mf    = umac[lev][0];
        const MultiFab& utrans_mf  = utrans[lev][0];
        const MultiFab& vtrans_mf  = utrans[lev][1];
              MultiFab& vmac_mf    = umac[lev][1];
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
        MultiFab& wmac_mf    = umac[lev][2];
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

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)

#if (AMREX_SPACEDIM == 2)

#ifdef AMREX_USE_CUDA
        int* bc_f = prepare_bc(bcs_u[0].data(), 1);
#else
        const int* bc_f = bcs_u[0].data();
#endif
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
#pragma gpu box(obx)
                slopex_2d(AMREX_INT_ANYD(obx.loVect()),
                       AMREX_INT_ANYD(obx.hiVect()),
                       BL_TO_FORTRAN_ANYD(utilde_mf[mfi]),
                       utilde_mf.nComp(),
                       BL_TO_FORTRAN_ANYD(Ipu[mfi]),Ipu.nComp(),
                       AMREX_INT_ANYD(domainBox.loVect()),
                       AMREX_INT_ANYD(domainBox.hiVect()),
                       2,bc_f,AMREX_SPACEDIM,1);

            } else {

#pragma gpu box(obx)
                ppm_2d(AMREX_INT_ANYD(obx.loVect()),
                       AMREX_INT_ANYD(obx.hiVect()),
                       BL_TO_FORTRAN_ANYD(utilde_mf[mfi]),
                       utilde_mf.nComp(),
                       BL_TO_FORTRAN_ANYD(u_mf[mfi]),
                       BL_TO_FORTRAN_ANYD(v_mf[mfi]),
                       BL_TO_FORTRAN_ANYD(Ipu[mfi]),
                       BL_TO_FORTRAN_ANYD(Imu[mfi]),
                       AMREX_INT_ANYD(domainBox.loVect()),
                       AMREX_INT_ANYD(domainBox.hiVect()),
                       bc_f, AMREX_REAL_ANYD(dx), dt, false,
                       1,1,AMREX_SPACEDIM);

                if (ppm_trace_forces == 1) {
#pragma gpu box(obx)
                    ppm_2d(AMREX_INT_ANYD(obx.loVect()),
                           AMREX_INT_ANYD(obx.hiVect()),
                           BL_TO_FORTRAN_ANYD(force_mf[mfi]),
                           force_mf.nComp(),
                           BL_TO_FORTRAN_ANYD(u_mf[mfi]),
                           BL_TO_FORTRAN_ANYD(v_mf[mfi]),
                           BL_TO_FORTRAN_ANYD(Ipfx[mfi]),
                           BL_TO_FORTRAN_ANYD(Imfx[mfi]),
                           AMREX_INT_ANYD(domainBox.loVect()),
                           AMREX_INT_ANYD(domainBox.hiVect()),
                           bc_f, AMREX_REAL_ANYD(dx), dt, false,
                           1,1,AMREX_SPACEDIM);
                }
            }

            if (ppm_type == 0) {
               // we're going to reuse Im here as slopey as it has the
               // correct number of ghost zones
#pragma gpu box(obx)
               slopey_2d(AMREX_INT_ANYD(obx.loVect()),
                      AMREX_INT_ANYD(obx.hiVect()),
                      BL_TO_FORTRAN_ANYD(utilde_mf[mfi]),
                      utilde_mf.nComp(),
                      BL_TO_FORTRAN_ANYD(Imv[mfi]),Imv.nComp(),
                      AMREX_INT_ANYD(domainBox.loVect()),
                      AMREX_INT_ANYD(domainBox.hiVect()),
                      2,bc_f,AMREX_SPACEDIM,1);

            } else {

#pragma gpu box(obx)
               ppm_2d(AMREX_INT_ANYD(obx.loVect()),
                      AMREX_INT_ANYD(obx.hiVect()),
                      BL_TO_FORTRAN_ANYD(utilde_mf[mfi]),
                      utilde_mf.nComp(),
                      BL_TO_FORTRAN_ANYD(u_mf[mfi]),
                      BL_TO_FORTRAN_ANYD(v_mf[mfi]),
                      BL_TO_FORTRAN_ANYD(Ipv[mfi]),
                      BL_TO_FORTRAN_ANYD(Imv[mfi]),
                      AMREX_INT_ANYD(domainBox.loVect()),
                      AMREX_INT_ANYD(domainBox.hiVect()),
                      bc_f, AMREX_REAL_ANYD(dx), dt, false,
                      2,2,AMREX_SPACEDIM);

              if (ppm_trace_forces == 1) {
#pragma gpu box(obx)
                ppm_2d(AMREX_INT_ANYD(obx.loVect()),
                    AMREX_INT_ANYD(obx.hiVect()),
                    BL_TO_FORTRAN_ANYD(force_mf[mfi]),
                    force_mf.nComp(),
                    BL_TO_FORTRAN_ANYD(u_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(v_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(Ipfy[mfi]),
                    BL_TO_FORTRAN_ANYD(Imfy[mfi]),
                    AMREX_INT_ANYD(domainBox.loVect()),
                    AMREX_INT_ANYD(domainBox.hiVect()),
                    bc_f, AMREX_REAL_ANYD(dx), dt, false,
                    2,2,AMREX_SPACEDIM);
              }
            }

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.

            // x-direction
#pragma gpu box(mxbx)
            velpred_interface_2d(AMREX_INT_ANYD(mxbx.loVect()), AMREX_INT_ANYD(mxbx.hiVect()),1,
                        AMREX_INT_ANYD(domainBox.loVect()), AMREX_INT_ANYD(domainBox.hiVect()),
                        BL_TO_FORTRAN_ANYD(utilde_mf[mfi]), utilde_mf.nComp(), utilde_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(ufull_mf[mfi]), ufull_mf.nComp(), ufull_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(utrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(Imu[mfi]),
                        BL_TO_FORTRAN_ANYD(Ipu[mfi]),
                        BL_TO_FORTRAN_ANYD(Imv[mfi]),
                        BL_TO_FORTRAN_ANYD(Ipv[mfi]),
                        BL_TO_FORTRAN_ANYD(ulx[mfi]),
                        BL_TO_FORTRAN_ANYD(urx[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhx[mfi]),
                        AMREX_REAL_ANYD(dx), dt, bc_f, phys_bc.dataPtr());

            // y-direction
#pragma gpu box(mybx)
            velpred_interface_2d(AMREX_INT_ANYD(mybx.loVect()), AMREX_INT_ANYD(mybx.hiVect()),2,
                        AMREX_INT_ANYD(domainBox.loVect()), AMREX_INT_ANYD(domainBox.hiVect()),
                        BL_TO_FORTRAN_ANYD(utilde_mf[mfi]), utilde_mf.nComp(), utilde_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(ufull_mf[mfi]), ufull_mf.nComp(), ufull_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(vtrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(Imu[mfi]),
                        BL_TO_FORTRAN_ANYD(Ipu[mfi]),
                        BL_TO_FORTRAN_ANYD(Imv[mfi]),
                        BL_TO_FORTRAN_ANYD(Ipv[mfi]),
                        BL_TO_FORTRAN_ANYD(uly[mfi]),
                        BL_TO_FORTRAN_ANYD(ury[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhy[mfi]),
                        AMREX_REAL_ANYD(dx), dt, bc_f, phys_bc.dataPtr());

            // x-direction
#pragma gpu box(xbx)
            velpred_2d(AMREX_INT_ANYD(xbx.loVect()), AMREX_INT_ANYD(xbx.hiVect()),lev,1,
                        AMREX_INT_ANYD(domainBox.loVect()), AMREX_INT_ANYD(domainBox.hiVect()),
                        BL_TO_FORTRAN_ANYD(utilde_mf[mfi]), utilde_mf.nComp(), utilde_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(utrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(vtrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(Imfx[mfi]),
                        BL_TO_FORTRAN_ANYD(Ipfx[mfi]),
                        BL_TO_FORTRAN_ANYD(ulx[mfi]),
                        BL_TO_FORTRAN_ANYD(urx[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhx[mfi]),
                        BL_TO_FORTRAN_ANYD(uly[mfi]),
                        BL_TO_FORTRAN_ANYD(ury[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhy[mfi]),
                        BL_TO_FORTRAN_ANYD(force_mf[mfi]), force_mf.nComp(), force_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(w0_mf[mfi]), 
                        AMREX_REAL_ANYD(dx), dt, bc_f, phys_bc.dataPtr());

            // y-direction
#pragma gpu box(ybx)
            velpred_2d(AMREX_INT_ANYD(ybx.loVect()), AMREX_INT_ANYD(ybx.hiVect()),lev,2,
                        AMREX_INT_ANYD(domainBox.loVect()), AMREX_INT_ANYD(domainBox.hiVect()),
                        BL_TO_FORTRAN_ANYD(utilde_mf[mfi]), utilde_mf.nComp(), utilde_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(utrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(vtrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(Imfy[mfi]),
                        BL_TO_FORTRAN_ANYD(Ipfy[mfi]),
                        BL_TO_FORTRAN_ANYD(ulx[mfi]),
                        BL_TO_FORTRAN_ANYD(urx[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhx[mfi]),
                        BL_TO_FORTRAN_ANYD(uly[mfi]),
                        BL_TO_FORTRAN_ANYD(ury[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhy[mfi]),
                        BL_TO_FORTRAN_ANYD(force_mf[mfi]), force_mf.nComp(), force_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(w0_mf[mfi]), 
                        AMREX_REAL_ANYD(dx), dt, bc_f, phys_bc.dataPtr());
        } // end MFIter loop

#ifdef AMREX_USE_CUDA
        clean_bc(bc_f);
#endif

#elif (AMREX_SPACEDIM == 3)

#ifdef AMREX_USE_CUDA
        int* bc_f = prepare_bc(bcs_u[0].data(), 1);
#else
        const int* bc_f = bcs_u[0].data();
#endif
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
#pragma gpu box(obx)
                slopex_2d(AMREX_INT_ANYD(obx.loVect()),
                       AMREX_INT_ANYD(obx.hiVect()),
                       BL_TO_FORTRAN_ANYD(utilde_mf[mfi]),
                       utilde_mf.nComp(),
                       BL_TO_FORTRAN_ANYD(Ipu[mfi]),Ipu.nComp(),
                       AMREX_INT_ANYD(domainBox.loVect()),
                       AMREX_INT_ANYD(domainBox.hiVect()),
                       3,bc_f,AMREX_SPACEDIM,1);

            } else {

#pragma gpu box(obx)
                ppm_3d(AMREX_INT_ANYD(obx.loVect()),
                       AMREX_INT_ANYD(obx.hiVect()),
                       BL_TO_FORTRAN_ANYD(utilde_mf[mfi]),
                       utilde_mf.nComp(),
                       BL_TO_FORTRAN_ANYD(u_mf[mfi]),
                       BL_TO_FORTRAN_ANYD(v_mf[mfi]),
                       BL_TO_FORTRAN_ANYD(w_mf[mfi]),
                       BL_TO_FORTRAN_ANYD(Ipu[mfi]),
                       BL_TO_FORTRAN_ANYD(Imu[mfi]),
                       BL_TO_FORTRAN_ANYD(u_mf[mfi]),
                       BL_TO_FORTRAN_ANYD(v_mf[mfi]),
                       AMREX_INT_ANYD(domainBox.loVect()),
                       AMREX_INT_ANYD(domainBox.hiVect()),
                       bc_f, AMREX_REAL_ANYD(dx), dt, false,
                       1,1,AMREX_SPACEDIM,false);

                if (ppm_trace_forces == 1) {
#pragma gpu box(obx)
                    ppm_3d(AMREX_INT_ANYD(obx.loVect()),
                           AMREX_INT_ANYD(obx.hiVect()),
                           BL_TO_FORTRAN_ANYD(force_mf[mfi]),
                           force_mf.nComp(),
                           BL_TO_FORTRAN_ANYD(u_mf[mfi]),
                           BL_TO_FORTRAN_ANYD(v_mf[mfi]),
                           BL_TO_FORTRAN_ANYD(w_mf[mfi]),
                           BL_TO_FORTRAN_ANYD(Ipfx[mfi]),
                           BL_TO_FORTRAN_ANYD(Imfx[mfi]),
                           BL_TO_FORTRAN_ANYD(u_mf[mfi]),
                           BL_TO_FORTRAN_ANYD(v_mf[mfi]),
                           AMREX_INT_ANYD(domainBox.loVect()),
                           AMREX_INT_ANYD(domainBox.hiVect()),
                           bc_f, AMREX_REAL_ANYD(dx), dt, false,
                           1,1,AMREX_SPACEDIM,false);
                }
            }

            // y-direction
            if (ppm_type == 0) {
               // we're going to reuse Imv here as slopey as it has the
               // correct number of ghost zones
#pragma gpu box(obx)
               slopey_2d(AMREX_INT_ANYD(obx.loVect()),
                      AMREX_INT_ANYD(obx.hiVect()),
                      BL_TO_FORTRAN_ANYD(utilde_mf[mfi]),
                      utilde_mf.nComp(),
                      BL_TO_FORTRAN_ANYD(Imv[mfi]),Imv.nComp(),
                      AMREX_INT_ANYD(domainBox.loVect()),
                      AMREX_INT_ANYD(domainBox.hiVect()),
                      3,bc_f,AMREX_SPACEDIM,1);

            } else {

#pragma gpu box(obx)
               ppm_3d(AMREX_INT_ANYD(obx.loVect()),
                      AMREX_INT_ANYD(obx.hiVect()),
                      BL_TO_FORTRAN_ANYD(utilde_mf[mfi]),
                      utilde_mf.nComp(),
                      BL_TO_FORTRAN_ANYD(u_mf[mfi]),
                      BL_TO_FORTRAN_ANYD(v_mf[mfi]),
                      BL_TO_FORTRAN_ANYD(w_mf[mfi]),
                      BL_TO_FORTRAN_ANYD(Ipv[mfi]),
                      BL_TO_FORTRAN_ANYD(Imv[mfi]),
                      BL_TO_FORTRAN_ANYD(u_mf[mfi]),
                      BL_TO_FORTRAN_ANYD(v_mf[mfi]),
                      AMREX_INT_ANYD(domainBox.loVect()),
                      AMREX_INT_ANYD(domainBox.hiVect()),
                      bc_f, AMREX_REAL_ANYD(dx), dt, false,
                      2,2,AMREX_SPACEDIM,false);

              if (ppm_trace_forces == 1) {
#pragma gpu box(obx)
                ppm_3d(AMREX_INT_ANYD(obx.loVect()),
                    AMREX_INT_ANYD(obx.hiVect()),
                    BL_TO_FORTRAN_ANYD(force_mf[mfi]),
                    force_mf.nComp(),
                    BL_TO_FORTRAN_ANYD(u_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(v_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(w_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(Ipfy[mfi]),
                    BL_TO_FORTRAN_ANYD(Imfy[mfi]),
                    BL_TO_FORTRAN_ANYD(u_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(v_mf[mfi]),
                    AMREX_INT_ANYD(domainBox.loVect()),
                    AMREX_INT_ANYD(domainBox.hiVect()),
                    bc_f, AMREX_REAL_ANYD(dx), dt, false,
                    2,2,AMREX_SPACEDIM,false);
              }
            }

            // z-direction
            if (ppm_type == 0) {
               // we're going to reuse Imw here as slopey as it has the
               // correct number of ghost zones
#pragma gpu box(obx)
               slopez_3d(AMREX_INT_ANYD(obx.loVect()),
                      AMREX_INT_ANYD(obx.hiVect()),
                      BL_TO_FORTRAN_ANYD(utilde_mf[mfi]),
                      utilde_mf.nComp(),
                      BL_TO_FORTRAN_ANYD(Imw[mfi]),Imw.nComp(),
                      AMREX_INT_ANYD(domainBox.loVect()),
                      AMREX_INT_ANYD(domainBox.hiVect()),
                      3,bc_f,AMREX_SPACEDIM,1);

            } else {

#pragma gpu box(obx)
               ppm_3d(AMREX_INT_ANYD(obx.loVect()),
                      AMREX_INT_ANYD(obx.hiVect()),
                      BL_TO_FORTRAN_ANYD(utilde_mf[mfi]),
                      utilde_mf.nComp(),
                      BL_TO_FORTRAN_ANYD(u_mf[mfi]),
                      BL_TO_FORTRAN_ANYD(v_mf[mfi]),
                      BL_TO_FORTRAN_ANYD(w_mf[mfi]),
                      BL_TO_FORTRAN_ANYD(Ipw[mfi]),
                      BL_TO_FORTRAN_ANYD(Imw[mfi]),
                      BL_TO_FORTRAN_ANYD(u_mf[mfi]),
                      BL_TO_FORTRAN_ANYD(v_mf[mfi]),
                      AMREX_INT_ANYD(domainBox.loVect()),
                      AMREX_INT_ANYD(domainBox.hiVect()),
                      bc_f, AMREX_REAL_ANYD(dx), dt, false,
                      3,3,AMREX_SPACEDIM,false);

              if (ppm_trace_forces == 1) {
#pragma gpu box(obx)
                ppm_3d(AMREX_INT_ANYD(obx.loVect()),
                    AMREX_INT_ANYD(obx.hiVect()),
                    BL_TO_FORTRAN_ANYD(force_mf[mfi]),
                    force_mf.nComp(),
                    BL_TO_FORTRAN_ANYD(u_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(v_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(w_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(Ipfz[mfi]),
                    BL_TO_FORTRAN_ANYD(Imfz[mfi]),
                    BL_TO_FORTRAN_ANYD(u_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(v_mf[mfi]),
                    AMREX_INT_ANYD(domainBox.loVect()),
                    AMREX_INT_ANYD(domainBox.hiVect()),
                    bc_f, AMREX_REAL_ANYD(dx), dt, false,
                    3,3,AMREX_SPACEDIM,false);
              }
            }

            // x-direction
#pragma gpu box(mxbx)
            velpred_interface_3d(AMREX_INT_ANYD(mxbx.loVect()), AMREX_INT_ANYD(mxbx.hiVect()),1,
                        AMREX_INT_ANYD(domainBox.loVect()), AMREX_INT_ANYD(domainBox.hiVect()),
                        BL_TO_FORTRAN_ANYD(utilde_mf[mfi]), utilde_mf.nComp(), utilde_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(ufull_mf[mfi]), ufull_mf.nComp(), ufull_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(utrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(Imu[mfi]),
                        BL_TO_FORTRAN_ANYD(Ipu[mfi]),
                        BL_TO_FORTRAN_ANYD(Imv[mfi]),
                        BL_TO_FORTRAN_ANYD(Ipv[mfi]),
                        BL_TO_FORTRAN_ANYD(Imw[mfi]),
                        BL_TO_FORTRAN_ANYD(Ipw[mfi]),
                        BL_TO_FORTRAN_ANYD(ulx[mfi]),
                        BL_TO_FORTRAN_ANYD(urx[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhx[mfi]),
                        AMREX_REAL_ANYD(dx), dt, phys_bc.dataPtr());

            // y-direction
#pragma gpu box(mybx)
            velpred_interface_3d(AMREX_INT_ANYD(mybx.loVect()), AMREX_INT_ANYD(mybx.hiVect()),2,
                        AMREX_INT_ANYD(domainBox.loVect()), AMREX_INT_ANYD(domainBox.hiVect()),
                        BL_TO_FORTRAN_ANYD(utilde_mf[mfi]), utilde_mf.nComp(), utilde_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(ufull_mf[mfi]), ufull_mf.nComp(), ufull_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(vtrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(Imu[mfi]),
                        BL_TO_FORTRAN_ANYD(Ipu[mfi]),
                        BL_TO_FORTRAN_ANYD(Imv[mfi]),
                        BL_TO_FORTRAN_ANYD(Ipv[mfi]),
                        BL_TO_FORTRAN_ANYD(Imw[mfi]),
                        BL_TO_FORTRAN_ANYD(Ipw[mfi]),
                        BL_TO_FORTRAN_ANYD(uly[mfi]),
                        BL_TO_FORTRAN_ANYD(ury[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhy[mfi]),
                        AMREX_REAL_ANYD(dx), dt, phys_bc.dataPtr());

            // z-direction
#pragma gpu box(mzbx)
            velpred_interface_3d(AMREX_INT_ANYD(mzbx.loVect()), AMREX_INT_ANYD(mzbx.hiVect()),3,
                        AMREX_INT_ANYD(domainBox.loVect()), AMREX_INT_ANYD(domainBox.hiVect()),
                        BL_TO_FORTRAN_ANYD(utilde_mf[mfi]), utilde_mf.nComp(), utilde_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(ufull_mf[mfi]), ufull_mf.nComp(), ufull_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(wtrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(Imu[mfi]),
                        BL_TO_FORTRAN_ANYD(Ipu[mfi]),
                        BL_TO_FORTRAN_ANYD(Imv[mfi]),
                        BL_TO_FORTRAN_ANYD(Ipv[mfi]),
                        BL_TO_FORTRAN_ANYD(Imw[mfi]),
                        BL_TO_FORTRAN_ANYD(Ipw[mfi]),
                        BL_TO_FORTRAN_ANYD(ulz[mfi]),
                        BL_TO_FORTRAN_ANYD(urz[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhz[mfi]),
                        AMREX_REAL_ANYD(dx), dt, phys_bc.dataPtr());

            // uimhyz, 1, 2
            Box imhbox = amrex::grow(mfi.tilebox(), 0, 1);
            imhbox = amrex::growHi(imhbox, 1, 1);
#pragma gpu box(imhbox)
            velpred_transverse_3d(AMREX_INT_ANYD(imhbox.loVect()), AMREX_INT_ANYD(imhbox.hiVect()),1,2,
                        AMREX_INT_ANYD(domainBox.loVect()), AMREX_INT_ANYD(domainBox.hiVect()),
                        BL_TO_FORTRAN_ANYD(utilde_mf[mfi]), utilde_mf.nComp(), utilde_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(utrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(vtrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(wtrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(ulx[mfi]),
                        BL_TO_FORTRAN_ANYD(urx[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhx[mfi]),
                        BL_TO_FORTRAN_ANYD(uly[mfi]),
                        BL_TO_FORTRAN_ANYD(ury[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhy[mfi]),
                        BL_TO_FORTRAN_ANYD(ulz[mfi]),
                        BL_TO_FORTRAN_ANYD(urz[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhz[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhyz[mfi]),
                        AMREX_REAL_ANYD(dx), dt, phys_bc.dataPtr());

            // uimhzy, 1, 3
            imhbox = amrex::grow(mfi.tilebox(), 0, 1);
            imhbox = amrex::growHi(imhbox, 2, 1);
#pragma gpu box(imhbox)
            velpred_transverse_3d(AMREX_INT_ANYD(imhbox.loVect()), AMREX_INT_ANYD(imhbox.hiVect()),1,3,
                        AMREX_INT_ANYD(domainBox.loVect()), AMREX_INT_ANYD(domainBox.hiVect()),
                        BL_TO_FORTRAN_ANYD(utilde_mf[mfi]), utilde_mf.nComp(), utilde_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(utrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(vtrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(wtrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(ulx[mfi]),
                        BL_TO_FORTRAN_ANYD(urx[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhx[mfi]),
                        BL_TO_FORTRAN_ANYD(uly[mfi]),
                        BL_TO_FORTRAN_ANYD(ury[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhy[mfi]),
                        BL_TO_FORTRAN_ANYD(ulz[mfi]),
                        BL_TO_FORTRAN_ANYD(urz[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhz[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhzy[mfi]),
                        AMREX_REAL_ANYD(dx), dt, phys_bc.dataPtr());

            // vimhxz, 2, 1
            imhbox = amrex::grow(mfi.tilebox(), 1, 1);
            imhbox = amrex::growHi(imhbox, 0, 1);
#pragma gpu box(imhbox)
            velpred_transverse_3d(AMREX_INT_ANYD(imhbox.loVect()), AMREX_INT_ANYD(imhbox.hiVect()),2,1,
                        AMREX_INT_ANYD(domainBox.loVect()), AMREX_INT_ANYD(domainBox.hiVect()),
                        BL_TO_FORTRAN_ANYD(utilde_mf[mfi]), utilde_mf.nComp(), utilde_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(utrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(vtrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(wtrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(ulx[mfi]),
                        BL_TO_FORTRAN_ANYD(urx[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhx[mfi]),
                        BL_TO_FORTRAN_ANYD(uly[mfi]),
                        BL_TO_FORTRAN_ANYD(ury[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhy[mfi]),
                        BL_TO_FORTRAN_ANYD(ulz[mfi]),
                        BL_TO_FORTRAN_ANYD(urz[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhz[mfi]),
                        BL_TO_FORTRAN_ANYD(vimhxz[mfi]),
                        AMREX_REAL_ANYD(dx), dt, phys_bc.dataPtr());

            // vimhxz, 2, 3
            imhbox = amrex::grow(mfi.tilebox(), 1, 1);
            imhbox = amrex::growHi(imhbox, 2, 1);
#pragma gpu box(imhbox)
            velpred_transverse_3d(AMREX_INT_ANYD(imhbox.loVect()), AMREX_INT_ANYD(imhbox.hiVect()),2,3,
                        AMREX_INT_ANYD(domainBox.loVect()), AMREX_INT_ANYD(domainBox.hiVect()),
                        BL_TO_FORTRAN_ANYD(utilde_mf[mfi]), utilde_mf.nComp(), utilde_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(utrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(vtrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(wtrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(ulx[mfi]),
                        BL_TO_FORTRAN_ANYD(urx[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhx[mfi]),
                        BL_TO_FORTRAN_ANYD(uly[mfi]),
                        BL_TO_FORTRAN_ANYD(ury[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhy[mfi]),
                        BL_TO_FORTRAN_ANYD(ulz[mfi]),
                        BL_TO_FORTRAN_ANYD(urz[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhz[mfi]),
                        BL_TO_FORTRAN_ANYD(vimhzx[mfi]),
                        AMREX_REAL_ANYD(dx), dt, phys_bc.dataPtr());

            // wimhxy, 3, 1
            imhbox = amrex::grow(mfi.tilebox(), 2, 1);
            imhbox = amrex::growHi(imhbox, 0, 1);
#pragma gpu box(imhbox)
            velpred_transverse_3d(AMREX_INT_ANYD(imhbox.loVect()), AMREX_INT_ANYD(imhbox.hiVect()),3,1,
                        AMREX_INT_ANYD(domainBox.loVect()), AMREX_INT_ANYD(domainBox.hiVect()),
                        BL_TO_FORTRAN_ANYD(utilde_mf[mfi]), utilde_mf.nComp(), utilde_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(utrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(vtrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(wtrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(ulx[mfi]),
                        BL_TO_FORTRAN_ANYD(urx[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhx[mfi]),
                        BL_TO_FORTRAN_ANYD(uly[mfi]),
                        BL_TO_FORTRAN_ANYD(ury[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhy[mfi]),
                        BL_TO_FORTRAN_ANYD(ulz[mfi]),
                        BL_TO_FORTRAN_ANYD(urz[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhz[mfi]),
                        BL_TO_FORTRAN_ANYD(wimhxy[mfi]),
                        AMREX_REAL_ANYD(dx), dt, phys_bc.dataPtr());

            // wimhyx, 3, 2
            imhbox = amrex::grow(mfi.tilebox(), 2, 1);
            imhbox = amrex::growHi(imhbox, 1, 1);
#pragma gpu box(imhbox)
            velpred_transverse_3d(AMREX_INT_ANYD(imhbox.loVect()), AMREX_INT_ANYD(imhbox.hiVect()),3,2,
                        AMREX_INT_ANYD(domainBox.loVect()), AMREX_INT_ANYD(domainBox.hiVect()),
                        BL_TO_FORTRAN_ANYD(utilde_mf[mfi]), utilde_mf.nComp(), utilde_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(utrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(vtrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(wtrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(ulx[mfi]),
                        BL_TO_FORTRAN_ANYD(urx[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhx[mfi]),
                        BL_TO_FORTRAN_ANYD(uly[mfi]),
                        BL_TO_FORTRAN_ANYD(ury[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhy[mfi]),
                        BL_TO_FORTRAN_ANYD(ulz[mfi]),
                        BL_TO_FORTRAN_ANYD(urz[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhz[mfi]),
                        BL_TO_FORTRAN_ANYD(wimhyx[mfi]),
                        AMREX_REAL_ANYD(dx), dt, phys_bc.dataPtr());

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.

            // x-direction
#pragma gpu box(xbx)
            velpred_3d(AMREX_INT_ANYD(xbx.loVect()),
                        AMREX_INT_ANYD(xbx.hiVect()),
                        lev, 1, AMREX_INT_ANYD(domainBox.loVect()),
                        AMREX_INT_ANYD(domainBox.hiVect()),
                        BL_TO_FORTRAN_ANYD(utilde_mf[mfi]),
                        utilde_mf.nComp(), utilde_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(utrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(vtrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(wtrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(wmac_mf[mfi]),
            			BL_TO_FORTRAN_ANYD(w0macx_mf[mfi]),
            			BL_TO_FORTRAN_ANYD(w0macy_mf[mfi]),
            			BL_TO_FORTRAN_ANYD(w0macz_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(Imfx[mfi]),
                        BL_TO_FORTRAN_ANYD(Ipfx[mfi]),
                        BL_TO_FORTRAN_ANYD(ulx[mfi]),
                        BL_TO_FORTRAN_ANYD(urx[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhyz[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhzy[mfi]),
                        BL_TO_FORTRAN_ANYD(vimhxz[mfi]),
                        BL_TO_FORTRAN_ANYD(vimhzx[mfi]),
                        BL_TO_FORTRAN_ANYD(wimhxy[mfi]),
                        BL_TO_FORTRAN_ANYD(wimhyx[mfi]),
                        BL_TO_FORTRAN_ANYD(force_mf[mfi]), force_mf.nComp(), force_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(w0_mf[mfi]), 
                        AMREX_REAL_ANYD(dx), dt, phys_bc.dataPtr());

            // y-direction
#pragma gpu box(ybx)
            velpred_3d(AMREX_INT_ANYD(ybx.loVect()),
                        AMREX_INT_ANYD(ybx.hiVect()),
                        lev, 2, AMREX_INT_ANYD(domainBox.loVect()),
                        AMREX_INT_ANYD(domainBox.hiVect()),
                        BL_TO_FORTRAN_ANYD(utilde_mf[mfi]),
                        utilde_mf.nComp(), utilde_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(utrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(vtrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(wtrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(wmac_mf[mfi]),
            			BL_TO_FORTRAN_ANYD(w0macx_mf[mfi]),
            			BL_TO_FORTRAN_ANYD(w0macy_mf[mfi]),
            			BL_TO_FORTRAN_ANYD(w0macz_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(Imfy[mfi]),
                        BL_TO_FORTRAN_ANYD(Ipfy[mfi]),
                        BL_TO_FORTRAN_ANYD(uly[mfi]),
                        BL_TO_FORTRAN_ANYD(ury[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhyz[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhzy[mfi]),
                        BL_TO_FORTRAN_ANYD(vimhxz[mfi]),
                        BL_TO_FORTRAN_ANYD(vimhzx[mfi]),
                        BL_TO_FORTRAN_ANYD(wimhxy[mfi]),
                        BL_TO_FORTRAN_ANYD(wimhyx[mfi]),
                        BL_TO_FORTRAN_ANYD(force_mf[mfi]), force_mf.nComp(), force_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(w0_mf[mfi]), 
                        AMREX_REAL_ANYD(dx), dt, phys_bc.dataPtr());

            // z-direction
#pragma gpu box(zbx)
            velpred_3d(AMREX_INT_ANYD(zbx.loVect()),
                        AMREX_INT_ANYD(zbx.hiVect()),
                        lev, 3, AMREX_INT_ANYD(domainBox.loVect()),
                        AMREX_INT_ANYD(domainBox.hiVect()),
                        BL_TO_FORTRAN_ANYD(utilde_mf[mfi]),
                        utilde_mf.nComp(), utilde_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(utrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(vtrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(wtrans_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(wmac_mf[mfi]),
            			BL_TO_FORTRAN_ANYD(w0macx_mf[mfi]),
            			BL_TO_FORTRAN_ANYD(w0macy_mf[mfi]),
            			BL_TO_FORTRAN_ANYD(w0macz_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(Imfz[mfi]),
                        BL_TO_FORTRAN_ANYD(Ipfz[mfi]),
                        BL_TO_FORTRAN_ANYD(ulz[mfi]),
                        BL_TO_FORTRAN_ANYD(urz[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhyz[mfi]),
                        BL_TO_FORTRAN_ANYD(uimhzy[mfi]),
                        BL_TO_FORTRAN_ANYD(vimhxz[mfi]),
                        BL_TO_FORTRAN_ANYD(vimhzx[mfi]),
                        BL_TO_FORTRAN_ANYD(wimhxy[mfi]),
                        BL_TO_FORTRAN_ANYD(wimhyx[mfi]),
                        BL_TO_FORTRAN_ANYD(force_mf[mfi]), force_mf.nComp(), force_mf.nGrow(),
                        BL_TO_FORTRAN_ANYD(w0_mf[mfi]), 
                        AMREX_REAL_ANYD(dx), dt, phys_bc.dataPtr());
        } // end MFIter loop

#ifdef AMREX_USE_CUDA
        clean_bc(bc_f);
#endif

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
        MultiFab& sedgex_mf = sedge[lev][0];
        const MultiFab& umac_mf   = umac[lev][0];
        MultiFab& sedgey_mf = sedge[lev][1];
        const MultiFab& vmac_mf   = umac[lev][1];

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
        MultiFab& sedgez_mf = sedge[lev][2];
        const MultiFab& wmac_mf   = umac[lev][2];

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

        MultiFab sl, sr;
        sl.define(grids[lev],dmap[lev],1,1);
        sr.define(grids[lev],dmap[lev],1,1);
        sl.setVal(0.);
        sr.setVal(0);


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

#ifdef AMREX_USE_CUDA
        int* bc_f = prepare_bc(bcs[0].data(), 1);
#else
        const int* bc_f = bcs[0].data();
#endif
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
            const Box& xbx = mfi.nodaltilebox(0);
            const Box& ybx = mfi.nodaltilebox(1);
            const Box& mxbx = amrex::growLo(obx,0, -1);
            const Box& mybx = amrex::growLo(obx,1, -1);

            // Be careful to pass in comp+1 for fortran indexing
            for (int scomp = start_scomp+1; scomp <= start_scomp + num_comp; ++scomp) {

                int vcomp = scomp - start_scomp - 1;
                int bccomp = start_bccomp + scomp - start_scomp;

                // x-direction
                if (ppm_type == 0) {
                    // we're going to reuse Ip here as slopex and Im as slopey
                    // as they have the correct number of ghost zones

                    // x-direction
#pragma gpu box(obx)
                    slopex_2d(AMREX_INT_ANYD(obx.loVect()),
                              AMREX_INT_ANYD(obx.hiVect()),
                              BL_TO_FORTRAN_ANYD(vec_scal_mf[vcomp][mfi]),
                              vec_scal_mf[vcomp].nComp(),
                              BL_TO_FORTRAN_ANYD(Ip[mfi]),Ip.nComp(),
                              AMREX_INT_ANYD(domainBox.loVect()),
                              AMREX_INT_ANYD(domainBox.hiVect()),
                              1,bc_f,nbccomp,bccomp);

                    // y-direction
#pragma gpu box(obx)
                    slopey_2d(AMREX_INT_ANYD(obx.loVect()),
                              AMREX_INT_ANYD(obx.hiVect()),
                              BL_TO_FORTRAN_ANYD(vec_scal_mf[vcomp][mfi]),
                              vec_scal_mf[vcomp].nComp(),
                              BL_TO_FORTRAN_ANYD(Im[mfi]),Im.nComp(),
                              AMREX_INT_ANYD(domainBox.loVect()),
                              AMREX_INT_ANYD(domainBox.hiVect()),
                              1,bc_f,nbccomp,bccomp);

                } else {
#pragma gpu box(obx)
                    ppm_2d(AMREX_INT_ANYD(obx.loVect()),
                           AMREX_INT_ANYD(obx.hiVect()),
                           BL_TO_FORTRAN_ANYD(scal_mf[mfi]),
                           scal_mf.nComp(),
                           BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
                           BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
                           BL_TO_FORTRAN_ANYD(Ip[mfi]),
                           BL_TO_FORTRAN_ANYD(Im[mfi]),
                           AMREX_INT_ANYD(domainBox.loVect()),
                           AMREX_INT_ANYD(domainBox.hiVect()),
                           bc_f, AMREX_REAL_ANYD(dx), dt, true,
                           scomp, bccomp, nbccomp);

                    if (ppm_trace_forces == 1) {
#pragma gpu box(obx)
                        ppm_2d(AMREX_INT_ANYD(obx.loVect()),
                               AMREX_INT_ANYD(obx.hiVect()),
                               BL_TO_FORTRAN_ANYD(force_mf[mfi]),
                               force_mf.nComp(),
                               BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
                               BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
                               BL_TO_FORTRAN_ANYD(Ipf[mfi]),
                               BL_TO_FORTRAN_ANYD(Imf[mfi]),
                               AMREX_INT_ANYD(domainBox.loVect()),
                               AMREX_INT_ANYD(domainBox.hiVect()),
                               bc_f, AMREX_REAL_ANYD(dx), dt, true,
                               scomp, bccomp, nbccomp);

                    }
                }

                // call fortran subroutine
                // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
                // lo/hi coordinates (including ghost cells), and/or the # of components
                // We will also pass "validBox", which specifies the "valid" region.

                // x-direction
#pragma gpu box(mxbx)
                make_edge_scal_predictor_2d(
                    AMREX_INT_ANYD(mxbx.loVect()), AMREX_INT_ANYD(mxbx.hiVect()), 1,
                    AMREX_INT_ANYD(domainBox.loVect()), AMREX_INT_ANYD(domainBox.hiVect()),
                    BL_TO_FORTRAN_ANYD(scal_mf[mfi]), scal_mf.nComp(), scal_mf.nGrow(),
                    BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(Ip[mfi]),
                    BL_TO_FORTRAN_ANYD(Im[mfi]),
                    BL_TO_FORTRAN_ANYD(slx[mfi]),
                    BL_TO_FORTRAN_ANYD(srx[mfi]),
                    BL_TO_FORTRAN_ANYD(simhx[mfi]),
                    AMREX_REAL_ANYD(dx), dt, is_vel, bc_f,
                    nbccomp, scomp, bccomp);

                // y-direction
#pragma gpu box(mybx)
                make_edge_scal_predictor_2d(
                    AMREX_INT_ANYD(mybx.loVect()), AMREX_INT_ANYD(mybx.hiVect()), 2,
                    AMREX_INT_ANYD(domainBox.loVect()), AMREX_INT_ANYD(domainBox.hiVect()),
                    BL_TO_FORTRAN_ANYD(scal_mf[mfi]), scal_mf.nComp(), scal_mf.nGrow(),
                    BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(Ip[mfi]),
                    BL_TO_FORTRAN_ANYD(Im[mfi]),
                    BL_TO_FORTRAN_ANYD(sly[mfi]),
                    BL_TO_FORTRAN_ANYD(sry[mfi]),
                    BL_TO_FORTRAN_ANYD(simhy[mfi]),
                    AMREX_REAL_ANYD(dx), dt, is_vel, bc_f,
                    nbccomp, scomp, bccomp);

                // x-direction
#pragma gpu box(xbx)
                make_edge_scal_2d(
                    AMREX_INT_ANYD(xbx.loVect()), AMREX_INT_ANYD(xbx.hiVect()),1,
                    AMREX_INT_ANYD(domainBox.loVect()), AMREX_INT_ANYD(domainBox.hiVect()),
                    BL_TO_FORTRAN_ANYD(scal_mf[mfi]), scal_mf.nComp(), scal_mf.nGrow(),
                    BL_TO_FORTRAN_ANYD(sedgex_mf[mfi]), sedgex_mf.nComp(),
                    BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(Ipf[mfi]),
                    BL_TO_FORTRAN_ANYD(Imf[mfi]),
                    BL_TO_FORTRAN_ANYD(slx[mfi]),
                    BL_TO_FORTRAN_ANYD(srx[mfi]),
                    BL_TO_FORTRAN_ANYD(simhy[mfi]),
                    BL_TO_FORTRAN_ANYD(force_mf[mfi]), force_mf.nComp(),
                    AMREX_REAL_ANYD(dx), dt, is_vel, bc_f,
                    nbccomp, scomp, bccomp, is_conservative);

                // y-direction
#pragma gpu box(ybx)
                make_edge_scal_2d(
                    AMREX_INT_ANYD(ybx.loVect()), AMREX_INT_ANYD(ybx.hiVect()),2,
                    AMREX_INT_ANYD(domainBox.loVect()), AMREX_INT_ANYD(domainBox.hiVect()),
                    BL_TO_FORTRAN_ANYD(scal_mf[mfi]), scal_mf.nComp(), scal_mf.nGrow(),
                    BL_TO_FORTRAN_ANYD(sedgey_mf[mfi]), sedgey_mf.nComp(),
                    BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(Ipf[mfi]),
                    BL_TO_FORTRAN_ANYD(Imf[mfi]),
                    BL_TO_FORTRAN_ANYD(sly[mfi]),
                    BL_TO_FORTRAN_ANYD(sry[mfi]),
                    BL_TO_FORTRAN_ANYD(simhx[mfi]),
                    BL_TO_FORTRAN_ANYD(force_mf[mfi]), force_mf.nComp(),
                    AMREX_REAL_ANYD(dx), dt, is_vel, bc_f,
                    nbccomp, scomp, bccomp, is_conservative);
            } // end loop over components
        } // end MFIter loop

#ifdef AMREX_USE_CUDA
        clean_bc(bc_f);
#endif

#elif (AMREX_SPACEDIM == 3)

#ifdef AMREX_USE_CUDA
        int* bc_f = prepare_bc(bcs[0].data(), 1);
#else
        const int* bc_f = bcs[0].data();
#endif
        Vector<MultiFab> vec_scal_mf(num_comp);
        for (int comp=0; comp < num_comp; ++comp) {
            vec_scal_mf[comp].define(grids[lev],dmap[lev],1,scal_mf.nGrow());
            vec_scal_mf[comp].setVal(0.);

            MultiFab::Copy(vec_scal_mf[comp], scal_mf, start_scomp+comp, 0, 1, scal_mf.nGrow());
        }

        // Be careful to pass in comp+1 for fortran indexing
        for (int scomp = start_scomp+1; scomp <= start_scomp + num_comp; ++scomp) {

            int vcomp = scomp - start_scomp - 1;

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

                // make divu 
                if (is_conservative == 1) {

                    Array4<Real> const divu_arr = divu.array(mfi);
                    Array4<Real> const umac_arr = (umac[lev][0]).array(mfi);
                    Array4<Real> const vmac_arr = (umac[lev][1]).array(mfi);
                    Array4<Real> const wmac_arr = (umac[lev][2]).array(mfi);

                    AMREX_PARALLEL_FOR_3D(obx, i, j, k, 
                    {
                        divu_arr(i,j,k) = 
                            umac_arr(i+1,j,k) - umac_arr(i,j,k) +
                            vmac_arr(i,j+1,k) - vmac_arr(i,j,k) +
                            wmac_arr(i,j,k+1) - wmac_arr(i,j,k);
                        divu_arr(i,j,k) /= dx[0];
                    });
                }
                          
                // x-direction
                if (ppm_type == 0) {
                    // we're going to reuse Ip here as slopex and Im as slopey
                    // as they have the correct number of ghost zones

                    // x-direction
#pragma gpu box(obx)
                    slopex_2d(AMREX_INT_ANYD(obx.loVect()),
                              AMREX_INT_ANYD(obx.hiVect()),
                              BL_TO_FORTRAN_ANYD(vec_scal_mf[vcomp][mfi]),
                              vec_scal_mf[vcomp].nComp(),
                              BL_TO_FORTRAN_ANYD(Ip[mfi]),Ip.nComp(),
                              AMREX_INT_ANYD(domainBox.loVect()),
                              AMREX_INT_ANYD(domainBox.hiVect()),
                              1,bc_f,nbccomp,bccomp);

                    // y-direction
#pragma gpu box(obx)
                    slopey_2d(AMREX_INT_ANYD(obx.loVect()),
                              AMREX_INT_ANYD(obx.hiVect()),
                              BL_TO_FORTRAN_ANYD(vec_scal_mf[vcomp][mfi]),
                              vec_scal_mf[vcomp].nComp(),
                              BL_TO_FORTRAN_ANYD(Im[mfi]),Im.nComp(),
                              AMREX_INT_ANYD(domainBox.loVect()),
                              AMREX_INT_ANYD(domainBox.hiVect()),
                              1,bc_f,nbccomp,bccomp);

                    // z-direction
#pragma gpu box(obx)
                    slopez_3d(AMREX_INT_ANYD(obx.loVect()),
                              AMREX_INT_ANYD(obx.hiVect()),
                              BL_TO_FORTRAN_ANYD(vec_scal_mf[vcomp][mfi]),
                              vec_scal_mf[vcomp].nComp(),
                              BL_TO_FORTRAN_ANYD(slopez[mfi]),slopez.nComp(),
                              AMREX_INT_ANYD(domainBox.loVect()),
                              AMREX_INT_ANYD(domainBox.hiVect()),
                              1,bc_f,nbccomp,bccomp);


                } else {

                    Array4<Real> const scal_arr = state[lev].array(mfi);

                    Array4<Real> const umac_arr = (umac[lev][0]).array(mfi);
                    Array4<Real> const vmac_arr = (umac[lev][1]).array(mfi);
                    Array4<Real> const wmac_arr = (umac[lev][2]).array(mfi);

                    Array4<Real> const Ip_arr = Ip.array(mfi);
                    Array4<Real> const Im_arr = Im.array(mfi);
                    Real rel_eps;
                    get_rel_eps(&rel_eps);

                    Array4<Real> const sp_arr = sl.array(mfi);
                    Array4<Real> const sm_arr = sr.array(mfi);

                    PPM_3d(obx, scal_arr, 
                           umac_arr, vmac_arr, wmac_arr,
                           Ip_arr, Im_arr, 
                           sp_arr, sm_arr, 
                           domainBox, bcs, dx, 
                           true, scomp-1, bccomp-1, rel_eps);

#pragma gpu box(obx)
                    ppm_3d(AMREX_INT_ANYD(obx.loVect()),
                           AMREX_INT_ANYD(obx.hiVect()),
                           BL_TO_FORTRAN_ANYD(scal_mf[mfi]),
                           scal_mf.nComp(),
                           BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
                           BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
                           BL_TO_FORTRAN_ANYD(wmac_mf[mfi]),
                           BL_TO_FORTRAN_ANYD(Ip[mfi]),
                           BL_TO_FORTRAN_ANYD(Im[mfi]),
                           BL_TO_FORTRAN_ANYD(sl[mfi]),
                           BL_TO_FORTRAN_ANYD(sr[mfi]),
                           AMREX_INT_ANYD(domainBox.loVect()),
                           AMREX_INT_ANYD(domainBox.hiVect()),
                           bc_f, AMREX_REAL_ANYD(dx), dt, true,
                           scomp, bccomp, nbccomp, true);

                    if (ppm_trace_forces == 1) {
#pragma gpu box(obx)
                        ppm_3d(AMREX_INT_ANYD(obx.loVect()),
                               AMREX_INT_ANYD(obx.hiVect()),
                               BL_TO_FORTRAN_ANYD(force_mf[mfi]),
                               force_mf.nComp(),
                               BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
                               BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
                               BL_TO_FORTRAN_ANYD(wmac_mf[mfi]),
                               BL_TO_FORTRAN_ANYD(Ipf[mfi]),
                               BL_TO_FORTRAN_ANYD(Imf[mfi]),
                               BL_TO_FORTRAN_ANYD(sl[mfi]),
                               BL_TO_FORTRAN_ANYD(sr[mfi]),
                               AMREX_INT_ANYD(domainBox.loVect()),
                               AMREX_INT_ANYD(domainBox.hiVect()),
                               bc_f, AMREX_REAL_ANYD(dx), dt, true,
                               scomp, bccomp, nbccomp, false);
                    }
                }
            }


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

                Array4<Real> const scal_arr = state[lev].array(mfi);

                Array4<Real> const umac_arr = (umac[lev][0]).array(mfi);
                Array4<Real> const vmac_arr = (umac[lev][1]).array(mfi);
                Array4<Real> const wmac_arr = (umac[lev][2]).array(mfi);

                Array4<Real> const slx_arr = slx.array(mfi);
                Array4<Real> const srx_arr = srx.array(mfi);
                Array4<Real> const sly_arr = sly.array(mfi);
                Array4<Real> const sry_arr = sry.array(mfi);
                Array4<Real> const slz_arr = slz.array(mfi);
                Array4<Real> const srz_arr = srz.array(mfi);

                Array4<Real> const simhx_arr = simhx.array(mfi);
                Array4<Real> const simhy_arr = simhy.array(mfi);
                Array4<Real> const simhz_arr = simhz.array(mfi);

                Array4<Real> const Ip_arr = Ip.array(mfi);
                Array4<Real> const Im_arr = Im.array(mfi);
                Array4<Real> const slopez_arr = slopez.array(mfi);

                Array4<Real> const sl_arr = sl.array(mfi);
                Array4<Real> const sr_arr = sr.array(mfi);

                Real ppm_type_local = ppm_type;
                Real rel_eps;
                get_rel_eps(&rel_eps);

                Real hx = dx[0];
                Real hy = dx[1];
                Real hz = dx[2];

                // loop over appropriate x-faces
                AMREX_PARALLEL_FOR_3D(mxbx, i, j, k, 
                {
                    if (ppm_type_local == 0) {
                        slx_arr(i,j,k) = scal_arr(i-1,j,k,scomp-1) + 
                            0.5 * (1.0 - dt * umac_arr(i,j,k) / hx) * Ip_arr(i-1,j,k,0);
                        srx_arr(i,j,k) = scal_arr(i,j,k,scomp-1) - 
                            0.5 * (1.0 + dt * umac_arr(i,j,k) / hx) * Ip_arr(i,j,k,0);
                    } else if (ppm_type_local == 1 || ppm_type_local == 2) {
                        slx_arr(i,j,k) = Ip_arr(i-1,j,k,0);
                        srx_arr(i,j,k) = Im_arr(i,j,k,0);
                    }
                });

                // loop over appropriate y-faces
                AMREX_PARALLEL_FOR_3D(mybx, i, j, k, 
                {
                    if (ppm_type_local == 0) {
                        sly_arr(i,j,k) = scal_arr(i,j-1,k,scomp-1) + 
                            0.5 * (1.0 - dt * vmac_arr(i,j,k) / hy) * Ip_arr(i,j-1,k,0);
                        sry_arr(i,j,k) = scal_arr(i,j,k,scomp-1) - 
                            0.5 * (1.0 + dt * vmac_arr(i,j,k) / hy) * Ip_arr(i,j,k,0);
                    } else if (ppm_type_local == 1 || ppm_type_local == 2) {
                        sly_arr(i,j,k) = Ip_arr(i,j-1,k,1);
                        sry_arr(i,j,k) = Im_arr(i,j,k,1);
                    }
                });

                // loop over appropriate z-faces
                AMREX_PARALLEL_FOR_3D(mzbx, i, j, k, 
                {
                    if (ppm_type_local == 0) {
                        slz_arr(i,j,k) = scal_arr(i,j,k-1,scomp-1) + 
                            0.5 * (1.0 - dt * wmac_arr(i,j,k) / hz) * slopez_arr(i,j,k-1);
                        srz_arr(i,j,k) = scal_arr(i,j,k,scomp-1) - 
                            0.5 * (1.0 + dt * wmac_arr(i,j,k) / hz) * slopez_arr(i,j,k);
                    } else if (ppm_type_local == 1 || ppm_type_local == 2) {
                        slz_arr(i,j,k) = Ip_arr(i,j,k-1,2);
                        srz_arr(i,j,k) = Im_arr(i,j,k,2);
                    }
                });
                
                int ilo = domainBox.loVect()[0];
                int ihi = domainBox.hiVect()[0];
                int bclo = bcs[bccomp-1].lo()[0];
                int bchi = bcs[bccomp-1].hi()[0];
                AMREX_PARALLEL_FOR_3D(mxbx, i, j, k, 
                {
                    // impose lo side bc's
                    if (i == ilo) {
                        if (bclo == EXT_DIR) {
                            slx_arr(i,j,k) = scal_arr(i-1,j,k,scomp-1);
                            srx_arr(i,j,k) = scal_arr(i-1,j,k,scomp-1);
                        } else if (bclo == FOEXTRAP || bclo == HOEXTRAP) {
                            if (is_vel == 1 && scomp-1 == 0) {
                                srx_arr(i,j,k) = min(srx_arr(i,j,k), 0.0);
                            }
                            slx_arr(i,j,k) = srx_arr(i,j,k);
                        } else if (bclo == REFLECT_EVEN) {
                            slx_arr(i,j,k) = srx_arr(i,j,k);
                        } else if (bclo == REFLECT_ODD) {
                            slx_arr(i,j,k) = 0.0;
                            srx_arr(i,j,k) = 0.0;
                        }
                    // impose hi side bc's
                    } else if (i == ihi+1) {
                        if (bchi == EXT_DIR) {
                            slx_arr(i,j,k) = scal_arr(i,j,k,scomp-1);
                            srx_arr(i,j,k) = scal_arr(i,j,k,scomp-1);
                        } else if (bchi == FOEXTRAP || bchi == HOEXTRAP) {
                            if (is_vel == 1 && scomp-1 == 0) {
                                slx_arr(i,j,k) = max(slx_arr(i,j,k), 0.0);
                            }
                            srx_arr(i,j,k) = slx_arr(i,j,k);
                        } else if (bchi == REFLECT_EVEN) {
                            srx_arr(i,j,k) = slx_arr(i,j,k);
                        } else if (bchi == REFLECT_ODD) {
                            slx_arr(i,j,k) = 0.0;
                            srx_arr(i,j,k) = 0.0;
                        }

                    }
                });
                
                int jlo = domainBox.loVect()[1];
                int jhi = domainBox.hiVect()[1];
                bclo = bcs[bccomp-1].lo()[1];
                bchi = bcs[bccomp-1].hi()[1];
                AMREX_PARALLEL_FOR_3D(mybx, i, j, k, 
                {
                    // impose lo side bc's
                    if (j == jlo) {
                        if (bclo == EXT_DIR) {
                            sly_arr(i,j,k) = scal_arr(i,j-1,k,scomp-1);
                            sry_arr(i,j,k) = scal_arr(i,j-1,k,scomp-1);
                        } else if (bclo == FOEXTRAP || bclo == HOEXTRAP) {
                            if (is_vel == 1 && scomp-1 == 1) {
                                sry_arr(i,j,k) = min(sry_arr(i,j,k), 0.0);
                            }
                            sly_arr(i,j,k) = sry_arr(i,j,k);
                        } else if (bclo == REFLECT_EVEN) {
                            sly_arr(i,j,k) = sry_arr(i,j,k);
                        } else if (bclo == REFLECT_ODD) {
                            sly_arr(i,j,k) = 0.0;
                            sry_arr(i,j,k) = 0.0;
                        }
                    // impose hi side bc's
                    } else if (j == jhi+1) {
                        if (bchi == EXT_DIR) {
                            sly_arr(i,j,k) = scal_arr(i,j,k,scomp-1);
                            sry_arr(i,j,k) = scal_arr(i,j,k,scomp-1);
                        } else if (bchi == FOEXTRAP || bchi == HOEXTRAP) {
                            if (is_vel == 1 && scomp-1 == 1) {
                                sly_arr(i,j,k) = max(sly_arr(i,j,k), 0.0);
                            }
                            sry_arr(i,j,k) = sly_arr(i,j,k);
                        } else if (bchi == REFLECT_EVEN) {
                            sry_arr(i,j,k) = sly_arr(i,j,k);
                        } else if (bchi == REFLECT_ODD) {
                            sly_arr(i,j,k) = 0.0;
                            sry_arr(i,j,k) = 0.0;
                        }

                    }
                });
                
                int klo = domainBox.loVect()[2];
                int khi = domainBox.hiVect()[2];
                bclo = bcs[bccomp-1].lo()[2];
                bchi = bcs[bccomp-1].hi()[2];
                AMREX_PARALLEL_FOR_3D(mzbx, i, j, k, 
                {
                    // impose lo side bc's
                    if (k == klo) {
                        if (bclo == EXT_DIR) {
                            slz_arr(i,j,k) = scal_arr(i,j,k-1,scomp-1);
                            srz_arr(i,j,k) = scal_arr(i,j,k-1,scomp-1);
                        } else if (bclo == FOEXTRAP || bclo == HOEXTRAP) {
                            if (is_vel == 1 && scomp-1 == 2) {
                                srz_arr(i,j,k) = min(srz_arr(i,j,k), 0.0);
                            }
                            slz_arr(i,j,k) = srz_arr(i,j,k);
                        } else if (bclo == REFLECT_EVEN) {
                            slz_arr(i,j,k) = srz_arr(i,j,k);
                        } else if (bclo == REFLECT_ODD) {
                            slz_arr(i,j,k) = 0.0;
                            srz_arr(i,j,k) = 0.0;
                        }
                    // impose hi side bc's
                    } else if (k == khi+1) {
                        if (bchi == EXT_DIR) {
                            slz_arr(i,j,k) = scal_arr(i,j,k,scomp-1);
                            srz_arr(i,j,k) = scal_arr(i,j,k,scomp-1);
                        } else if (bchi == FOEXTRAP || bchi == HOEXTRAP) {
                            if (is_vel == 1 && scomp-1 == 2) {
                                slz_arr(i,j,k) = max(slz_arr(i,j,k), 0.0);
                            }
                            srz_arr(i,j,k) = slz_arr(i,j,k);
                        } else if (bchi == REFLECT_EVEN) {
                            srz_arr(i,j,k) = slz_arr(i,j,k);
                        } else if (bchi == REFLECT_ODD) {
                            slz_arr(i,j,k) = 0.0;
                            srz_arr(i,j,k) = 0.0;
                        }

                    }
                });

                // make simhx by solving Riemann problem
                AMREX_PARALLEL_FOR_3D(mxbx, i, j, k, 
                {
                    simhx_arr(i,j,k) = (umac_arr(i,j,k) > 0.0) ? 
                        slx_arr(i,j,k) : srx_arr(i,j,k);
                    simhx_arr(i,j,k) = (fabs(umac_arr(i,j,k)) > 0.0) ? 
                        simhx_arr(i,j,k) : 0.5 * (slx_arr(i,j,k) + srx_arr(i,j,k));
                });

                // make simhy by solving Riemann problem
                AMREX_PARALLEL_FOR_3D(mybx, i, j, k, 
                {
                    simhy_arr(i,j,k) = (vmac_arr(i,j,k) > 0.0) ? 
                        sly_arr(i,j,k) : sry_arr(i,j,k);
                    simhy_arr(i,j,k) = (fabs(vmac_arr(i,j,k)) > 0.0) ? 
                        simhy_arr(i,j,k) : 0.5 * (sly_arr(i,j,k) + sry_arr(i,j,k));
                });

                // make simhz by solving Riemann problem
                AMREX_PARALLEL_FOR_3D(mzbx, i, j, k, 
                {
                    simhz_arr(i,j,k) = (wmac_arr(i,j,k) > 0.0) ? 
                        slz_arr(i,j,k) : srz_arr(i,j,k);
                    simhz_arr(i,j,k) = (fabs(wmac_arr(i,j,k)) > 0.0) ?
                        simhz_arr(i,j,k) : 0.5 * (slz_arr(i,j,k) + srz_arr(i,j,k));
                });

                ////////////////////////////////////////////////////////
                // Create transverse terms, s_{\i-\half\e_x}^{x|y}, etc.
                ////////////////////////////////////////////////////////
                Real dt2 = 0.5 * dt;
                Real dt3 = dt / 3.0;
                Real dt4 = 0.25 * dt;
                Real dt6 = dt / 6.0;

                Array4<Real> const divu_arr = divu.array(mfi);
                Array4<Real> const simhxy_arr = simhxy.array(mfi);

                // simhxy
                // Box imhbox = amrex::grow(mfi.tilebox(), 2, 1);
                // imhbox = amrex::growHi(imhbox, 0, 1);
                Box imhbox = mfi.grownnodaltilebox(0, amrex::IntVect(0,0,1)); 
                bclo = bcs[bccomp-1].lo()[0];
                bchi = bcs[bccomp-1].hi()[0];
                AMREX_PARALLEL_FOR_3D(imhbox, i, j, k, 
                {
                    Real slxy = 0.0;
                    Real srxy = 0.0;
                    
                    // loop over appropriate xy faces
                    if (is_conservative == 1) {
                        // make slxy, srxy by updating 1D extrapolation
                        slxy = slx_arr(i,j,k) 
                            - (dt3/hy) * (simhy_arr(i-1,j+1,k)*vmac_arr(i-1,j+1,k) 
                            - simhy_arr(i-1,j,k)*vmac_arr(i-1,j,k)) 
                            - dt3*scal_arr(i-1,j,k,scomp-1)*divu_arr(i-1,j,k) 
                            + (dt3/hy)*scal_arr(i-1,j,k,scomp-1)*
                            (vmac_arr(i-1,j+1,k)-vmac_arr(i-1,j,k));
                        srxy = srx_arr(i,j,k) 
                            - (dt3/hy)*(simhy_arr(i,j+1,k)*vmac_arr(i,j+1,k)
                            - simhy_arr(i,j,k)*vmac_arr(i,j,k)) 
                            - dt3*scal_arr(i,j,k,scomp-1)*divu_arr(i,j,k) 
                            + (dt3/hy)*scal_arr(i,j,k,scomp-1)*
                            (vmac_arr(i,j+1,k)-vmac_arr(i,j,k));
                    } else {
                        // make slxy, srxy by updating 1D extrapolation
                        slxy = slx_arr(i,j,k) 
                            - (dt6/hy)*(vmac_arr(i-1,j+1,k)+vmac_arr(i-1,j,k)) 
                            *(simhy_arr(i-1,j+1,k)-simhy_arr(i-1,j,k));
                        srxy = srx_arr(i,j,k) 
                            - (dt6/hy)*(vmac_arr(i,j+1,k)+vmac_arr(i,j,k))
                            *(simhy_arr(i,j+1,k)-simhy_arr(i,j,k));
                    }

                    // impose lo side bc's
                    if (i == ilo) {
                        if (bclo == EXT_DIR) {
                            slxy = scal_arr(i-1,j,k,scomp-1);
                            srxy = scal_arr(i-1,j,k,scomp-1);
                        } else if (bclo ==  FOEXTRAP || bclo == HOEXTRAP) {
                            if (is_vel == 1 && scomp-1 == 0) {
                                srxy = min(srxy,0.0);
                            }
                            slxy = srxy;
                        } else if (bclo == REFLECT_EVEN) {
                            slxy = srxy;
                        } else if (bclo ==  REFLECT_ODD) {
                            slxy = 0.0;
                            srxy = 0.0;
                        }

                    // impose hi side bc's
                    } else if (i == ihi+1) {
                        if (bchi == EXT_DIR) {
                            slxy = scal_arr(i,j,k,scomp-1);
                            srxy = scal_arr(i,j,k,scomp-1);
                        } else if (bchi == FOEXTRAP || bchi == HOEXTRAP) {
                            if (is_vel == 1 && scomp-1 == 0) {
                                slxy = max(slxy,0.0);
                            }
                            srxy = slxy;
                        } else if (bchi == REFLECT_EVEN) {
                            srxy = slxy;
                        } else if (bchi == REFLECT_ODD) {
                            slxy = 0.0;
                            srxy = 0.0;
                        }
                    }

                    // make simhxy by solving Riemann problem
                    simhxy_arr(i,j,k) = (umac_arr(i,j,k) > 0.0) ?
                        slxy : srxy;
                    simhxy_arr(i,j,k) = (fabs(umac_arr(i,j,k)) > rel_eps) ?
                        simhxy_arr(i,j,k) : 0.5 * (slxy + srxy);
     
                });

                // simhxz
                // imhbox = amrex::grow(mfi.tilebox(), 1, 1);
                // imhbox = amrex::growHi(imhbox, 0, 1);
                imhbox = mfi.grownnodaltilebox(0, amrex::IntVect(0,1,0));

                Array4<Real> const simhxz_arr = simhxz.array(mfi);

                AMREX_PARALLEL_FOR_3D(imhbox, i, j, k, 
                {
                    Real slxz = 0.0;
                    Real srxz = 0.0;
                    // loop over appropriate xz faces
                    if (is_conservative == 1) {
                        // make slxz, srxz by updating 1D extrapolation
                        slxz = slx_arr(i,j,k) 
                            - (dt3/hz) * (simhz_arr(i-1,j,k+1)*wmac_arr(i-1,j,k+1) 
                            - simhz_arr(i-1,j,k)*wmac_arr(i-1,j,k)) 
                            - dt3*scal_arr(i-1,j,k,scomp-1)*divu_arr(i-1,j,k) 
                            + (dt3/hz)*scal_arr(i-1,j,k,scomp-1)*
                            (wmac_arr(i-1,j,k+1)-wmac_arr(i-1,j,k));
                        srxz = srx_arr(i,j,k) 
                            - (dt3/hz)*(simhz_arr(i,j,k+1)*wmac_arr(i,j,k+1)
                            - simhz_arr(i,j,k)*wmac_arr(i,j,k)) 
                            - dt3*scal_arr(i,j,k,scomp-1)*divu_arr(i,j,k) 
                            + (dt3/hz)*scal_arr(i,j,k,scomp-1)*
                            (wmac_arr(i,j,k+1)-wmac_arr(i,j,k));
                    } else {
                        // make slxz, srxz by updating 1D extrapolation
                        slxz = slx_arr(i,j,k) 
                            - (dt6/hz)*(wmac_arr(i-1,j,k+1)+wmac_arr(i-1,j,k)) 
                            *(simhz_arr(i-1,j,k+1)-simhz_arr(i-1,j,k));
                        srxz = srx_arr(i,j,k) 
                            - (dt6/hz)*(wmac_arr(i,j,k+1)+wmac_arr(i,j,k)) 
                            *(simhz_arr(i,j,k+1)-simhz_arr(i,j,k));
                    }

                    // impose lo side bc's
                    if (i == ilo) {
                        if (bclo == EXT_DIR) {
                            slxz = scal_arr(i-1,j,k,scomp-1);
                            srxz = scal_arr(i-1,j,k,scomp-1);
                        } else if (bclo ==  FOEXTRAP || bclo == HOEXTRAP) {
                            if (is_vel == 1 && scomp-1 == 0) {
                                srxz = min(srxz,0.0);
                            }
                            slxz = srxz;
                        } else if (bclo == REFLECT_EVEN) {
                            slxz = srxz;
                        } else if (bclo ==  REFLECT_ODD) {
                            slxz = 0.0;
                            srxz = 0.0;
                        }

                    // impose hi side bc's
                    } else if (i == ihi+1) {
                        if (bchi == EXT_DIR) {
                            slxz = scal_arr(i,j,k,scomp-1);
                            srxz = scal_arr(i,j,k,scomp-1);
                        } else if (bchi == FOEXTRAP || bchi == HOEXTRAP) {
                            if (is_vel == 1 && scomp-1 == 0) {
                                slxz = max(slxz,0.0);
                            }
                            srxz = slxz;
                        } else if (bchi == REFLECT_EVEN) {
                            srxz = slxz;
                        } else if (bchi == REFLECT_ODD) {
                            slxz = 0.0;
                            srxz = 0.0;
                        }
                    }

                    // make simhxy by solving Riemann problem
                    simhxz_arr(i,j,k) = (umac_arr(i,j,k) > 0.0) ?
                        slxz : srxz;
                    simhxz_arr(i,j,k) = (fabs(umac_arr(i,j,k)) > rel_eps) ?
                        simhxz_arr(i,j,k) : 0.5 * (slxz + srxz);
     
                });

                // simhyx
                imhbox = mfi.grownnodaltilebox(1, amrex::IntVect(0,0,1));
                // imhbox = amrex::grow(mfi.tilebox(), 2, 1);
                // imhbox = amrex::growHi(imhbox, 1, 1);

                Array4<Real> const simhyx_arr = simhyx.array(mfi);
                
                bclo = bcs[bccomp-1].lo()[1];
                bchi = bcs[bccomp-1].hi()[1];

                AMREX_PARALLEL_FOR_3D(imhbox, i, j, k, 
                {
                    Real slyx = 0.0;
                    Real sryx = 0.0;
                    // loop over appropriate yx faces
                    if (is_conservative == 1) {
                        // make slyx, sryx by updating 1D extrapolation
                        slyx = sly_arr(i,j,k) 
                            - (dt3/hx) * (simhx_arr(i+1,j-1,k)*umac_arr(i+1,j-1,k) 
                            - simhx_arr(i,j-1,k)*umac_arr(i,j-1,k)) 
                            - dt3*scal_arr(i,j-1,k,scomp-1)*divu_arr(i,j-1,k) 
                            + (dt3/hx)*scal_arr(i,j-1,k,scomp-1)*
                            (umac_arr(i+1,j-1,k)-umac_arr(i,j-1,k));
                        sryx = sry_arr(i,j,k) 
                            - (dt3/hx)*(simhx_arr(i+1,j,k)*umac_arr(i+1,j,k)
                            - simhx_arr(i,j,k)*umac_arr(i,j,k)) 
                            - dt3*scal_arr(i,j,k,scomp-1)*divu_arr(i,j,k) 
                            + (dt3/hx)*scal_arr(i,j,k,scomp-1)*
                            (umac_arr(i+1,j,k)-umac_arr(i,j,k));
                    } else {
                        // make slyx, sryx by updating 1D extrapolation
                        slyx = sly_arr(i,j,k) 
                            - (dt6/hx)*(umac_arr(i+1,j-1,k)+umac_arr(i,j-1,k)) 
                            *(simhx_arr(i+1,j-1,k)-simhx_arr(i,j-1,k));
                        sryx = sry_arr(i,j,k) 
                            - (dt6/hx)*(umac_arr(i+1,j,k)+umac_arr(i,j,k)) 
                            *(simhx_arr(i+1,j,k)-simhx_arr(i,j,k));
                    }

                    // impose lo side bc's
                    if (j == jlo) {
                        if (bclo == EXT_DIR) {
                            slyx = scal_arr(i,j-1,k,scomp-1);
                            sryx = scal_arr(i,j-1,k,scomp-1);
                        } else if (bclo ==  FOEXTRAP || bclo == HOEXTRAP) {
                            if (is_vel == 1 && scomp-1 == 1) {
                                sryx = min(sryx,0.0);
                            }
                            slyx = sryx;
                        } else if (bclo == REFLECT_EVEN) {
                            slyx = sryx;
                        } else if (bclo ==  REFLECT_ODD) {
                            slyx = 0.0;
                            sryx = 0.0;
                        }

                    // impose hi side bc's
                    } else if (j == jhi+1) {
                        if (bchi == EXT_DIR) {
                            slyx = scal_arr(i,j,k,scomp-1);
                            sryx = scal_arr(i,j,k,scomp-1);
                        } else if (bchi == FOEXTRAP || bchi == HOEXTRAP) {
                            if (is_vel == 1 && scomp-1 == 1) {
                                slyx = max(slyx,0.0);
                            }
                            sryx = slyx;
                        } else if (bchi == REFLECT_EVEN) {
                            sryx = slyx;
                        } else if (bchi == REFLECT_ODD) {
                            slyx = 0.0;
                            sryx = 0.0;
                        }
                    }

                    // make simhxy by solving Riemann problem
                    simhyx_arr(i,j,k) = (vmac_arr(i,j,k) > 0.0) ?
                        slyx : sryx;
                    simhyx_arr(i,j,k) = (fabs(vmac_arr(i,j,k)) > rel_eps) ?
                        simhyx_arr(i,j,k) : 0.5 * (slyx + sryx);
     
                });

                // simhyz
                imhbox = mfi.grownnodaltilebox(1, amrex::IntVect(1,0,0)); 
                // imhbox = amrex::grow(mfi.tilebox(), 0, 1);
                // imhbox = amrex::growHi(imhbox, 1, 1);

                Array4<Real> const simhyz_arr = simhyz.array(mfi);

                AMREX_PARALLEL_FOR_3D(imhbox, i, j, k, 
                {
                    Real slyz = 0.0;
                    Real sryz = 0.0;
                    // loop over appropriate yz faces
                    if (is_conservative == 1) {
                        // make slyz, sryz by updating 1D extrapolation
                        slyz = sly_arr(i,j,k) 
                            - (dt3/hz) * (simhz_arr(i,j-1,k+1)*wmac_arr(i,j-1,k+1) 
                            - simhz_arr(i,j-1,k)*umac_arr(i,j-1,k)) 
                            - dt3*scal_arr(i,j-1,k,scomp-1)*divu_arr(i,j-1,k) 
                            + (dt3/hz)*scal_arr(i,j-1,k,scomp-1)*
                            (wmac_arr(i,j-1,k+1)-wmac_arr(i,j-1,k));
                        sryz = sry_arr(i,j,k) 
                            - (dt3/hz)*(simhz_arr(i,j,k+1)*wmac_arr(i,j,k+1)
                            - simhz_arr(i,j,k)*wmac_arr(i,j,k)) 
                            - dt3*scal_arr(i,j,k,scomp-1)*divu_arr(i,j,k) 
                            + (dt3/hz)*scal_arr(i,j,k,scomp-1)*
                            (wmac_arr(i,j,k+1)-wmac_arr(i,j,k));
                    } else {
                        // make slyz, sryz by updating 1D extrapolation
                        slyz = sly_arr(i,j,k) 
                            - (dt6/hz)*(wmac_arr(i,j-1,k+1)+wmac_arr(i,j-1,k)) 
                            *(simhz_arr(i,j-1,k+1)-simhz_arr(i,j-1,k));
                        sryz = sry_arr(i,j,k) 
                            - (dt6/hz)*(wmac_arr(i,j,k+1)+wmac_arr(i,j,k)) 
                            *(simhz_arr(i,j,k+1)-simhz_arr(i,j,k));
                    }

                    // impose lo side bc's
                    if (j == jlo) {
                        if (bclo == EXT_DIR) {
                            slyz = scal_arr(i,j-1,k,scomp-1);
                            sryz = scal_arr(i,j-1,k,scomp-1);
                        } else if (bclo ==  FOEXTRAP || bclo == HOEXTRAP) {
                            if (is_vel == 1 && scomp-1 == 1) {
                                sryz = min(sryz,0.0);
                            }
                            slyz = sryz;
                        } else if (bclo == REFLECT_EVEN) {
                            slyz = sryz;
                        } else if (bclo ==  REFLECT_ODD) {
                            slyz = 0.0;
                            sryz = 0.0;
                        }

                    // impose hi side bc's
                    } else if (j == jhi+1) {
                        if (bchi == EXT_DIR) {
                            slyz = scal_arr(i,j,k,scomp-1);
                            sryz = scal_arr(i,j,k,scomp-1);
                        } else if (bchi == FOEXTRAP || bchi == HOEXTRAP) {
                            if (is_vel == 1 && scomp-1 == 1) {
                                slyz = max(slyz,0.0);
                            }
                            sryz = slyz;
                        } else if (bchi == REFLECT_EVEN) {
                            sryz = slyz;
                        } else if (bchi == REFLECT_ODD) {
                            slyz = 0.0;
                            sryz = 0.0;
                        }
                    }

                    // make simhyz by solving Riemann problem
                    simhyz_arr(i,j,k) = (vmac_arr(i,j,k) > 0.0) ?
                        slyz : sryz;
                    simhyz_arr(i,j,k) = (fabs(vmac_arr(i,j,k)) > rel_eps) ?
                        simhyz_arr(i,j,k) : 0.5 * (slyz + sryz);
     
                });

                // simhzx
                imhbox = mfi.grownnodaltilebox(2, amrex::IntVect(0,1,0));
                // imhbox = amrex::grow(mfi.tilebox(), 1, 1);
                // imhbox = amrex::growHi(imhbox, 2, 1);

                Array4<Real> const simhzx_arr = simhzx.array(mfi);

                bclo = bcs[bccomp-1].lo()[2];
                bchi = bcs[bccomp-1].hi()[2];

                AMREX_PARALLEL_FOR_3D(imhbox, i, j, k, 
                {
                    Real slzx = 0.0;
                    Real srzx = 0.0;
                    // loop over appropriate zx faces
                    if (is_conservative == 1) {
                        // make slzx, srzx by updating 1D extrapolation
                        slzx = slz_arr(i,j,k) 
                            - (dt3/hx) * (simhx_arr(i+1,j,k-1)*umac_arr(i+1,j,k-1) 
                            - simhx_arr(i,j,k-1)*umac_arr(i,j,k-1)) 
                            - dt3*scal_arr(i,j,k-1,scomp-1)*divu_arr(i,j,k-1) 
                            + (dt3/hx)*scal_arr(i,j,k-1,scomp-1)*
                            (umac_arr(i+1,j,k-1)-umac_arr(i,j,k-1));
                        srzx = srz_arr(i,j,k) 
                            - (dt3/hx)*(simhx_arr(i+1,j,k)*umac_arr(i+1,j,k)
                            - simhx_arr(i,j,k)*umac_arr(i,j,k)) 
                            - dt3*scal_arr(i,j,k,scomp-1)*divu_arr(i,j,k) 
                            + (dt3/hx)*scal_arr(i,j,k,scomp-1)*
                            (umac_arr(i+1,j,k)-umac_arr(i,j,k));
                    } else {
                        // make slzx, srzx by updating 1D extrapolation
                        slzx = slz_arr(i,j,k) 
                            - (dt6/hx)*(umac_arr(i+1,j,k-1)+umac_arr(i,j,k-1)) 
                            *(simhx_arr(i+1,j,k-1)-simhx_arr(i,j,k-1));
                        srzx = srz_arr(i,j,k) 
                            - (dt6/hx)*(umac_arr(i+1,j,k)+umac_arr(i,j,k)) 
                            *(simhx_arr(i+1,j,k)-simhx_arr(i,j,k));
                    }

                    // impose lo side bc's
                    if (k == klo) {
                        if (bclo == EXT_DIR) {
                            slzx = scal_arr(i,j,k-1,scomp-1);
                            srzx = scal_arr(i,j,k-1,scomp-1);
                        } else if (bclo ==  FOEXTRAP || bclo == HOEXTRAP) {
                            if (is_vel == 1 && scomp-1 == 2) {
                                srzx = min(srzx,0.0);
                            }
                            slzx = srzx;
                        } else if (bclo == REFLECT_EVEN) {
                            slzx = srzx;
                        } else if (bclo ==  REFLECT_ODD) {
                            slzx = 0.0;
                            srzx = 0.0;
                        }

                    // impose hi side bc's
                    } else if (k == khi+1) {
                        if (bchi == EXT_DIR) {
                            slzx = scal_arr(i,j,k,scomp-1);
                            srzx = scal_arr(i,j,k,scomp-1);
                        } else if (bchi == FOEXTRAP || bchi == HOEXTRAP) {
                            if (is_vel == 1 && scomp-1 == 2) {
                                slzx = max(slzx,0.0);
                            }
                            srzx = slzx;
                        } else if (bchi == REFLECT_EVEN) {
                            srzx = slzx;
                        } else if (bchi == REFLECT_ODD) {
                            slzx = 0.0;
                            srzx = 0.0;
                        }
                    }

                    // make simhzx by solving Riemann problem
                    simhzx_arr(i,j,k) = (wmac_arr(i,j,k) > 0.0) ?
                        slzx : srzx;
                    simhzx_arr(i,j,k) = (fabs(wmac_arr(i,j,k)) > rel_eps) ?
                        simhzx_arr(i,j,k) : 0.5 * (slzx + srzx);
     
                });

                // simhzy
                imhbox = mfi.grownnodaltilebox(2, IntVect(1,0,0));
                // imhbox = amrex::grow(mfi.tilebox(), 0, 1);
                // imhbox = amrex::growHi(imhbox, 2, 1);

                Array4<Real> const simhzy_arr = simhzy.array(mfi);

                AMREX_PARALLEL_FOR_3D(imhbox, i, j, k, 
                {
                    Real slzy = 0.0;
                    Real srzy = 0.0;
                    // loop over appropriate zy faces
                    if (is_conservative == 1) {
                        // make slzy, srzy by updating 1D extrapolation
                        slzy = slz_arr(i,j,k) 
                            - (dt3/hy) * (simhy_arr(i,j+1,k-1)*vmac_arr(i,j+1,k-1) 
                            - simhy_arr(i,j,k-1)*vmac_arr(i,j,k-1)) 
                            - dt3*scal_arr(i,j,k-1,scomp-1)*divu_arr(i,j,k-1) 
                            + (dt3/hy)*scal_arr(i,j,k-1,scomp-1)*
                            (umac_arr(i,j+1,k-1)-vmac_arr(i,j,k-1));
                        srzy = srz_arr(i,j,k) 
                            - (dt3/hy)*(simhy_arr(i,j+1,k)*vmac_arr(i,j+1,k)
                            - simhy_arr(i,j,k)*vmac_arr(i,j,k)) 
                            - dt3*scal_arr(i,j,k,scomp-1)*divu_arr(i,j,k) 
                            + (dt3/hy)*scal_arr(i,j,k,scomp-1)*
                            (vmac_arr(i,j+1,k)-vmac_arr(i,j,k));
                    } else {
                        // make slzy, srzy by updating 1D extrapolation
                        slzy = slz_arr(i,j,k) 
                            - (dt6/hy)*(vmac_arr(i,j+1,k-1)+vmac_arr(i,j,k-1)) 
                            *(simhy_arr(i,j+1,k-1)-simhy_arr(i,j,k-1));
                        srzy = srz_arr(i,j,k) 
                            - (dt6/hy)*(vmac_arr(i,j+1,k)+vmac_arr(i,j,k)) 
                            *(simhy_arr(i,j+1,k)-simhy_arr(i,j,k));
                    }

                    // impose lo side bc's
                    if (k == klo) {
                        if (bclo == EXT_DIR) {
                            slzy = scal_arr(i,j,k-1,scomp-1);
                            srzy = scal_arr(i,j,k-1,scomp-1);
                        } else if (bclo ==  FOEXTRAP || bclo == HOEXTRAP) {
                            if (is_vel == 1 && scomp-1 == 2) {
                                srzy = min(srzy,0.0);
                            }
                            slzy = srzy;
                        } else if (bclo == REFLECT_EVEN) {
                            slzy = srzy;
                        } else if (bclo ==  REFLECT_ODD) {
                            slzy = 0.0;
                            srzy = 0.0;
                        }

                    // impose hi side bc's
                    } else if (k == khi+1) {
                        if (bchi == EXT_DIR) {
                            slzy = scal_arr(i,j,k,scomp-1);
                            srzy = scal_arr(i,j,k,scomp-1);
                        } else if (bchi == FOEXTRAP || bchi == HOEXTRAP) {
                            if (is_vel == 1 && scomp-1 == 2) {
                                slzy = max(slzy,0.0);
                            }
                            srzy = slzy;
                        } else if (bchi == REFLECT_EVEN) {
                            srzy = slzy;
                        } else if (bchi == REFLECT_ODD) {
                            slzy = 0.0;
                            srzy = 0.0;
                        }
                    }

                    // make simhzy by solving Riemann problem
                    simhzy_arr(i,j,k) = (wmac_arr(i,j,k) > 0.0) ?
                        slzy : srzy;
                    simhzy_arr(i,j,k) = (fabs(wmac_arr(i,j,k)) > rel_eps) ?
                        simhzy_arr(i,j,k) : 0.5 * (slzy + srzy);
     
                });

                ///////////////////////////////////////////////
                // Create sedgelx, etc.
                ///////////////////////////////////////////////

                Array4<Real> const Ipf_arr = Ipf.array(mfi);
                Array4<Real> const Imf_arr = Imf.array(mfi);
                Array4<Real> const force_arr = force[lev].array(mfi);

                Array4<Real> const sedgex_arr = sedge[lev][0].array(mfi);
                Array4<Real> const sedgey_arr = sedge[lev][1].array(mfi);
                Array4<Real> const sedgez_arr = sedge[lev][2].array(mfi);

                int ppm_trace_forces_local = ppm_trace_forces;

                // x-direction
                bclo = bcs[bccomp-1].lo()[0];
                bchi = bcs[bccomp-1].hi()[0];
                AMREX_PARALLEL_FOR_3D(xbx, i, j, k, 
                {
                    Real sedgelx = 0.0;
                    Real sedgerx = 0.0;

                    Real fl = (ppm_trace_forces_local == 0) ? 
                        force_arr(i-1,j,k,scomp-1) : Ipf_arr(i-1,j,k,0);
                    Real fr = (ppm_trace_forces_local == 0) ? 
                        force_arr(i,j,k,scomp-1) : Imf_arr(i,j,k,0);

                    // make sedgelx, sedgerx
                    if (is_conservative == 1) {
                        sedgelx = slx_arr(i,j,k) 
                            - (dt2/hy)*(simhyz_arr(i-1,j+1,k  )*vmac_arr(i-1,j+1,k) 
                            - simhyz_arr(i-1,j,k)*vmac_arr(i-1,j,k)) 
                            - (dt2/hz)*(simhzy_arr(i-1,j,k+1)*wmac_arr(i-1,j,k+1) 
                            - simhzy_arr(i-1,j,k)*wmac_arr(i-1,j,k)) 
                            - (dt2/hx)*scal_arr(i-1,j,k,scomp-1)*(umac_arr(i,j,k)-umac_arr(i-1,j,k)) 
                            + dt2*fl;

                        sedgerx = srx_arr(i,j,k) 
                            - (dt2/hy)*(simhyz_arr(i,j+1,k)*vmac_arr(i,j+1,k) 
                            - simhyz_arr(i,j,k)*vmac_arr(i,j,k)) 
                            - (dt2/hz)*(simhzy_arr(i,j,k+1)*wmac_arr(i,j,k+1) 
                            - simhzy_arr(i,j,k)*wmac_arr(i,j,k)) 
                            - (dt2/hx)*scal_arr(i,j,k,scomp-1)*(umac_arr(i+1,j,k)-umac_arr(i,j,k)) 
                            + dt2*fr;
                    } else {
                        sedgelx = slx_arr(i,j,k) 
                            - (dt4/hy)*(vmac_arr(i-1,j+1,k)+vmac_arr(i-1,j,k))* 
                            (simhyz_arr(i-1,j+1,k)-simhyz_arr(i-1,j,k)) 
                            - (dt4/hz)*(wmac_arr(i-1,j,k+1)+wmac_arr(i-1,j,k))* 
                            (simhzy_arr(i-1,j,k+1)-simhzy_arr(i-1,j,k)) 
                            + dt2*fl;

                        sedgerx = srx_arr(i,j,k) 
                            - (dt4/hy)*(vmac_arr(i,j+1,k)+vmac_arr(i,j,k))* 
                            (simhyz_arr(i,j+1,k)-simhyz_arr(i,j,k)) 
                            - (dt4/hz)*(wmac_arr(i,j,k+1)+wmac_arr(i,j,k))* 
                            (simhzy_arr(i,j,k+1)-simhzy_arr(i,j,k)) 
                            + dt2*fr;
                    } 

                    // make sedgex by solving Riemann problem
                    // boundary conditions enforced outside of i,j,k loop
                    sedgex_arr(i,j,k,scomp-1) = (umac_arr(i,j,k) > 0.0) ? 
                        sedgelx : sedgerx;
                    sedgex_arr(i,j,k,scomp-1) = (fabs(umac_arr(i,j,k))  > rel_eps) ? 
                        sedgex_arr(i,j,k,scomp-1) : 0.5*(sedgelx+sedgerx);

                    // impose lo side bc's
                    if (i == ilo) {
                        if (bclo == EXT_DIR) {
                            sedgex_arr(i,j,k,scomp-1) = scal_arr(i-1,j,k,scomp-1);
                        } else if (bclo == FOEXTRAP || bclo == HOEXTRAP) {
                            if (is_vel == 1 && scomp-1 == 0) {
                                sedgex_arr(i,j,k,scomp-1) = min(sedgerx,0.0);
                            } else {
                                sedgex_arr(i,j,k,scomp-1) = sedgerx;
                            } 
                        } else if (bclo == REFLECT_EVEN) {
                            sedgex_arr(i,j,k,scomp-1) = sedgerx;
                        } else if (bclo == REFLECT_ODD) {
                            sedgex_arr(i,j,k,scomp-1) = 0.0;
                        }

                    // impose hi side bc's
                    } else if (i == ihi+1) {
                        if (bchi == EXT_DIR) {
                            sedgex_arr(i,j,k,scomp-1) = scal_arr(i,j,k,scomp-1);
                        } else if (bchi == FOEXTRAP || bchi == HOEXTRAP) {
                            if (is_vel == 1 && scomp-1 == 0) {
                                sedgex_arr(i,j,k,scomp-1) = max(sedgelx,0.0);
                            } else {
                                sedgex_arr(i,j,k,scomp-1) = sedgelx;
                            } 
                        } else if (bchi == REFLECT_EVEN) {
                            sedgex_arr(i,j,k,scomp-1) = sedgelx;
                        } else if (bchi == REFLECT_ODD) {
                            sedgex_arr(i,j,k,scomp-1) = 0.0;
                        } 
                    } 
     
                });

                // y-direction
                bclo = bcs[bccomp-1].lo()[1];
                bchi = bcs[bccomp-1].hi()[1];
                AMREX_PARALLEL_FOR_3D(ybx, i, j, k, 
                {
                    Real sedgely = 0.0;
                    Real sedgery = 0.0;

                    Real fl = (ppm_trace_forces_local == 0) ? 
                        force_arr(i,j-1,k,scomp-1) : Ipf_arr(i,j-1,k,1);
                    Real fr = (ppm_trace_forces_local == 0) ? 
                        force_arr(i,j,k,scomp-1) : Imf_arr(i,j,k,1);

                    // make sedgely, sedgery
                    if (is_conservative == 1) {
                        sedgely = sly_arr(i,j,k) 
                            - (dt2/hx)*(simhxz_arr(i+1,j-1,k  )*umac_arr(i+1,j-1,k) 
                            - simhxz_arr(i,j-1,k)*umac_arr(i,j-1,k)) 
                            - (dt2/hz)*(simhzx_arr(i,j-1,k+1)*wmac_arr(i,j-1,k+1) 
                            - simhzx_arr(i,j-1,k)*wmac_arr(i,j-1,k)) 
                            - (dt2/hy)*scal_arr(i,j-1,k,scomp-1)*(vmac_arr(i,j,k)-vmac_arr(i,j-1,k)) 
                            + dt2*fl;

                        sedgery = sry_arr(i,j,k) 
                            - (dt2/hx)*(simhxz_arr(i+1,j,k)*umac_arr(i+1,j,k) 
                            - simhxz_arr(i,j,k)*umac_arr(i,j,k)) 
                            - (dt2/hz)*(simhzx_arr(i,j,k+1)*wmac_arr(i,j,k+1) 
                            - simhzx_arr(i,j,k)*wmac_arr(i,j,k)) 
                            - (dt2/hy)*scal_arr(i,j,k,scomp-1)*(vmac_arr(i,j+1,k)-vmac_arr(i,j,k)) 
                            + dt2*fr;
                    } else {
                        sedgely = sly_arr(i,j,k) 
                            - (dt4/hx)*(umac_arr(i+1,j-1,k)+umac_arr(i,j-1,k))* 
                            (simhxz_arr(i+1,j-1,k)-simhxz_arr(i,j-1,k)) 
                            - (dt4/hz)*(wmac_arr(i,j-1,k+1)+wmac_arr(i,j-1,k))* 
                            (simhzx_arr(i,j-1,k+1)-simhzx_arr(i,j-1,k)) 
                            + dt2*fl;

                        sedgery = sry_arr(i,j,k) 
                            - (dt4/hx)*(umac_arr(i+1,j,k)+umac_arr(i,j,k))* 
                            (simhxz_arr(i+1,j,k)-simhxz_arr(i,j,k)) 
                            - (dt4/hz)*(wmac_arr(i,j,k+1)+wmac_arr(i,j,k))* 
                            (simhzx_arr(i,j,k+1)-simhzx_arr(i,j,k)) 
                            + dt2*fr;
                    } 

                    // make sedgey by solving Riemann problem
                    // boundary conditions enforced outside of i,j,k loop
                    sedgey_arr(i,j,k,scomp-1) = (vmac_arr(i,j,k) > 0.0) ? 
                        sedgely : sedgery;
                    sedgey_arr(i,j,k,scomp-1) = (fabs(vmac_arr(i,j,k))  > rel_eps) ? 
                        sedgey_arr(i,j,k,scomp-1) : 0.5*(sedgely+sedgery);

                    // impose lo side bc's
                    if (j == jlo) {
                        if (bclo == EXT_DIR) {
                            sedgey_arr(i,j,k,scomp-1) = scal_arr(i,j-1,k,scomp-1);
                        } else if (bclo == FOEXTRAP || bclo == HOEXTRAP) {
                            if (is_vel == 1 && scomp-1 == 1) {
                                sedgey_arr(i,j,k,scomp-1) = min(sedgery,0.0);
                            } else {
                                sedgey_arr(i,j,k,scomp-1) = sedgery;
                            } 
                        } else if (bclo == REFLECT_EVEN) {
                            sedgey_arr(i,j,k,scomp-1) = sedgery;
                        } else if (bclo == REFLECT_ODD) {
                            sedgey_arr(i,j,k,scomp-1) = 0.0;
                        }

                    // impose hi side bc's
                    } else if (j == jhi+1) {
                        if (bchi == EXT_DIR) {
                            sedgey_arr(i,j,k,scomp-1) = scal_arr(i,j,k,scomp-1);
                        } else if (bchi == FOEXTRAP || bchi == HOEXTRAP) {
                            if (is_vel == 1 && scomp-1 == 1) {
                                sedgey_arr(i,j,k,scomp-1) = max(sedgely,0.0);
                            } else {
                                sedgey_arr(i,j,k,scomp-1) = sedgely;
                            } 
                        } else if (bchi == REFLECT_EVEN) {
                            sedgey_arr(i,j,k,scomp-1) = sedgely;
                        } else if (bchi == REFLECT_ODD) {
                            sedgey_arr(i,j,k,scomp-1) = 0.0;
                        } 
                    } 
     
                });

                // z-direction
                bclo = bcs[bccomp-1].lo()[2];
                bchi = bcs[bccomp-1].hi()[2];
                AMREX_PARALLEL_FOR_3D(zbx, i, j, k, 
                {
                    Real sedgelz = 0.0;
                    Real sedgerz = 0.0;

                    Real fl = (ppm_trace_forces_local == 0) ? 
                        force_arr(i,j,k-1,scomp-1) : Ipf_arr(i,j,k-1,2);
                    Real fr = (ppm_trace_forces_local == 0) ? 
                        force_arr(i,j,k,scomp-1) : Imf_arr(i,j,k,2);

                    // make sedgelz, sedgerz
                    if (is_conservative == 1) {
                        sedgelz = slz_arr(i,j,k) 
                            - (dt2/hx)*(simhxy_arr(i+1,j,k-1)*umac_arr(i+1,j,k-1) 
                            - simhxy_arr(i,j,k-1)*umac_arr(i,j,k-1)) 
                            - (dt2/hy)*(simhyx_arr(i,j+1,k-1)*vmac_arr(i,j+1,k-1) 
                            - simhyx_arr(i,j,k-1)*vmac_arr(i,j,k-1)) 
                            - (dt2/hz)*scal_arr(i,j,k-1,scomp-1)*(wmac_arr(i,j,k)-wmac_arr(i,j,k-1)) 
                            + dt2*fl;

                        sedgerz = srz_arr(i,j,k) 
                            - (dt2/hx)*(simhxy_arr(i+1,j,k)*umac_arr(i+1,j,k) 
                            - simhxy_arr(i,j,k)*umac_arr(i,j,k)) 
                            - (dt2/hy)*(simhyx_arr(i,j+1,k)*vmac_arr(i,j+1,k) 
                            - simhyx_arr(i,j,k)*vmac_arr(i,j,k)) 
                            - (dt2/hz)*scal_arr(i,j,k,scomp-1)*(wmac_arr(i,j,k+1)-wmac_arr(i,j,k)) 
                            + dt2*fr;
                    } else {
                        sedgelz = slz_arr(i,j,k) 
                            - (dt4/hx)*(umac_arr(i+1,j,k-1)+umac_arr(i,j,k-1))* 
                            (simhxy_arr(i+1,j,k-1)-simhxy_arr(i,j,k-1)) 
                            - (dt4/hy)*(vmac_arr(i,j+1,k-1)+vmac_arr(i,j,k-1))* 
                            (simhyx_arr(i,j+1,k-1)-simhyx_arr(i,j,k-1)) 
                            + dt2*fl;

                        sedgerz = srz_arr(i,j,k) 
                            - (dt4/hx)*(umac_arr(i+1,j,k)+umac_arr(i,j,k))* 
                            (simhxy_arr(i+1,j,k)-simhxy_arr(i,j,k)) 
                            - (dt4/hy)*(vmac_arr(i,j+1,k)+vmac_arr(i,j,k))* 
                            (simhyx_arr(i,j+1,k)-simhyx_arr(i,j,k)) 
                            + dt2*fr;
                    } 

                    sl_arr(i,j,k) = sedgelz;
                    sr_arr(i,j,k) = sedgerz;

                    // make sedgez by solving Riemann problem
                    // boundary conditions enforced outside of i,j,k loop
                    sedgez_arr(i,j,k,scomp-1) = (wmac_arr(i,j,k) > 0.0) ? 
                        sedgelz : sedgerz;
                    sedgez_arr(i,j,k,scomp-1) = (fabs(wmac_arr(i,j,k))  > rel_eps) ? 
                        sedgez_arr(i,j,k,scomp-1) : 0.5*(sedgelz+sedgerz);

                    // impose lo side bc's
                    if (k == klo) {
                        if (bclo == EXT_DIR) {
                            sedgez_arr(i,j,k,scomp-1) = scal_arr(i,j,k-1,scomp-1);
                        } else if (bclo == FOEXTRAP || bclo == HOEXTRAP) {
                            if (is_vel == 1 && scomp-1 == 2) {
                                sedgez_arr(i,j,k,scomp-1) = min(sedgerz,0.0);
                            } else {
                                sedgez_arr(i,j,k,scomp-1) = sedgerz;
                            } 
                        } else if (bclo == REFLECT_EVEN) {
                            sedgez_arr(i,j,k,scomp-1) = sedgerz;
                        } else if (bclo == REFLECT_ODD) {
                            sedgez_arr(i,j,k,scomp-1) = 0.0;
                        }

                    // impose hi side bc's
                    } else if (k == khi+1) {
                        if (bchi == EXT_DIR) {
                            sedgez_arr(i,j,k,scomp-1) = scal_arr(i,j,k,scomp-1);
                        } else if (bchi == FOEXTRAP || bchi == HOEXTRAP) {
                            if (is_vel == 1 && scomp-1 == 2) {
                                sedgez_arr(i,j,k,scomp-1) = max(sedgelz,0.0);
                            } else {
                                sedgez_arr(i,j,k,scomp-1) = sedgelz;
                            } 
                        } else if (bchi == REFLECT_EVEN) {
                            sedgez_arr(i,j,k,scomp-1) = sedgelz;
                        } else if (bchi == REFLECT_ODD) {
                            sedgez_arr(i,j,k,scomp-1) = 0.0;
                        } 
                    } 
     
                });


            } // end MFIter loop
        } // end loop over components

#ifdef AMREX_USE_CUDA
        clean_bc(bc_f);
#endif

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

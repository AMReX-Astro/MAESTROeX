
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
                       AMREX_INT_ANYD(domainBox.loVect()),
                       AMREX_INT_ANYD(domainBox.hiVect()),
                       bc_f, AMREX_REAL_ANYD(dx), dt, false,
                       1,1,AMREX_SPACEDIM);
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
                       AMREX_INT_ANYD(domainBox.loVect()),
                       AMREX_INT_ANYD(domainBox.hiVect()),
                       bc_f, AMREX_REAL_ANYD(dx), dt, false,
                       2,2,AMREX_SPACEDIM);
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
                        AMREX_INT_ANYD(domainBox.loVect()),
                        AMREX_INT_ANYD(domainBox.hiVect()),
                        bc_f, AMREX_REAL_ANYD(dx), dt, false,
                        3,3,AMREX_SPACEDIM);
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
                       AMREX_INT_ANYD(domainBox.loVect()),
                       AMREX_INT_ANYD(domainBox.hiVect()),
                       bc_f, AMREX_REAL_ANYD(dx), dt, false,
                       1,1,AMREX_SPACEDIM);

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
                           AMREX_INT_ANYD(domainBox.loVect()),
                           AMREX_INT_ANYD(domainBox.hiVect()),
                           bc_f, AMREX_REAL_ANYD(dx), dt, false,
                           1,1,AMREX_SPACEDIM);
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
                      AMREX_INT_ANYD(domainBox.loVect()),
                      AMREX_INT_ANYD(domainBox.hiVect()),
                      bc_f, AMREX_REAL_ANYD(dx), dt, false,
                      2,2,AMREX_SPACEDIM);

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
                    AMREX_INT_ANYD(domainBox.loVect()),
                    AMREX_INT_ANYD(domainBox.hiVect()),
                    bc_f, AMREX_REAL_ANYD(dx), dt, false,
                    2,2,AMREX_SPACEDIM);
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
                      AMREX_INT_ANYD(domainBox.loVect()),
                      AMREX_INT_ANYD(domainBox.hiVect()),
                      bc_f, AMREX_REAL_ANYD(dx), dt, false,
                      3,3,AMREX_SPACEDIM);

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
                    AMREX_INT_ANYD(domainBox.loVect()),
                    AMREX_INT_ANYD(domainBox.hiVect()),
                    bc_f, AMREX_REAL_ANYD(dx), dt, false,
                    3,3,AMREX_SPACEDIM);
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
                       const Vector<MultiFab>& force,
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

#pragma gpu box(obx)
                make_divu(AMREX_INT_ANYD(obx.loVect()),
                          AMREX_INT_ANYD(obx.hiVect()),
                          BL_TO_FORTRAN_ANYD(divu[mfi]),
                          BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
                          BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
                          BL_TO_FORTRAN_ANYD(wmac_mf[mfi]),
                          AMREX_REAL_ANYD(dx), is_conservative);
                          
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
                           AMREX_INT_ANYD(domainBox.loVect()),
                           AMREX_INT_ANYD(domainBox.hiVect()),
                           bc_f, AMREX_REAL_ANYD(dx), dt, true,
                           scomp, bccomp, nbccomp);

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
                               AMREX_INT_ANYD(domainBox.loVect()),
                               AMREX_INT_ANYD(domainBox.hiVect()),
                               bc_f, AMREX_REAL_ANYD(dx), dt, true,
                               scomp, bccomp, nbccomp);
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

                Real ppm_type_local = ppm_type;

                // loop over appropriate x-faces
                AMREX_PARALLEL_FOR_3D(mxbx, i, j, k, 
                {
                    if (ppm_type_local == 0) {
                        slx_arr(i,j,k) = scal_arr(i-1,j,k,scomp-1) + 0.5 * (1.0 - dt * umac_arr(i,j,k) / dx[0]) * Ip_arr(i-1,j,k,0);
                        srx_arr(i,j,k) = scal_arr(i,j,k,scomp-1) - 0.5 * (1.0 - dt * umac_arr(i,j,k) / dx[0]) * Ip_arr(i,j,k,0);
                    } else if (ppm_type_local == 1 || ppm_type_local == 2) {
                        slx_arr(i,j,k) = Ip_arr(i-1,j,k,0);
                        srx_arr(i,j,k) = Im_arr(i,j,k,0);
                    }
                });

                // loop over appropriate y-faces
                AMREX_PARALLEL_FOR_3D(mybx, i, j, k, 
                {
                    if (ppm_type_local == 0) {
                        sly_arr(i,j,k) = scal_arr(i,j-1,k,scomp-1) + 0.5 * (1.0 - dt * vmac_arr(i,j,k) / dx[1]) * Ip_arr(i,j-1,k,0);
                        sry_arr(i,j,k) = scal_arr(i,j,k,scomp-1) - 0.5 * (1.0 - dt * vmac_arr(i,j,k) / dx[1]) * Ip_arr(i,j,k,0);
                    } else if (ppm_type_local == 1 || ppm_type_local == 2) {
                        sly_arr(i,j,k) = Ip_arr(i,j-1,k,1);
                        sry_arr(i,j,k) = Im_arr(i,j,k,1);
                    }
                });

                // loop over appropriate z-faces
                AMREX_PARALLEL_FOR_3D(mzbx, i, j, k, 
                {
                    if (ppm_type_local == 0) {
                        slz_arr(i,j,k) = scal_arr(i,j,k-1,scomp-1) + 0.5 * (1.0 - dt * wmac_arr(i,j,k) / dx[2]) * Ip_arr(i,j,k-1,0);
                        srz_arr(i,j,k) = scal_arr(i,j,k,scomp-1) - 0.5 * (1.0 - dt * wmac_arr(i,j,k) / dx[2]) * Ip_arr(i,j,k,0);
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
                            slx_arr(i,j,k) = scal_arr(i,j-1,k,scomp-1);
                            srx_arr(i,j,k) = scal_arr(i,j-1,k,scomp-1);
                        } else if (bclo == FOEXTRAP || bclo == HOEXTRAP) {
                            if (is_vel == 1 && scomp-1 == 1) {
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
                    } else if (j == jhi+1) {
                        if (bchi == EXT_DIR) {
                            slx_arr(i,j,k) = scal_arr(i,j,k,scomp-1);
                            srx_arr(i,j,k) = scal_arr(i,j,k,scomp-1);
                        } else if (bchi == FOEXTRAP || bchi == HOEXTRAP) {
                            if (is_vel == 1 && scomp-1 == 1) {
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
                
                int klo = domainBox.loVect()[2];
                int khi = domainBox.hiVect()[2];
                bclo = bcs[bccomp-1].lo()[2];
                bchi = bcs[bccomp-1].hi()[2];
                AMREX_PARALLEL_FOR_3D(mzbx, i, j, k, 
                {
                    // impose lo side bc's
                    if (k == klo) {
                        // if (bclo == EXT_DIR) {
                        //     slx_arr(i,j,k) = scal_arr(i,j,k-1,scomp-1);
                        //     srx_arr(i,j,k) = scal_arr(i,j,k-1,scomp-1);
                        // } else if (bclo == FOEXTRAP || bclo == HOEXTRAP) {
                        //     if (is_vel == 1 && scomp-1 == 2) {
                        //         srx_arr(i,j,k) = min(srx_arr(i,j,k), 0.0);
                        //     }
                        //     slx_arr(i,j,k) = srx_arr(i,j,k);
                        // } else if (bclo == REFLECT_EVEN) {
                        //     slx_arr(i,j,k) = srx_arr(i,j,k);
                        // } else if (bclo == REFLECT_ODD) {
                        //     slx_arr(i,j,k) = 0.0;
                        //     srx_arr(i,j,k) = 0.0;
                        // }
                    // impose hi side bc's
                    } else if (k == khi+1) {
                        if (bchi == EXT_DIR) {
                            slx_arr(i,j,k) = scal_arr(i,j,k,scomp-1);
                            srx_arr(i,j,k) = scal_arr(i,j,k,scomp-1);
                        } else if (bchi == FOEXTRAP || bchi == HOEXTRAP) {
                            if (is_vel == 1 && scomp-1 == 2) {
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

                // x-direction
// #pragma gpu box(mxbx)
//                 make_edge_scal_predictor_3d(AMREX_INT_ANYD(mxbx.loVect()),
//                                             AMREX_INT_ANYD(mxbx.hiVect()), 1,
//                                             AMREX_INT_ANYD(domainBox.loVect()),
//                                             AMREX_INT_ANYD(domainBox.hiVect()),
//                                             BL_TO_FORTRAN_ANYD(scal_mf[mfi]), scal_mf.nComp(),
//                                             BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
//                                             BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
//                                             BL_TO_FORTRAN_ANYD(wmac_mf[mfi]),
//                                             BL_TO_FORTRAN_ANYD(Ip[mfi]),
//                                             BL_TO_FORTRAN_ANYD(Im[mfi]),
//                                             BL_TO_FORTRAN_ANYD(slopez[mfi]),
//                                             BL_TO_FORTRAN_ANYD(slx[mfi]),
//                                             BL_TO_FORTRAN_ANYD(srx[mfi]),
//                                             AMREX_REAL_ANYD(dx), dt, is_vel, bc_f,
//                                             nbccomp, scomp, bccomp);

                // y-direction
// #pragma gpu box(mybx)
//                 make_edge_scal_predictor_3d(AMREX_INT_ANYD(mybx.loVect()),
//                                             AMREX_INT_ANYD(mybx.hiVect()), 2,
//                                             AMREX_INT_ANYD(domainBox.loVect()),
//                                             AMREX_INT_ANYD(domainBox.hiVect()),
//                                             BL_TO_FORTRAN_ANYD(scal_mf[mfi]), scal_mf.nComp(),
//                                             BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
//                                             BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
//                                             BL_TO_FORTRAN_ANYD(wmac_mf[mfi]),
//                                             BL_TO_FORTRAN_ANYD(Ip[mfi]),
//                                             BL_TO_FORTRAN_ANYD(Im[mfi]),
//                                             BL_TO_FORTRAN_ANYD(slopez[mfi]),
//                                             BL_TO_FORTRAN_ANYD(sly[mfi]),
//                                             BL_TO_FORTRAN_ANYD(sry[mfi]),
//                                             AMREX_REAL_ANYD(dx), dt, is_vel, bc_f,
//                                             nbccomp, scomp, bccomp);

                // z-direction
#pragma gpu box(mzbx)
                make_edge_scal_predictor_3d(AMREX_INT_ANYD(mzbx.loVect()),
                                            AMREX_INT_ANYD(mzbx.hiVect()), 3,
                                            AMREX_INT_ANYD(domainBox.loVect()),
                                            AMREX_INT_ANYD(domainBox.hiVect()),
                                            BL_TO_FORTRAN_ANYD(scal_mf[mfi]), scal_mf.nComp(),
                                            BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
                                            BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
                                            BL_TO_FORTRAN_ANYD(wmac_mf[mfi]),
                                            BL_TO_FORTRAN_ANYD(Ip[mfi]),
                                            BL_TO_FORTRAN_ANYD(Im[mfi]),
                                            BL_TO_FORTRAN_ANYD(slopez[mfi]),
                                            BL_TO_FORTRAN_ANYD(slz[mfi]),
                                            BL_TO_FORTRAN_ANYD(srz[mfi]),
                                            AMREX_REAL_ANYD(dx), dt, is_vel, bc_f,
                                            nbccomp, scomp, bccomp);

                AMREX_PARALLEL_FOR_3D(mxbx, i, j, k, 
                {
                    simhx_arr(i,j,k) = (umac_arr(i,j,k) > 0.0) ? slx_arr(i,j,k) : srx_arr(i,j,k);
                    simhx_arr(i,j,k) = (abs(umac_arr(i,j,k)) > 0.0) ? simhx_arr(i,j,k) : 0.5 * (slx_arr(i,j,k) + srx_arr(i,j,k));
                });

                AMREX_PARALLEL_FOR_3D(mybx, i, j, k, 
                {
                    simhy_arr(i,j,k) = (vmac_arr(i,j,k) > 0.0) ? sly_arr(i,j,k) : sry_arr(i,j,k);
                    simhy_arr(i,j,k) = (abs(vmac_arr(i,j,k)) > 0.0) ? simhy_arr(i,j,k) : 0.5 * (sly_arr(i,j,k) + sry_arr(i,j,k));
                });

                AMREX_PARALLEL_FOR_3D(mzbx, i, j, k, 
                {
                    simhz_arr(i,j,k) = (wmac_arr(i,j,k) > 0.0) ? slz_arr(i,j,k) : srz_arr(i,j,k);
                    simhz_arr(i,j,k) = (abs(wmac_arr(i,j,k)) > 0.0) ? simhz_arr(i,j,k) : 0.5 * (slz_arr(i,j,k) + srz_arr(i,j,k));
                });

//                 // x-direction
// #pragma gpu box(mxbx)
//                 make_edge_scal_merge_3d(AMREX_INT_ANYD(mxbx.loVect()),
//                                             AMREX_INT_ANYD(mxbx.hiVect()), 
//                                             BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
//                                             BL_TO_FORTRAN_ANYD(slx[mfi]),
//                                             BL_TO_FORTRAN_ANYD(srx[mfi]),
//                                             BL_TO_FORTRAN_ANYD(simhx[mfi]));

//                 // y-direction
// #pragma gpu box(mybx)
//                 make_edge_scal_merge_3d(AMREX_INT_ANYD(mybx.loVect()),
//                                             AMREX_INT_ANYD(mybx.hiVect()), 
//                                             BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
//                                             BL_TO_FORTRAN_ANYD(sly[mfi]),
//                                             BL_TO_FORTRAN_ANYD(sry[mfi]),
//                                             BL_TO_FORTRAN_ANYD(simhy[mfi]));

//                 // z-direction
// #pragma gpu box(mzbx)
//                 make_edge_scal_merge_3d(AMREX_INT_ANYD(mzbx.loVect()),
//                                             AMREX_INT_ANYD(mzbx.hiVect()),
//                                             BL_TO_FORTRAN_ANYD(wmac_mf[mfi]),
//                                             BL_TO_FORTRAN_ANYD(slz[mfi]),
//                                             BL_TO_FORTRAN_ANYD(srz[mfi]),
//                                             BL_TO_FORTRAN_ANYD(simhz[mfi]));

                // simhxy
                Box imhbox = amrex::grow(mfi.tilebox(), 2, 1);
                imhbox = amrex::growHi(imhbox, 0, 1);
                // Box imhbox = mfi.grownnodaltilebox(0, amrex::IntVect(0,0,1)); 
#pragma gpu box(imhbox)
                make_edge_scal_transverse_3d(
                    AMREX_INT_ANYD(imhbox.loVect()), AMREX_INT_ANYD(imhbox.hiVect()),1,2,
                    AMREX_INT_ANYD(domainBox.loVect()), AMREX_INT_ANYD(domainBox.hiVect()),
                    BL_TO_FORTRAN_ANYD(scal_mf[mfi]), scal_mf.nComp(),
                    BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(wmac_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(divu[mfi]),
                    BL_TO_FORTRAN_ANYD(slx[mfi]),
                    BL_TO_FORTRAN_ANYD(srx[mfi]),
                    BL_TO_FORTRAN_ANYD(simhx[mfi]),
                    BL_TO_FORTRAN_ANYD(sly[mfi]),
                    BL_TO_FORTRAN_ANYD(sry[mfi]),
                    BL_TO_FORTRAN_ANYD(simhy[mfi]),
                    BL_TO_FORTRAN_ANYD(slz[mfi]),
                    BL_TO_FORTRAN_ANYD(srz[mfi]),
                    BL_TO_FORTRAN_ANYD(simhz[mfi]),
                    BL_TO_FORTRAN_ANYD(simhxy[mfi]),
                    AMREX_REAL_ANYD(dx), dt, is_vel, bc_f,
                    nbccomp, scomp, bccomp, is_conservative);

                // simhxz
                imhbox = amrex::grow(mfi.tilebox(), 1, 1);
                imhbox = amrex::growHi(imhbox, 0, 1);
                // imhbox = mfi.grownnodaltilebox(0, amrex::IntVect(0,1,0)); 
#pragma gpu box(imhbox)
                make_edge_scal_transverse_3d(
                    AMREX_INT_ANYD(imhbox.loVect()), AMREX_INT_ANYD(imhbox.hiVect()),1,3,
                    AMREX_INT_ANYD(domainBox.loVect()), AMREX_INT_ANYD(domainBox.hiVect()),
                    BL_TO_FORTRAN_ANYD(scal_mf[mfi]), scal_mf.nComp(),
                    BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(wmac_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(divu[mfi]),
                    BL_TO_FORTRAN_ANYD(slx[mfi]),
                    BL_TO_FORTRAN_ANYD(srx[mfi]),
                    BL_TO_FORTRAN_ANYD(simhx[mfi]),
                    BL_TO_FORTRAN_ANYD(sly[mfi]),
                    BL_TO_FORTRAN_ANYD(sry[mfi]),
                    BL_TO_FORTRAN_ANYD(simhy[mfi]),
                    BL_TO_FORTRAN_ANYD(slz[mfi]),
                    BL_TO_FORTRAN_ANYD(srz[mfi]),
                    BL_TO_FORTRAN_ANYD(simhz[mfi]),
                    BL_TO_FORTRAN_ANYD(simhxz[mfi]),
                    AMREX_REAL_ANYD(dx), dt, is_vel, bc_f,
                    nbccomp, scomp, bccomp, is_conservative);

                // simhyx
                // imhbox = mfi.grownnodaltilebox(1, amrex::IntVect(0,0,1));
                imhbox = amrex::grow(mfi.tilebox(), 2, 1);
                imhbox = amrex::growHi(imhbox, 1, 1);
#pragma gpu box(imhbox)
                make_edge_scal_transverse_3d(
                    AMREX_INT_ANYD(imhbox.loVect()), AMREX_INT_ANYD(imhbox.hiVect()),2,1,
                    AMREX_INT_ANYD(domainBox.loVect()), AMREX_INT_ANYD(domainBox.hiVect()),
                    BL_TO_FORTRAN_ANYD(scal_mf[mfi]), scal_mf.nComp(),
                    BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(wmac_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(divu[mfi]),
                    BL_TO_FORTRAN_ANYD(slx[mfi]),
                    BL_TO_FORTRAN_ANYD(srx[mfi]),
                    BL_TO_FORTRAN_ANYD(simhx[mfi]),
                    BL_TO_FORTRAN_ANYD(sly[mfi]),
                    BL_TO_FORTRAN_ANYD(sry[mfi]),
                    BL_TO_FORTRAN_ANYD(simhy[mfi]),
                    BL_TO_FORTRAN_ANYD(slz[mfi]),
                    BL_TO_FORTRAN_ANYD(srz[mfi]),
                    BL_TO_FORTRAN_ANYD(simhz[mfi]),
                    BL_TO_FORTRAN_ANYD(simhyx[mfi]),
                    AMREX_REAL_ANYD(dx), dt, is_vel, bc_f,
                    nbccomp, scomp, bccomp, is_conservative);

                // simhyz
                // imhbox = mfi.grownnodaltilebox(1, amrex::IntVect(1,0,0)); 
                imhbox = amrex::grow(mfi.tilebox(), 0, 1);
                imhbox = amrex::growHi(imhbox, 1, 1);
#pragma gpu box(imhbox)
                make_edge_scal_transverse_3d(
                    AMREX_INT_ANYD(imhbox.loVect()), AMREX_INT_ANYD(imhbox.hiVect()),2,3,
                    AMREX_INT_ANYD(domainBox.loVect()), AMREX_INT_ANYD(domainBox.hiVect()),
                    BL_TO_FORTRAN_ANYD(scal_mf[mfi]), scal_mf.nComp(),
                    BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(wmac_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(divu[mfi]),
                    BL_TO_FORTRAN_ANYD(slx[mfi]),
                    BL_TO_FORTRAN_ANYD(srx[mfi]),
                    BL_TO_FORTRAN_ANYD(simhx[mfi]),
                    BL_TO_FORTRAN_ANYD(sly[mfi]),
                    BL_TO_FORTRAN_ANYD(sry[mfi]),
                    BL_TO_FORTRAN_ANYD(simhy[mfi]),
                    BL_TO_FORTRAN_ANYD(slz[mfi]),
                    BL_TO_FORTRAN_ANYD(srz[mfi]),
                    BL_TO_FORTRAN_ANYD(simhz[mfi]),
                    BL_TO_FORTRAN_ANYD(simhyz[mfi]),
                    AMREX_REAL_ANYD(dx), dt, is_vel, bc_f,
                    nbccomp, scomp, bccomp, is_conservative);

                // simhzx
                // imhbox = mfi.grownnodaltilebox(2, amrex::IntVect(0,1,0));
                imhbox = amrex::grow(mfi.tilebox(), 1, 1);
                imhbox = amrex::growHi(imhbox, 2, 1);
#pragma gpu box(imhbox)
                make_edge_scal_transverse_3d(
                    AMREX_INT_ANYD(imhbox.loVect()), AMREX_INT_ANYD(imhbox.hiVect()),3,1,
                    AMREX_INT_ANYD(domainBox.loVect()), AMREX_INT_ANYD(domainBox.hiVect()),
                    BL_TO_FORTRAN_ANYD(scal_mf[mfi]), scal_mf.nComp(),
                    BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(wmac_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(divu[mfi]),
                    BL_TO_FORTRAN_ANYD(slx[mfi]),
                    BL_TO_FORTRAN_ANYD(srx[mfi]),
                    BL_TO_FORTRAN_ANYD(simhx[mfi]),
                    BL_TO_FORTRAN_ANYD(sly[mfi]),
                    BL_TO_FORTRAN_ANYD(sry[mfi]),
                    BL_TO_FORTRAN_ANYD(simhy[mfi]),
                    BL_TO_FORTRAN_ANYD(slz[mfi]),
                    BL_TO_FORTRAN_ANYD(srz[mfi]),
                    BL_TO_FORTRAN_ANYD(simhz[mfi]),
                    BL_TO_FORTRAN_ANYD(simhzx[mfi]),
                    AMREX_REAL_ANYD(dx), dt, is_vel, bc_f,
                    nbccomp, scomp, bccomp, is_conservative);

                // simhzy
                // imhbox = mfi.grownnodaltilebox(2, IntVect(1,0,0));
                imhbox = amrex::grow(mfi.tilebox(), 0, 1);
                imhbox = amrex::growHi(imhbox, 2, 1);
#pragma gpu box(imhbox)
                make_edge_scal_transverse_3d(
                    AMREX_INT_ANYD(imhbox.loVect()), AMREX_INT_ANYD(imhbox.hiVect()),3,2,
                    AMREX_INT_ANYD(domainBox.loVect()), AMREX_INT_ANYD(domainBox.hiVect()),
                    BL_TO_FORTRAN_ANYD(scal_mf[mfi]), scal_mf.nComp(),
                    BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(wmac_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(divu[mfi]),
                    BL_TO_FORTRAN_ANYD(slx[mfi]),
                    BL_TO_FORTRAN_ANYD(srx[mfi]),
                    BL_TO_FORTRAN_ANYD(simhx[mfi]),
                    BL_TO_FORTRAN_ANYD(sly[mfi]),
                    BL_TO_FORTRAN_ANYD(sry[mfi]),
                    BL_TO_FORTRAN_ANYD(simhy[mfi]),
                    BL_TO_FORTRAN_ANYD(slz[mfi]),
                    BL_TO_FORTRAN_ANYD(srz[mfi]),
                    BL_TO_FORTRAN_ANYD(simhz[mfi]),
                    BL_TO_FORTRAN_ANYD(simhzy[mfi]),
                    AMREX_REAL_ANYD(dx), dt, is_vel, bc_f,
                    nbccomp, scomp, bccomp, is_conservative);

                // x-direction
#pragma gpu box(xbx)
                make_edge_scal_3d(AMREX_INT_ANYD(xbx.loVect()),
                                  AMREX_INT_ANYD(xbx.hiVect()),1,
                                  AMREX_INT_ANYD(domainBox.loVect()),
                                  AMREX_INT_ANYD(domainBox.hiVect()),
                                  BL_TO_FORTRAN_ANYD(scal_mf[mfi]), scal_mf.nComp(),
                                  BL_TO_FORTRAN_ANYD(sedgex_mf[mfi]), sedgex_mf.nComp(),
                                  BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
                                  BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
                                  BL_TO_FORTRAN_ANYD(wmac_mf[mfi]),
                                  BL_TO_FORTRAN_ANYD(Ipf[mfi]),
                                  BL_TO_FORTRAN_ANYD(Imf[mfi]),
                                  BL_TO_FORTRAN_ANYD(slx[mfi]),
                                  BL_TO_FORTRAN_ANYD(srx[mfi]),
                                  BL_TO_FORTRAN_ANYD(simhxy[mfi]),
                                  BL_TO_FORTRAN_ANYD(simhxz[mfi]),
                                  BL_TO_FORTRAN_ANYD(simhyx[mfi]),
                                  BL_TO_FORTRAN_ANYD(simhyz[mfi]),
                                  BL_TO_FORTRAN_ANYD(simhzx[mfi]),
                                  BL_TO_FORTRAN_ANYD(simhzy[mfi]),
                                  BL_TO_FORTRAN_ANYD(force_mf[mfi]), force_mf.nComp(),
                                  AMREX_REAL_ANYD(dx), dt, is_vel, bc_f,
                                  nbccomp, scomp, bccomp, is_conservative);

                // y-direction
#pragma gpu box(ybx)
                make_edge_scal_3d(AMREX_INT_ANYD(ybx.loVect()),
                                  AMREX_INT_ANYD(ybx.hiVect()),2,
                                  AMREX_INT_ANYD(domainBox.loVect()),
                                  AMREX_INT_ANYD(domainBox.hiVect()),
                                  BL_TO_FORTRAN_ANYD(scal_mf[mfi]), scal_mf.nComp(),
                                  BL_TO_FORTRAN_ANYD(sedgey_mf[mfi]), sedgey_mf.nComp(),
                                  BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
                                  BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
                                  BL_TO_FORTRAN_ANYD(wmac_mf[mfi]),
                                  BL_TO_FORTRAN_ANYD(Ipf[mfi]),
                                  BL_TO_FORTRAN_ANYD(Imf[mfi]),
                                  BL_TO_FORTRAN_ANYD(sly[mfi]),
                                  BL_TO_FORTRAN_ANYD(sry[mfi]),
                                  BL_TO_FORTRAN_ANYD(simhxy[mfi]),
                                  BL_TO_FORTRAN_ANYD(simhxz[mfi]),
                                  BL_TO_FORTRAN_ANYD(simhyx[mfi]),
                                  BL_TO_FORTRAN_ANYD(simhyz[mfi]),
                                  BL_TO_FORTRAN_ANYD(simhzx[mfi]),
                                  BL_TO_FORTRAN_ANYD(simhzy[mfi]),
                                  BL_TO_FORTRAN_ANYD(force_mf[mfi]), force_mf.nComp(),
                                  AMREX_REAL_ANYD(dx), dt, is_vel, bc_f,
                                  nbccomp, scomp, bccomp, is_conservative);

                // z-direction
#pragma gpu box(zbx)
                make_edge_scal_3d(AMREX_INT_ANYD(zbx.loVect()),
                                  AMREX_INT_ANYD(zbx.hiVect()),3,
                                  AMREX_INT_ANYD(domainBox.loVect()),
                                  AMREX_INT_ANYD(domainBox.hiVect()),
                                  BL_TO_FORTRAN_ANYD(scal_mf[mfi]), scal_mf.nComp(),
                                  BL_TO_FORTRAN_ANYD(sedgez_mf[mfi]), sedgez_mf.nComp(),
                                  BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
                                  BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
                                  BL_TO_FORTRAN_ANYD(wmac_mf[mfi]),
                                  BL_TO_FORTRAN_ANYD(Ipf[mfi]),
                                  BL_TO_FORTRAN_ANYD(Imf[mfi]),
                                  BL_TO_FORTRAN_ANYD(slz[mfi]),
                                  BL_TO_FORTRAN_ANYD(srz[mfi]),
                                  BL_TO_FORTRAN_ANYD(simhxy[mfi]),
                                  BL_TO_FORTRAN_ANYD(simhxz[mfi]),
                                  BL_TO_FORTRAN_ANYD(simhyx[mfi]),
                                  BL_TO_FORTRAN_ANYD(simhyz[mfi]),
                                  BL_TO_FORTRAN_ANYD(simhzx[mfi]),
                                  BL_TO_FORTRAN_ANYD(simhzy[mfi]),
                                  BL_TO_FORTRAN_ANYD(force_mf[mfi]), force_mf.nComp(),
                                  AMREX_REAL_ANYD(dx), dt, is_vel, bc_f,
                                  nbccomp, scomp, bccomp, is_conservative);

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

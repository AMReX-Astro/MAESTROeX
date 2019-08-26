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

#ifdef AMREX_USE_CUDA
    auto not_launched = Gpu::notInLaunchRegion();
    // turn on GPU
    if (not_launched) Gpu::setLaunchRegion(true);
#endif

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

#ifdef AMREX_USE_CUDA
    // turn off GPU
    if (not_launched) Gpu::setLaunchRegion(false);
#endif

}

void
Maestro::Put1dArrayOnCart (int level,
                           const RealVector& s0,
                           Vector<MultiFab>& s0_cart,
                           int is_input_edge_centered,
                           int is_output_a_vector,
                           const Vector<BCRec>& bcs,
                           int sbccomp)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Put1dArrayOnCart_lev()",Put1dArrayOnCart);

#ifdef AMREX_USE_CUDA
    auto not_launched = Gpu::notInLaunchRegion();
    // turn on GPU
    if (not_launched) Gpu::setLaunchRegion(true);
#endif

    // get references to the MultiFabs at level lev
    MultiFab& s0_cart_mf = s0_cart[level];
    MultiFab& cc_to_r = cell_cc_to_r[level];

    // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
    for ( MFIter mfi(s0_cart_mf, true); mfi.isValid(); ++mfi ) {

    	// Get the index space of the valid region
    	const Box& tileBox = mfi.tilebox();
    	const Real* dx = geom[level].CellSize();

    	// call fortran subroutine
    	// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
    	// lo/hi coordinates (including ghost cells), and/or the # of components
    	// We will also pass "validBox", which specifies the "valid" region.
    	if (spherical == 0) {
#pragma gpu box(tileBox)
    	    put_1d_array_on_cart(AMREX_INT_ANYD(tileBox.loVect()),
                     AMREX_INT_ANYD(tileBox.hiVect()),level,
    				 BL_TO_FORTRAN_ANYD(s0_cart_mf[mfi]), s0_cart_mf.nComp(),
    				 s0.dataPtr(), is_input_edge_centered, is_output_a_vector);
    	} else {
#pragma gpu box(tileBox)
    	    put_1d_array_on_cart_sphr(AMREX_INT_ANYD(tileBox.loVect()),
                          AMREX_INT_ANYD(tileBox.hiVect()),
    				      BL_TO_FORTRAN_ANYD(s0_cart_mf[mfi]), s0_cart_mf.nComp(),
    				      s0.dataPtr(), AMREX_REAL_ANYD(dx),
    				      is_input_edge_centered, is_output_a_vector,
    				      r_cc_loc.dataPtr(), r_edge_loc.dataPtr(),
    				      BL_TO_FORTRAN_ANYD(cc_to_r[mfi]));
    	}
    }

#ifdef AMREX_USE_CUDA
    // turn off GPU
    if (not_launched) Gpu::setLaunchRegion(false);
#endif

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
        // NOTE: don't use tiling - put OpenMP into F90
        for ( MFIter mfi(sold_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            if (spherical == 0) {

                addw0(&lev,ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                      BL_TO_FORTRAN_3D(uedge_mf[mfi]),
                      BL_TO_FORTRAN_3D(vedge_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
                      BL_TO_FORTRAN_3D(wedge_mf[mfi]),
#endif
                      w0.dataPtr(),&mult);

            } else {

#if (AMREX_SPACEDIM == 3)
                addw0_sphr(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                           BL_TO_FORTRAN_3D(uedge_mf[mfi]),
                           BL_TO_FORTRAN_3D(vedge_mf[mfi]),
                           BL_TO_FORTRAN_3D(wedge_mf[mfi]),
                           BL_TO_FORTRAN_3D(w0macx_mf[mfi]),
                           BL_TO_FORTRAN_3D(w0macy_mf[mfi]),
                           BL_TO_FORTRAN_3D(w0macz_mf[mfi]),
                           &mult);
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

    // Construct a cartesian version of w0
    Vector<MultiFab> w0_cart(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        w0_cart[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 2);
        w0_cart[lev].setVal(0.);
    }

    if (w0mac_interp_type == 1) {
        Put1dArrayOnCart(w0, w0_cart, 1, 1, bcs_u, 0, 1);
    }

    if (w0mac[0][0].nGrow() != 1) {
        Abort("Error: MakeW0mac assumes one ghost cell");
    }

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& w0macx_mf = w0mac[lev][0];
        MultiFab& w0macy_mf = w0mac[lev][1];
        MultiFab& w0macz_mf = w0mac[lev][2];
        MultiFab& w0cart_mf = w0_cart[lev];
        const Real* dx = geom[lev].CellSize();

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(w0cart_mf, true); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            make_w0mac_sphr(ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
                            w0.dataPtr(),
                            BL_TO_FORTRAN_3D(w0macx_mf[mfi]),
                            BL_TO_FORTRAN_3D(w0macy_mf[mfi]),
                            BL_TO_FORTRAN_3D(w0macz_mf[mfi]),
                            BL_TO_FORTRAN_FAB(w0cart_mf[mfi]),
                            dx, r_edge_loc.dataPtr());
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

    // Construct a cartesian version of w0
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
        for ( MFIter mfi(s0cart_mf, true); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
	    if (use_exact_base_state) {
		make_s0mac_sphr_irreg(ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
				      s0.dataPtr(),
				      BL_TO_FORTRAN_3D(s0macx_mf[mfi]),
				      BL_TO_FORTRAN_3D(s0macy_mf[mfi]),
				      BL_TO_FORTRAN_3D(s0macz_mf[mfi]),
				      BL_TO_FORTRAN_3D(s0cart_mf[mfi]),
				      dx, r_cc_loc.dataPtr());
	    } else {
		make_s0mac_sphr(ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
				s0.dataPtr(),
				BL_TO_FORTRAN_3D(s0macx_mf[mfi]),
				BL_TO_FORTRAN_3D(s0macy_mf[mfi]),
				BL_TO_FORTRAN_3D(s0macz_mf[mfi]),
				BL_TO_FORTRAN_3D(s0cart_mf[mfi]),
				dx, r_cc_loc.dataPtr());
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
        for ( MFIter mfi(normal_mf, true); mfi.isValid(); ++mfi ) {

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            make_normal(BL_TO_FORTRAN_3D(normal_mf[mfi]), dx);
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
        for ( MFIter mfi(scc_mf, true); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            put_data_on_faces(ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
                              BL_TO_FORTRAN_3D(scc_mf[mfi]),
                              BL_TO_FORTRAN_3D(facex_mf[mfi]),
                              BL_TO_FORTRAN_3D(facey_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
                              BL_TO_FORTRAN_3D(facez_mf[mfi]),
#endif
                              &harmonic_avg);
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
        for ( MFIter mfi(cc_to_r); mfi.isValid(); ++mfi ) {

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            init_base_state_map_sphr(BL_TO_FORTRAN_3D(cc_to_r[mfi]),
				     ZFILL(dx_fine),
				     ZFILL(dx));
        }
    }

}

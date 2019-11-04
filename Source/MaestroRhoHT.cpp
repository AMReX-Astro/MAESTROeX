
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void
Maestro::TfromRhoH (Vector<MultiFab>& scal,
                    const RealVector& p0)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::TfromRhoH()",TfromRhoH);

#ifdef AMREX_USE_CUDA
    auto not_launched = Gpu::notInLaunchRegion();
    // turn on GPU
    if (not_launched) Gpu::setLaunchRegion(true);
#endif

    Vector<MultiFab> p0_cart(finest_level+1);

    for (int lev=0; lev<=finest_level; ++lev) {
        p0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
        p0_cart[lev].setVal(0.);
    }
    Put1dArrayOnCart(p0,p0_cart,0,0,bcs_f,0);

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& scal_mf = scal[lev];
        const MultiFab& p0_mf = p0_cart[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(scal_mf, true); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.

#pragma gpu box(tileBox)
            makeTfromRhoH(AMREX_INT_ANYD(tileBox.loVect()),
                                AMREX_INT_ANYD(tileBox.hiVect()),
                                BL_TO_FORTRAN_ANYD(scal_mf[mfi]),
                                BL_TO_FORTRAN_ANYD(p0_mf[mfi]));
        }

    }

    // average down and fill ghost cells
    AverageDown(scal,Temp,1);
    FillPatch(t_old,scal,scal,scal,Temp,Temp,1,Temp,bcs_s);

#ifdef AMREX_USE_CUDA
    // turn off GPU
    if (not_launched) Gpu::setLaunchRegion(false);
#endif
}

void
Maestro::TfromRhoP (Vector<MultiFab>& scal,
                    const RealVector& p0,
                    int updateRhoH)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::TfromRhoP()",TfromRhoP);

#ifdef AMREX_USE_CUDA
    auto not_launched = Gpu::notInLaunchRegion();
    // turn on GPU
    if (not_launched) Gpu::setLaunchRegion(true);
#endif

    Vector<MultiFab> p0_cart(finest_level+1);

    for (int lev=0; lev<=finest_level; ++lev) {
        p0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
        p0_cart[lev].setVal(0.);
    }
    Put1dArrayOnCart(p0,p0_cart,0,0,bcs_f,0);

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& scal_mf = scal[lev];
        const MultiFab& p0_mf = p0_cart[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(scal_mf, true); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.

#pragma gpu box(tileBox)
            makeTfromRhoP(AMREX_INT_ANYD(tileBox.loVect()),
                                AMREX_INT_ANYD(tileBox.hiVect()),
                                BL_TO_FORTRAN_ANYD(scal_mf[mfi]),
                                BL_TO_FORTRAN_ANYD(p0_mf[mfi]), updateRhoH);
        }
    }

    // average down and fill ghost cells (Temperature)
    AverageDown(scal,Temp,1);
    FillPatch(t_old,scal,scal,scal,Temp,Temp,1,Temp,bcs_s);

    // average down and fill ghost cells (Enthalpy)
    if (updateRhoH == 1) {
        AverageDown(scal,RhoH,1);
        FillPatch(t_old,scal,scal,scal,RhoH,RhoH,1,RhoH,bcs_s);
    }

#ifdef AMREX_USE_CUDA
    // turn off GPU
    if (not_launched) Gpu::setLaunchRegion(false);
#endif
}

void
Maestro::PfromRhoH (const Vector<MultiFab>& state,
                    const Vector<MultiFab>& s_old,
                    Vector<MultiFab>& peos)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::PfromRhoH()",PfromRhoH);

#ifdef AMREX_USE_CUDA
    auto not_launched = Gpu::notInLaunchRegion();
    // turn on GPU
    if (not_launched) Gpu::setLaunchRegion(true);
#endif

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& state_mf = state[lev];
        const MultiFab& sold_mf = s_old[lev];
        MultiFab& peos_mf = peos[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(state_mf, true); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
#pragma gpu box(tileBox)
            makePfromRhoH(AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),
                          BL_TO_FORTRAN_ANYD(state_mf[mfi]),
                          sold_mf[mfi].dataPtr(Temp),
                          AMREX_INT_ANYD(sold_mf[mfi].loVect()), AMREX_INT_ANYD(sold_mf[mfi].hiVect()),
                          BL_TO_FORTRAN_ANYD(peos_mf[mfi]));
        }

    }

    // average down and fill ghost cells
    AverageDown(peos,0,1);
    FillPatch(t_old,peos,peos,peos,0,0,1,0,bcs_f);

#ifdef AMREX_USE_CUDA
    // turn off GPU
    if (not_launched) Gpu::setLaunchRegion(false);
#endif
}

void
Maestro::MachfromRhoH (const Vector<MultiFab>& scal,
                           const Vector<MultiFab>& vel,
                           const RealVector& p0,
                           const Vector<MultiFab>& w0cart,
                           Vector<MultiFab>& mach)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MachfromRhoH()",MachfromRhoH);

    Vector<MultiFab> p0_cart(finest_level+1);

    for (int lev=0; lev<=finest_level; ++lev) {
        p0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
        p0_cart[lev].setVal(0.);
    }
    Put1dArrayOnCart(p0,p0_cart,0,0,bcs_f,0);

#ifdef AMREX_USE_CUDA
    auto not_launched = Gpu::notInLaunchRegion();
    // turn on GPU
    if (not_launched) Gpu::setLaunchRegion(true);
#endif

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& scal_mf = scal[lev];
        const MultiFab& vel_mf = vel[lev];
        const MultiFab& w0cart_mf = w0cart[lev];
        MultiFab& mach_mf = mach[lev];
        const MultiFab& p0_mf = p0_cart[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(scal_mf, true); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const Real* dx = geom[lev].CellSize();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
#pragma gpu box(tileBox)
            makeMachfromRhoH(AMREX_INT_ANYD(tileBox.loVect()),
                                  AMREX_INT_ANYD(tileBox.hiVect()),
                                  BL_TO_FORTRAN_ANYD(scal_mf[mfi]),
                                  BL_TO_FORTRAN_ANYD(vel_mf[mfi]),
                                  BL_TO_FORTRAN_ANYD(p0_mf[mfi]),
                                  BL_TO_FORTRAN_ANYD(w0cart_mf[mfi]),
                                  BL_TO_FORTRAN_ANYD(mach_mf[mfi]));
        }

    }

    // average down and fill ghost cells
    AverageDown(mach,0,1);
    FillPatch(t_old,mach,mach,mach,0,0,1,0,bcs_f);

#ifdef AMREX_USE_CUDA
    // turn off GPU
    if (not_launched) Gpu::setLaunchRegion(false);
#endif
}

void
Maestro::CsfromRhoH (const Vector<MultiFab>& scal,
                     const RealVector& p0,
                     const Vector<MultiFab>& p0cart,
                     Vector<MultiFab>& cs)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::CsfromRhoH()",CsfromRhoH);

#ifdef AMREX_USE_CUDA
    auto not_launched = Gpu::notInLaunchRegion();
    // turn on GPU
    if (not_launched) Gpu::setLaunchRegion(true);
#endif

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& scal_mf = scal[lev];
        MultiFab& cs_mf = cs[lev];

        const MultiFab& p0cart_mf = p0cart[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(scal_mf, true); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
#pragma gpu box(tileBox)
            makeCsfromRhoH(AMREX_INT_ANYD(tileBox.loVect()),
                                AMREX_INT_ANYD(tileBox.hiVect()),
                                BL_TO_FORTRAN_ANYD(scal_mf[mfi]),
                                BL_TO_FORTRAN_ANYD(p0cart_mf[mfi]),
                                BL_TO_FORTRAN_ANYD(cs_mf[mfi]));
        }
    }

    // average down and fill ghost cells
    AverageDown(cs,0,1);
    FillPatch(t_old,cs,cs,cs,0,0,1,0,bcs_f);

#ifdef AMREX_USE_CUDA
    // turn off GPU
    if (not_launched) Gpu::setLaunchRegion(false);
#endif
}

void
Maestro::HfromRhoTedge (Vector<std::array< MultiFab, AMREX_SPACEDIM > >& sedge,
			const RealVector& rho0_edge_old,
			const RealVector& rhoh0_edge_old,
			const RealVector& rho0_edge_new,
			const RealVector& rhoh0_edge_new)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::HfromRhoTedge()",HfromRhoTedge);

// #ifdef AMREX_USE_CUDA
//     auto not_launched = Gpu::notInLaunchRegion();
//     // turn on GPU
//     if (not_launched) Gpu::setLaunchRegion(true);
// #endif

    Vector<MultiFab> rho0_cart(finest_level+1);
    Vector<MultiFab> rhoh0_cart(finest_level+1);
    Vector<MultiFab> tempbar_cart(finest_level+1);

    for (int lev=0; lev<=finest_level; ++lev) {
        rho0_cart[lev].define(grids[lev], dmap[lev], 1, 2);
        rhoh0_cart[lev].define(grids[lev], dmap[lev], 1, 2);
        tempbar_cart[lev].define(grids[lev], dmap[lev], 1, 2);
        rho0_cart[lev].setVal(0.);
	rhoh0_cart[lev].setVal(0.);
	tempbar_cart[lev].setVal(0.);
    }

    RealVector rho0_halftime( (max_radial_level+1)*nr_fine );
    RealVector rhoh0_halftime( (max_radial_level+1)*nr_fine );
    rho0_halftime.shrink_to_fit();
    rhoh0_halftime.shrink_to_fit();

    for (int i = 0; i < (max_radial_level+1)*nr_fine; ++i) {
	rho0_halftime[i] = 0.5*(rho0_old[i] + rho0_new[i]);
	rhoh0_halftime[i] = 0.5*(rhoh0_old[i] + rhoh0_new[i]);
    }
    
    Put1dArrayOnCart(rho0_halftime,rho0_cart,0,0,bcs_s,Rho);
    Put1dArrayOnCart(rhoh0_halftime,rhoh0_cart,0,0,bcs_s,RhoH);
    Put1dArrayOnCart(tempbar,tempbar_cart,0,0,bcs_s,Temp);

    // edge variables    
    RealVector rho0_edge( (max_radial_level+1)*(nr_fine+1) );
    RealVector rhoh0_edge( (max_radial_level+1)*(nr_fine+1) );
    RealVector tempbar_edge( (max_radial_level+1)*(nr_fine+1) );
    rho0_edge.shrink_to_fit();
    rhoh0_edge.shrink_to_fit();
    tempbar_edge.shrink_to_fit();

    if (spherical == 0) {
	cell_to_edge(tempbar.dataPtr(), tempbar_edge.dataPtr());
	for (int i = 0; i < (max_radial_level+1)*(nr_fine+1); ++i) {
	    rho0_edge[i] = 0.5*(rho0_edge_old[i] + rho0_edge_new[i]);
	    rhoh0_edge[i] = 0.5*(rhoh0_edge_old[i] + rhoh0_edge_new[i]);
	}
    }
    
    Vector<MultiFab> rho0_edge_cart(finest_level+1);
    Vector<MultiFab> rhoh0_edge_cart(finest_level+1);
    Vector<MultiFab> tempbar_edge_cart(finest_level+1);

    for (int lev=0; lev<=finest_level; ++lev) {
        rho0_edge_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        rhoh0_edge_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        tempbar_edge_cart[lev].define(grids[lev], dmap[lev], 1, 1);
	rho0_edge_cart[lev].setVal(0.);
	rhoh0_edge_cart[lev].setVal(0.);
	tempbar_edge_cart[lev].setVal(0.);
    }

    if (spherical == 0) {
        Put1dArrayOnCart(rho0_edge,rho0_edge_cart,1,0,bcs_s,Rho);
        Put1dArrayOnCart(rhoh0_edge,rhoh0_edge_cart,1,0,bcs_s,RhoH);
        Put1dArrayOnCart(tempbar_edge,tempbar_edge_cart,1,0,bcs_s,Temp);
    }


    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
	const MultiFab& scal_mf = sold[lev];
        MultiFab& sedgex_mf = sedge[lev][0];
        MultiFab& sedgey_mf = sedge[lev][1];
#if (AMREX_SPACEDIM == 3)
        MultiFab& sedgez_mf = sedge[lev][2];
#endif
        const MultiFab& rho0_mf = rho0_cart[lev];
        const MultiFab& rhoh0_mf = rhoh0_cart[lev];
        const MultiFab& tempbar_mf = tempbar_cart[lev];
        const MultiFab& rho0_edge_mf = rho0_edge_cart[lev];
        const MultiFab& rhoh0_edge_mf = rhoh0_edge_cart[lev];
        const MultiFab& tempbar_edge_mf = tempbar_edge_cart[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(scal_mf, true); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.

	    if (spherical == 0) {
		makeHfromRhoT_edge((tileBox.loVect()),
				   (tileBox.hiVect()),
				   BL_TO_FORTRAN_ANYD(sedgex_mf[mfi]),
				   BL_TO_FORTRAN_ANYD(sedgey_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
				   BL_TO_FORTRAN_ANYD(sedgez_mf[mfi]),
#endif
				   BL_TO_FORTRAN_ANYD(rho0_mf[mfi]),
				   BL_TO_FORTRAN_ANYD(rhoh0_mf[mfi]),
				   BL_TO_FORTRAN_ANYD(tempbar_mf[mfi]),
				   BL_TO_FORTRAN_ANYD(rho0_edge_mf[mfi]),
				   BL_TO_FORTRAN_ANYD(rhoh0_edge_mf[mfi]),
				   BL_TO_FORTRAN_ANYD(tempbar_edge_mf[mfi]));
	    } else {
#if (AMREX_SPACEDIM == 3)
		makeHfromRhoT_edge_sphr(tileBox.loVect(),
					tileBox.hiVect(),
					BL_TO_FORTRAN_ANYD(sedgex_mf[mfi]),
					BL_TO_FORTRAN_ANYD(sedgey_mf[mfi]),
					BL_TO_FORTRAN_ANYD(sedgez_mf[mfi]),
					BL_TO_FORTRAN_ANYD(rho0_mf[mfi]),
					BL_TO_FORTRAN_ANYD(rhoh0_mf[mfi]),
					BL_TO_FORTRAN_ANYD(tempbar_mf[mfi]));
#endif
		
	    }
	}

    }

// #ifdef AMREX_USE_CUDA
//     // turn off GPU
//     if (not_launched) Gpu::setLaunchRegion(false);
// #endif
}

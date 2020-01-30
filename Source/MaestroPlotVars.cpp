#include <Maestro.H>
#include <Maestro_F.H>
#include <MaestroPlot.H>
#include <AMReX_buildInfo.H>
#include <iterator>     // std::istream_iterator

using namespace amrex;

void
Maestro::MakeMagvel (const Vector<MultiFab>& vel,
                     Vector<MultiFab>& magvel)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeMagvel()",MakeMagvel);

    Vector<std::array< MultiFab, AMREX_SPACEDIM > > w0mac(finest_level+1);

#if (AMREX_SPACEDIM == 3)
    if (spherical == 1) {
        for (int lev=0; lev<=finest_level; ++lev) {
            w0mac[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);
            w0mac[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
            w0mac[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1);
        }
        MakeW0mac(w0mac);
    }
#endif

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& vel_mf = vel[lev];
        MultiFab& magvel_mf = magvel[lev];
        const MultiFab& w0_mf = w0_cart[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        if (spherical == 0) {
#ifdef _OPENMP
#pragma omp parallel
#endif
            for ( MFIter mfi(vel_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

                // Get the index space of the valid region
                const Box& tileBox = mfi.tilebox();

                // call fortran subroutine
                // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
                // lo/hi coordinates (including ghost cells), and/or the # of components
                // We will also pass "validBox", which specifies the "valid" region.
#pragma gpu box(tileBox)
                make_magvel(AMREX_INT_ANYD(tileBox.loVect()),
                            AMREX_INT_ANYD(tileBox.hiVect()),
                            BL_TO_FORTRAN_ANYD(vel_mf[mfi]),
                            BL_TO_FORTRAN_ANYD(w0_mf[mfi]), w0_mf.nComp(),
                            BL_TO_FORTRAN_ANYD(magvel_mf[mfi]));
            }

        } else {

#ifdef _OPENMP
#pragma omp parallel
#endif
            for ( MFIter mfi(vel_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

                // Get the index space of the valid region
                const Box& tileBox = mfi.tilebox();

                MultiFab& w0macx_mf = w0mac[lev][0];
                MultiFab& w0macy_mf = w0mac[lev][1];
                MultiFab& w0macz_mf = w0mac[lev][2];

                // call fortran subroutine
                // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
                // lo/hi coordinates (including ghost cells), and/or the # of components
                // We will also pass "validBox", which specifies the "valid" region.
#pragma gpu box(tileBox)
                make_magvel_sphr(AMREX_INT_ANYD(tileBox.loVect()),
                                 AMREX_INT_ANYD(tileBox.hiVect()),
                                 BL_TO_FORTRAN_ANYD(vel_mf[mfi]),
                                 BL_TO_FORTRAN_ANYD(w0macx_mf[mfi]),
                                 BL_TO_FORTRAN_ANYD(w0macy_mf[mfi]),
                                 BL_TO_FORTRAN_ANYD(w0macz_mf[mfi]),
                                 BL_TO_FORTRAN_ANYD(magvel_mf[mfi]));
            }
        }
    }

    // average down and fill ghost cells
    AverageDown(magvel,0,1);
    FillPatch(t_old,magvel,magvel,magvel,0,0,1,0,bcs_f);
}


void
Maestro::MakeVelrc (const Vector<MultiFab>& vel,
                    const Vector<MultiFab>& w0rcart,
                    Vector<MultiFab>& rad_vel,
                    Vector<MultiFab>& circ_vel)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeVelrc()",MakeVelrc);

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& vel_mf = vel[lev];
        MultiFab& radvel_mf = rad_vel[lev];
        MultiFab& circvel_mf = circ_vel[lev];
        const MultiFab& w0rcart_mf = w0rcart[lev];
        const MultiFab& normal_mf = normal[lev];

#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(vel_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
#pragma gpu box(tileBox)
            make_velrc(AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),
                       BL_TO_FORTRAN_ANYD(vel_mf[mfi]),
                       BL_TO_FORTRAN_ANYD(w0rcart_mf[mfi]),
                       BL_TO_FORTRAN_ANYD(normal_mf[mfi]),
                       BL_TO_FORTRAN_ANYD(radvel_mf[mfi]),
                       BL_TO_FORTRAN_ANYD(circvel_mf[mfi]));
        }
    }

    // average down and fill ghost cells
    AverageDown(rad_vel,0,1);
    FillPatch(t_old,rad_vel,rad_vel,rad_vel,0,0,1,0,bcs_f);
    AverageDown(circ_vel,0,1);
    FillPatch(t_old,circ_vel,circ_vel,circ_vel,0,0,1,0,bcs_f);
}


void
Maestro::MakeAdExcess (const Vector<MultiFab>& state,
                       Vector<MultiFab>& ad_excess)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeAdExcess()",MakeAdExcess);

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& state_mf = state[lev];
        MultiFab& ad_excess_mf = ad_excess[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        if (spherical == 0) {
#ifdef _OPENMP
#pragma omp parallel
#endif
            for ( MFIter mfi(state_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

                // Get the index space of the valid region
                const Box& tileBox = mfi.tilebox();

                // call fortran subroutine
                // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
                // lo/hi coordinates (including ghost cells), and/or the # of components
                // We will also pass "validBox", which specifies the "valid" region.
#pragma gpu box(tileBox)
                make_ad_excess(AMREX_INT_ANYD(tileBox.loVect()),
                               AMREX_INT_ANYD(tileBox.hiVect()),
                               BL_TO_FORTRAN_ANYD(state_mf[mfi]), state_mf.nComp(),
                               BL_TO_FORTRAN_ANYD(ad_excess_mf[mfi]));
            }

        } else {

            const MultiFab& normal_mf = normal[lev];

#ifdef _OPENMP
#pragma omp parallel
#endif
            for ( MFIter mfi(state_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

                // Get the index space of the valid region
                const Box& tileBox = mfi.tilebox();

                // call fortran subroutine
                // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
                // lo/hi coordinates (including ghost cells), and/or the # of components
                // We will also pass "validBox", which specifies the "valid" region.
#pragma gpu box(tileBox)
                make_ad_excess_sphr(AMREX_INT_ANYD(tileBox.loVect()),
                                    AMREX_INT_ANYD(tileBox.hiVect()),
                                    BL_TO_FORTRAN_ANYD(state_mf[mfi]), state_mf.nComp(),
                                    BL_TO_FORTRAN_ANYD(normal_mf[mfi]),
                                    BL_TO_FORTRAN_ANYD(ad_excess_mf[mfi]));
            }
        }
    }

    // average down and fill ghost cells
    AverageDown(ad_excess,0,1);
    FillPatch(t_old,ad_excess,ad_excess,ad_excess,0,0,1,0,bcs_f);
}


void
Maestro::MakeGrav (const RealVector& rho0,
                   Vector<MultiFab>& grav)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeGrav()",MakeGrav);

    RealVector grav_cell( (max_radial_level+1)*nr_fine );
    grav_cell.shrink_to_fit();

    make_grav_cell(grav_cell.dataPtr(),
                   rho0.dataPtr(),
                   r_cc_loc.dataPtr(),
                   r_edge_loc.dataPtr());

    Put1dArrayOnCart(grav_cell,grav,0,0,bcs_f,0);

    // average down and fill ghost cells
    AverageDown(grav,0,1);
    FillPatch(t_old,grav,grav,grav,0,0,1,0,bcs_f);
}


void
Maestro::MakeVorticity (const Vector<MultiFab>& vel,
                        Vector<MultiFab>& vorticity)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeVorticity()",MakeVorticity);

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& vel_mf = vel[lev];
        MultiFab& vorticity_mf = vorticity[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(vel_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const Real* dx = geom[lev].CellSize();
            const Box& domainBox = geom[lev].Domain();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.

            // NOTE: do not offload to the gpu as make_vorticity_3d contains nested
            // functions
            make_vorticity(ARLIM_3D(tileBox.loVect()),
                           ARLIM_3D(tileBox.hiVect()),
                           ARLIM_3D(domainBox.loVect()),
                           ARLIM_3D(domainBox.hiVect()),
                           BL_TO_FORTRAN_3D(vel_mf[mfi]), ZFILL(dx),
                           BL_TO_FORTRAN_3D(vorticity_mf[mfi]), phys_bc.dataPtr());
        }


    }

    // average down and fill ghost cells
    AverageDown(vorticity,0,1);
    FillPatch(t_old,vorticity,vorticity,vorticity,0,0,1,0,bcs_f);
}

void
Maestro::MakeDeltaGamma (const Vector<MultiFab>& state,
                         const RealVector& p0,
                         const Vector<MultiFab>& p0_cart,
                         const RealVector& gamma1bar,
                         const Vector<MultiFab>& gamma1bar_cart,
                         Vector<MultiFab>& deltagamma)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeDeltaGamma()",MakeDeltaGamma);

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& state_mf = state[lev];
        MultiFab& deltagamma_mf = deltagamma[lev];

        const MultiFab& p0cart_mf = p0_cart[lev];
        const MultiFab& gamma1barcart_mf = gamma1bar_cart[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(state_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
#pragma gpu box(tileBox)
            make_deltagamma(AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),
                            BL_TO_FORTRAN_ANYD(state_mf[mfi]), state_mf.nComp(),
                            BL_TO_FORTRAN_ANYD(p0cart_mf[mfi]), BL_TO_FORTRAN_ANYD(gamma1barcart_mf[mfi]),
                            BL_TO_FORTRAN_ANYD(deltagamma_mf[mfi]));

        }
    }

    // average down and fill ghost cells
    AverageDown(deltagamma,0,1);
    FillPatch(t_old,deltagamma,deltagamma,deltagamma,0,0,1,0,bcs_f);
}

void
Maestro::MakeEntropy (const Vector<MultiFab>& state,
                      Vector<MultiFab>& entropy)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeEntropy()",MakeEntropy);

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& state_mf = state[lev];
        MultiFab& entropy_mf = entropy[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(state_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
#pragma gpu box(tileBox)
            make_entropy(AMREX_INT_ANYD(tileBox.loVect()),
                         AMREX_INT_ANYD(tileBox.hiVect()),lev,
                         BL_TO_FORTRAN_ANYD(state_mf[mfi]), state_mf.nComp(),
                         BL_TO_FORTRAN_ANYD(entropy_mf[mfi]));
        }

    }

    // average down and fill ghost cells
    AverageDown(entropy,0,1);
    FillPatch(t_old,entropy,entropy,entropy,0,0,1,0,bcs_f);
}

void
Maestro::MakeDivw0 (const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& w0mac,
                    Vector<MultiFab>& divw0)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeDivw0()",MakeDivw0);

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& divw0_mf = divw0[lev];
        const MultiFab& w0_mf = w0_cart[lev];

        if (spherical == 0) {

            // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
            for ( MFIter mfi(divw0_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

                // Get the index space of the valid region
                const Box& tileBox = mfi.tilebox();
                const Real* dx = geom[lev].CellSize();

                // call fortran subroutine
                // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
                // lo/hi coordinates (including ghost cells), and/or the # of components
                // We will also pass "validBox", which specifies the "valid" region.
#pragma gpu box(tileBox)
                make_divw0(AMREX_INT_ANYD(tileBox.loVect()),
                           AMREX_INT_ANYD(tileBox.hiVect()),
                           BL_TO_FORTRAN_ANYD(w0_mf[mfi]), w0_mf.nComp(),
                           AMREX_REAL_ANYD(dx),
                           BL_TO_FORTRAN_ANYD(divw0_mf[mfi]));
            }

        } else {

            const MultiFab& w0macx_mf = w0mac[lev][0];
            const MultiFab& w0macy_mf = w0mac[lev][1];
            const MultiFab& w0macz_mf = w0mac[lev][2];

            // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
            for ( MFIter mfi(divw0_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

                // Get the index space of the valid region
                const Box& tileBox = mfi.tilebox();
                const Real* dx = geom[lev].CellSize();

                // call fortran subroutine
                // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
                // lo/hi coordinates (including ghost cells), and/or the # of components
                // We will also pass "validBox", which specifies the "valid" region.
#pragma gpu box(tileBox)
                make_divw0_sphr(AMREX_INT_ANYD(tileBox.loVect()),
                                AMREX_INT_ANYD(tileBox.hiVect()),
                                BL_TO_FORTRAN_ANYD(w0macx_mf[mfi]),
                                BL_TO_FORTRAN_ANYD(w0macy_mf[mfi]),
                                BL_TO_FORTRAN_ANYD(w0macz_mf[mfi]),
                                AMREX_REAL_ANYD(dx),
                                BL_TO_FORTRAN_ANYD(divw0_mf[mfi]));
            }
        }
    }

    // average down and fill ghost cells
    AverageDown(divw0,0,1);
    FillPatch(t_old,divw0,divw0,divw0,0,0,1,0,bcs_f);
}

void
Maestro::MakePiDivu (const Vector<MultiFab>& vel,
                     const Vector<MultiFab>& state,
                     Vector<MultiFab>& pidivu)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakePiDivu()",MakePiDivu);

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& vel_mf = vel[lev];
        const MultiFab& state_mf = state[lev];
        MultiFab& pidivu_mf = pidivu[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(pidivu_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const Real* dx = geom[lev].CellSize();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
#pragma gpu box(tileBox)
            make_pidivu(AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),
                        BL_TO_FORTRAN_ANYD(vel_mf[mfi]),
                        AMREX_REAL_ANYD(dx),
                        BL_TO_FORTRAN_ANYD(state_mf[mfi]), state_mf.nComp(),
                        BL_TO_FORTRAN_ANYD(pidivu_mf[mfi]));
        }

    }

    // average down and fill ghost cells
    AverageDown(pidivu,0,1);
    FillPatch(t_old,pidivu,pidivu,pidivu,0,0,1,0,bcs_f);
}

void
Maestro::MakeAbar (const Vector<MultiFab>& state,
                   Vector<MultiFab>& abar)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeAbar()",MakeAbar);

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& state_mf = state[lev];
        MultiFab& abar_mf = abar[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(abar_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
#pragma gpu box(tileBox)
            make_abar(AMREX_INT_ANYD(tileBox.loVect()),
                      AMREX_INT_ANYD(tileBox.hiVect()),
                      BL_TO_FORTRAN_ANYD(state_mf[mfi]), state_mf.nComp(),
                      BL_TO_FORTRAN_ANYD(abar_mf[mfi]));
        }

    }

    // average down and fill ghost cells
    AverageDown(abar,0,1);
    FillPatch(t_old,abar,abar,abar,0,0,1,0,bcs_f);
}

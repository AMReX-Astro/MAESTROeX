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

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(vel[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const Array4<const Real> vel_arr = vel[lev].array(mfi);
            const Array4<Real> magvel_arr = magvel[lev].array(mfi);

            if (spherical == 0) {
                const Array4<const Real> w0_arr = w0_cart[lev].array(mfi);

                AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
#if (AMREX_SPACEDIM == 2)
                    Real vertical_vel = vel_arr(i,j,k,1) 
                        + 0.5*(w0_arr(i,j,k,1) + w0_arr(i,j+1,k,1));

                    magvel_arr(i,j,k) = sqrt(vel_arr(i,j,k,0)*vel_arr(i,j,k,0) 
                        + vertical_vel*vertical_vel);
#elif (AMREX_SPACEDIM == 3)
                    Real vertical_vel = vel_arr(i,j,k,2) 
                        + 0.5*(w0_arr(i,j,k,2) + w0_arr(i,j,k+1,2));

                    magvel_arr(i,j,k) = sqrt(vel_arr(i,j,k,0)*vel_arr(i,j,k,0) 
                        + vel_arr(i,j,k,1)*vel_arr(i,j,k,1)  
                        + vertical_vel*vertical_vel);
#endif
                });
            } else {
                const Array4<const Real> w0macx = w0mac[lev][0].array(mfi);
                const Array4<const Real> w0macy = w0mac[lev][1].array(mfi);
                const Array4<const Real> w0macz = w0mac[lev][2].array(mfi);

                AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                    Real velx = vel_arr(i,j,k,0) 
                        + 0.5*(w0macx(i,j,k) + w0macx(i+1,j,k));
                    Real vely = vel_arr(i,j,k,1) 
                        + 0.5*(w0macy(i,j,k) + w0macy(i,j+1,k));
                    Real velz = vel_arr(i,j,k,2) 
                        + 0.5*(w0macz(i,j,k) + w0macz(i,j,k+1));

                    magvel_arr(i,j,k) = sqrt(velx*velx 
                        + vely*vely + velz*velz);
                });
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

#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(vel[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const Array4<const Real> vel_arr = vel[lev].array(mfi);
            const Array4<Real> radvel_arr = rad_vel[lev].array(mfi);
            const Array4<Real> circvel_arr = circ_vel[lev].array(mfi);
            const Array4<const Real> w0rcart_arr = w0rcart[lev].array(mfi);
            const Array4<const Real> normal_arr = normal[lev].array(mfi);

            AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                circvel_arr(i,j,k) = 0.0;
                radvel_arr(i,j,k) = 0.0;

                for (int n=0; n<AMREX_SPACEDIM; ++n) {
                    radvel_arr(i,j,k) += vel_arr(i,j,k,n) * normal_arr(i,j,k,n);
                }

                for (int n=0; n<AMREX_SPACEDIM; ++n) {
                    Real v = vel_arr(i,j,k,n) - radvel_arr(i,j,k)*normal_arr(i,j,k,n);
                    circvel_arr(i,j,k) += v*v;
                }

                circvel_arr(i,j,k) = sqrt(circvel_arr(i,j,k));

                // add base state vel to get full radial velocity
                radvel_arr(i,j,k) += w0rcart_arr(i,j,k);
            });
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

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(divw0[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const auto dx = geom[lev].CellSizeArray();
            const Array4<Real> divw0_arr = divw0[lev].array(mfi);

            if (spherical == 0) {

                const Array4<const Real> w0_arr = w0_cart[lev].array(mfi);

                AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
#if (AMREX_SPACEDIM == 2)
                    divw0_arr(i,j,k) = (w0_arr(i,j+1,k,1) - w0_arr(i,j,k,1)) / dx[1];
#else
                    divw0_arr(i,j,k) = (w0_arr(i,j,k+1,2) - w0_arr(i,j,k,2)) / dx[2];
#endif                    
                });

            } else {

                const Array4<const Real> w0macx = w0mac[lev][0].array(mfi);
                const Array4<const Real> w0macy = w0mac[lev][1].array(mfi);
                const Array4<const Real> w0macz = w0mac[lev][2].array(mfi);

                AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                    divw0_arr(i,j,k) = (w0macx(i+1,j,k) - w0macx(i,j,k)) / dx[0]
                    + (w0macy(i,j+1,k) - w0macy(i,j,k)) / dx[1]
                    + (w0macz(i,j,k+1) - w0macz(i,j,k)) / dx[2];    
                });
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

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(pidivu[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const auto dx = geom[lev].CellSizeArray();

            const Array4<const Real> vel_arr = vel[lev].array(mfi);
            const Array4<const Real> pi_cc = state[lev].array(mfi, Pi);
            const Array4<Real> pidivu_arr = pidivu[lev].array(mfi);

            AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                pidivu_arr(i,j,k) = pi_cc(i,j,k) * 0.5 * (vel_arr(i+1,j,k,0) - vel_arr(i-1,j,k,0)) / dx[0];
                pidivu_arr(i,j,k) += pi_cc(i,j,k) * 0.5 * (vel_arr(i,j+1,k,1) - vel_arr(i,j-1,k,1)) / dx[1];
#if (AMREX_SPACEDIM == 3)
                pidivu_arr(i,j,k) += pi_cc(i,j,k) * 0.5 * (vel_arr(i,j,k+1,2) - vel_arr(i,j,k-1,2)) / dx[2];
#endif
            });
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


#include <Maestro.H>
#include <Maestro_F.H>

#ifdef AMREX_USE_CUDA
#include <cuda_runtime_api.h>
#include <AMReX_Arena.H>
#endif

using namespace amrex;

void
Maestro::EstDt ()
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::EstDt()", EstDt);

    dt = 1.e20;

    // allocate a dummy w0_force and set equal to zero
    RealVector w0_force_dummy( (max_radial_level+1)*nr_fine );
    w0_force_dummy.shrink_to_fit();
    std::fill(w0_force_dummy.begin(),w0_force_dummy.end(), 0.);

    // build dummy w0_force_cart and set equal to zero
    Vector<MultiFab> w0_force_cart_dummy(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        w0_force_cart_dummy[lev].define(grids[lev], dmap[lev], 3, 1);
        w0_force_cart_dummy[lev].setVal(0.);
    }

    // build a dummy umac and set equal to zero
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > umac_dummy(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        umac_dummy[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);
        umac_dummy[lev][0].setVal(0.);
        umac_dummy[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
        umac_dummy[lev][1].setVal(0.);
#if (AMREX_SPACEDIM == 3)
        umac_dummy[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1);
        umac_dummy[lev][2].setVal(0.);
#endif
    }

    // face-centered
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > w0mac(finest_level+1);

    if (spherical == 1) {
        // initialize
        for (int lev=0; lev<=finest_level; ++lev) {
            w0mac[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);
            w0mac[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
#if (AMREX_SPACEDIM == 3)
            w0mac[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1);
#endif
        }

        for (int lev=0; lev<=finest_level; ++lev) {
            for (int idim=0; idim<AMREX_SPACEDIM; ++idim) {
                w0mac[lev][idim].setVal(0.);
            }
        }

        if (evolve_base_state && (use_exact_base_state == 0 && average_base_state == 0)) {
            MakeW0mac(w0mac);
        }
    }

    // build and compute vel_force
    Vector<MultiFab> vel_force(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        vel_force[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        // needed to avoid NaNs in filling corner ghost cells with 2 physical boundaries
        vel_force[lev].setVal(0.);
    }

    int do_add_utilde_force = 0;
    MakeVelForce(vel_force,umac_dummy,sold,rho0_old,grav_cell_old,
                 w0_force_cart_dummy,do_add_utilde_force);

#if (AMREX_SPACEDIM == 3)
    // build and initialize grad_p0 for spherical case
    Vector<MultiFab> gp0_cart(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        gp0_cart[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        gp0_cart[lev].setVal(0.);
    }
    RealVector gp0( (max_radial_level+1)*(nr_fine+1) );
    gp0.shrink_to_fit();
    std::fill(gp0.begin(),gp0.end(), 0.);

    // divU constraint
    if (spherical == 1) {
      estdt_divu(gp0.dataPtr(), p0_old.dataPtr(), gamma1bar_old.dataPtr(),
                 r_cc_loc.dataPtr(), r_edge_loc.dataPtr());
    }
      
    Put1dArrayOnCart (gp0,gp0_cart,1,1,bcs_f,0);
#endif

    Vector<MultiFab> p0_cart(finest_level+1);
    Vector<MultiFab> gamma1bar_cart(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        p0_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        gamma1bar_cart[lev].define(grids[lev], dmap[lev], 1, 1);
    }

    Put1dArrayOnCart(p0_old,p0_cart,0,0,bcs_f,0);
    Put1dArrayOnCart(gamma1bar_old,gamma1bar_cart,0,0,bcs_f,0);

    Real umax = 0.;

    Real dt_lev = 1.e99;
    Real umax_lev = 0.;

    for (int lev = 0; lev <= finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& uold_mf = uold[lev];
        MultiFab& sold_mf = sold[lev];
        MultiFab& vel_force_mf = vel_force[lev];
        MultiFab& S_cc_old_mf = S_cc_old[lev];
        MultiFab& dSdt_mf = dSdt[lev];
#if (AMREX_SPACEDIM == 3)
        MultiFab& w0macx_mf = w0mac[lev][0];
        MultiFab& w0macy_mf = w0mac[lev][1];
        MultiFab& w0macz_mf = w0mac[lev][2];
        const MultiFab& gp0_cart_mf = gp0_cart[lev];
#endif
        const MultiFab& w0_mf = w0_cart[lev];
        const MultiFab& p0_mf = p0_cart[lev];
        const MultiFab& gamma1bar_mf = gamma1bar_cart[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel reduction(min:dt_lev) reduction(max:umax_lev)
#endif
        {

        Real dt_grid = 1.e99;
        Real umax_grid = 0.;

        for (MFIter mfi(uold_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const auto dx = geom[lev].CellSizeArray();

            const Array4<const Real> scal_arr = sold[lev].array(mfi);
            const Array4<const Real> u = uold[lev].array(mfi);
            const Array4<const Real> force = vel_force[lev].array(mfi);
            const Array4<const Real> S_cc_arr = S_cc_old[lev].array(mfi);
            const Array4<const Real> dSdt_arr = dSdt[lev].array(mfi);
            const Array4<const Real> w0_arr = w0_cart[lev].array(mfi);
            const Array4<const Real> p0_arr = p0_cart[lev].array(mfi);
            const Array4<const Real> gamma1bar_arr = gamma1bar_cart[lev].array(mfi);

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            if (spherical == 0) {

                const Real rho_min = 1.e-20;
                Real dt_temp = 1.e99;
                const Real eps = 1.e-8;

                Real spdx = 0.0;
                Real spdy = 0.0;
                Real spdz = 0.0;
                Real spdr = 0.0;

                // Limit dt based on velocity terms.
                AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                    spdx = amrex::max(spdx, fabs(u(i,j,k,0)));
#if (AMREX_SPACEDIM == 2)
                    spdy = amrex::max(spdy, fabs(u(i,j,k,1) + 0.5 * (w0_arr(i,j,k,1)+w0_arr(i,j+1,k,1))));
#else   
                    spdy = amrex::max(spdy, fabs(u(i,j,k,1)));
                    spdz = amrex::max(spdz, fabs(u(i,j,k,2) + 0.5 * (w0_arr(i,j,k,2)+w0_arr(i,j,k+1,2))));
#endif
                    spdr = amrex::max(spdr, fabs(w0_arr(i,j,k,AMREX_SPACEDIM-1)));

                    umax_grid = amrex::max(umax_grid, std::max(spdx,std::max(spdy,std::max(spdz,spdr))));
                });

                if (spdx > eps) dt_temp = amrex::min(dt_temp, dx[0] / spdx);
                if (spdy > eps) dt_temp = amrex::min(dt_temp, dx[1] / spdy);
                if (spdz > eps) dt_temp = amrex::min(dt_temp, dx[2] / spdz);
                if (spdr > eps) dt_temp = amrex::min(dt_temp, dx[AMREX_SPACEDIM-1] / spdr);

                dt_temp *= cfl;

                // Limit dt based on forcing terms
                Real fx = 0.0;
                Real fy = 0.0;
                Real fz = 0.0;

                AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                    fx = amrex::max(fx, fabs(force(i,j,k,0)));
                    fy = amrex::max(fy, fabs(force(i,j,k,1)));
#if (AMREX_SPACEDIM == 3)
                    fz = amrex::max(fz, fabs(force(i,j,k,2)));
#endif
                });

                if (fx > eps) dt_temp = amrex::min(dt_temp, std::sqrt(2.0 * dx[0] / fx));
                if (fy > eps) dt_temp = amrex::min(dt_temp, std::sqrt(2.0 * dx[1] / fy));
#if (AMREX_SPACEDIM == 3)
                if (fz > eps) dt_temp = amrex::min(dt_temp, std::sqrt(2.0 * dx[2] / fz));
#endif

                const auto nr_lev = nr[lev];

                // divU constraint
                AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                    Real gradp0 = 0.0;
#if (AMREX_SPACEDIM == 2)
                    if (j == 0) {
                        gradp0 = (p0_arr(i,j+1,k) - p0_arr(i,j,k)) / dx[1];
                    } else if (j == nr_lev-1) {
                        gradp0 = (p0_arr(i,j,k) - p0_arr(i,j-1,k)) / dx[1];
                    } else {
                        gradp0 = 0.5 * (p0_arr(i,j+1,k) - p0_arr(i,j-1,k)) / dx[1];
                    }
#else 
                    if (k == 0) {
                        gradp0 = (p0_arr(i,j,k+1) - p0_arr(i,j,k)) / dx[2];
                    } else if (k == nr_lev-1) {
                        gradp0 = (p0_arr(i,j,k) - p0_arr(i,j,k-1)) / dx[2];
                    } else {
                        gradp0 = 0.5 * (p0_arr(i,j,k+1) - p0_arr(i,j,k-1)) / dx[2];
                    }
#endif
                    Real denom = S_cc_arr(i,j,k) - u(i,j,k,AMREX_SPACEDIM-1) * gradp0 / (gamma1bar_arr(i,j,k)*p0_arr(i,j,k));

                    if (denom > 0.0) {
                        dt_temp = amrex::min(dt_temp, 0.4*(1.0 - rho_min / scal_arr(i,j,k,Rho)) / denom);
                    }
                });

                // An additional dS/dt timestep constraint originally
                // used in nova
                // solve the quadratic equation
                // (rho - rho_min)/(rho dt) = S + (dt/2)*(dS/dt)
                // which is equivalent to
                // (rho/2)*dS/dt*dt^2 + rho*S*dt + (rho_min-rho) = 0
                // which has solution dt = 2.0d0*c/(-b-sqrt(b**2-4.0d0*a*c))
                AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {

                    if (dSdt_arr(i,j,k) > 1.e-20) {
                        auto a = 0.5 * scal_arr(i,j,k,Rho) * dSdt_arr(i,j,k);
                        auto b = scal_arr(i,j,k,Rho) * S_cc_arr(i,j,k);
                        auto c = rho_min - scal_arr(i,j,k,Rho);
                        dt_temp = amrex::min(dt_temp, 0.4*2.0*c / (-b-std::sqrt(b*b-4.0*a*c)));
                    }
                });

                dt_grid = amrex::min(dt_grid, dt_temp);


// #pragma gpu box(tileBox)
//                 estdt(lev, AMREX_MFITER_REDUCE_MIN(&dt_grid),
//                       AMREX_MFITER_REDUCE_MAX(&umax_grid),
//                       AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),
//                       AMREX_REAL_ANYD(dx),
//                       BL_TO_FORTRAN_ANYD(sold_mf[mfi]), sold_mf[mfi].nCompPtr(),
//                       BL_TO_FORTRAN_ANYD(uold_mf[mfi]), uold_mf[mfi].nCompPtr(),
//                       BL_TO_FORTRAN_ANYD(vel_force_mf[mfi]), vel_force_mf[mfi].nCompPtr(),
//                       BL_TO_FORTRAN_ANYD(S_cc_old_mf[mfi]),
//                       BL_TO_FORTRAN_ANYD(dSdt_mf[mfi]),
//                       BL_TO_FORTRAN_ANYD(w0_mf[mfi]),
//                       BL_TO_FORTRAN_ANYD(p0_mf[mfi]),
//                       BL_TO_FORTRAN_ANYD(gamma1bar_mf[mfi]));
            } else {
#if (AMREX_SPACEDIM == 3)

// #pragma gpu box(tileBox)
//                 estdt_sphr(AMREX_MFITER_REDUCE_MIN(&dt_grid),
//                            AMREX_MFITER_REDUCE_MAX(&umax_grid),
//                            AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),
//                            AMREX_REAL_ANYD(dx),
//                            BL_TO_FORTRAN_ANYD(sold_mf[mfi]), sold_mf[mfi].nCompPtr(),
//                            BL_TO_FORTRAN_ANYD(uold_mf[mfi]), uold_mf[mfi].nCompPtr(),
//                            BL_TO_FORTRAN_ANYD(vel_force_mf[mfi]), vel_force_mf[mfi].nCompPtr(),
//                            BL_TO_FORTRAN_ANYD(S_cc_old_mf[mfi]),
//                            BL_TO_FORTRAN_ANYD(dSdt_mf[mfi]),
//                BL_TO_FORTRAN_ANYD(w0_mf[mfi]),
//                            BL_TO_FORTRAN_ANYD(w0macx_mf[mfi]),
//                            BL_TO_FORTRAN_ANYD(w0macy_mf[mfi]),
//                            BL_TO_FORTRAN_ANYD(w0macz_mf[mfi]),
//                            BL_TO_FORTRAN_ANYD(gp0_cart_mf[mfi]));

#else
                Abort("EstDt: Spherical is not valid for DIM < 3");
#endif
            }
        }

        dt_lev = std::min(dt_lev,dt_grid);
        umax_lev = std::max(umax_lev,umax_grid);

        } //end openmp

        // find the smallest dt over all processors
        ParallelDescriptor::ReduceRealMin(dt_lev);

        // find the largest umax over all processors
        ParallelDescriptor::ReduceRealMax(umax_lev);

        // update umax over all levels
        umax = std::max(umax,umax_lev);

        if (maestro_verbose > 0) {
            Print() << "Call to estdt for level " << lev << " gives dt_lev = " << dt_lev << std::endl;
        }

        // update dt over all levels
        dt = std::min(dt,dt_lev);

    }     // end loop over levels

    if (maestro_verbose > 0) {
        Print() << "Minimum estdt over all levels = " << dt << std::endl;
    }

    if (dt < small_dt) {
        Abort("EstDt: dt < small_dt");
    }

    if (dt > max_dt) {
        if (maestro_verbose > 0) {
            Print() << "max_dt limits the new dt = " << max_dt << std::endl;
        }
        dt = max_dt;
    }

    if (fixed_dt != -1.0) {
        // fixed dt
        dt = fixed_dt;
        if (maestro_verbose > 0) {
            Print() << "Setting fixed dt = " << dt << std::endl;
        }
    }

    // set rel_eps in fortran module
    umax *= 1.e-8;
    c_rel_eps = umax;
    set_rel_eps(&umax);
}


void
Maestro::FirstDt ()
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::FirstDt()",FirstDt);

    dt = 1.e20;

    // allocate a dummy w0_force and set equal to zero
    RealVector w0_force_dummy( (max_radial_level+1)*nr_fine );
    w0_force_dummy.shrink_to_fit();
    std::fill(w0_force_dummy.begin(),w0_force_dummy.end(), 0.);

    // build dummy w0_force_cart and set equal to zero
    Vector<MultiFab> w0_force_cart_dummy(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        w0_force_cart_dummy[lev].define(grids[lev], dmap[lev], 3, 1);
        w0_force_cart_dummy[lev].setVal(0.);
    }

    // build a dummy umac and set equal to zero
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > umac_dummy(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        umac_dummy[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);
        umac_dummy[lev][0].setVal(0.);
        umac_dummy[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
        umac_dummy[lev][1].setVal(0.);
#if (AMREX_SPACEDIM == 3)
        umac_dummy[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1);
        umac_dummy[lev][2].setVal(0.);
#endif
    }

    // build and compute vel_force
    Vector<MultiFab> vel_force(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        vel_force[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        // needed to avoid NaNs in filling corner ghost cells with 2 physical boundaries
        vel_force[lev].setVal(0.);
    }

    int do_add_utilde_force = 0;
    MakeVelForce(vel_force,umac_dummy,sold,rho0_old,grav_cell_old,
                 w0_force_cart_dummy,do_add_utilde_force);

#if (AMREX_SPACEDIM == 3)
    // build and initialize grad_p0 for spherical case
    Vector<MultiFab> gp0_cart(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        gp0_cart[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        gp0_cart[lev].setVal(0.);
    }
    RealVector gp0( (max_radial_level+1)*(nr_fine+1) );
    gp0.shrink_to_fit();
    std::fill(gp0.begin(),gp0.end(), 0.);

    // divU constraint
    if (use_divu_firstdt && spherical == 1) {
        estdt_divu(gp0.dataPtr(), p0_old.dataPtr(), gamma1bar_old.dataPtr(),
                   r_cc_loc.dataPtr(), r_edge_loc.dataPtr());
    }

    Put1dArrayOnCart (gp0,gp0_cart,1,1,bcs_f,0);
#endif

    Vector<MultiFab> p0_cart(finest_level+1);
    Vector<MultiFab> gamma1bar_cart(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        p0_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        gamma1bar_cart[lev].define(grids[lev], dmap[lev], 1, 1);
    }

    Put1dArrayOnCart(p0_old,p0_cart,0,0,bcs_f,0);
    Put1dArrayOnCart(gamma1bar_old,gamma1bar_cart,0,0,bcs_f,0);

    Real umax = 0.;

    Real dt_lev = 1.e99;
    Real umax_lev = 0.;

    for (int lev = 0; lev <= finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& uold_mf = uold[lev];
        const MultiFab& sold_mf = sold[lev];
        const MultiFab& vel_force_mf = vel_force[lev];
        const MultiFab& S_cc_old_mf = S_cc_old[lev];
#if (AMREX_SPACEDIM == 3)
        const MultiFab& gp0_cart_mf = gp0_cart[lev];
#endif
        const MultiFab& p0_mf = p0_cart[lev];
        const MultiFab& gamma1bar_mf = gamma1bar_cart[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel reduction(min:dt_lev) reduction(max:umax_lev)
#endif
        {

        Real dt_grid = 1.e99;
        Real umax_grid = 0.;

        for ( MFIter mfi(sold_mf,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const Real* dx = geom[lev].CellSize();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            if (spherical == 0 ) {

#pragma gpu box(tileBox)
                firstdt(lev,
                    AMREX_MFITER_REDUCE_MIN(&dt_grid),
                    AMREX_MFITER_REDUCE_MAX(&umax_grid),
                    AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),
                    AMREX_REAL_ANYD(dx),
                    BL_TO_FORTRAN_ANYD(sold_mf[mfi]), sold_mf[mfi].nCompPtr(),
                    BL_TO_FORTRAN_ANYD(uold_mf[mfi]), uold_mf[mfi].nCompPtr(),
                    BL_TO_FORTRAN_ANYD(vel_force_mf[mfi]), vel_force_mf[mfi].nCompPtr(),
                    BL_TO_FORTRAN_ANYD(S_cc_old_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(p0_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(gamma1bar_mf[mfi]));

            } else {
#if (AMREX_SPACEDIM == 3)

#pragma gpu box(tileBox)
                firstdt_sphr(AMREX_MFITER_REDUCE_MIN(&dt_grid),
                             AMREX_MFITER_REDUCE_MAX(&umax_grid),
                             AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),
                             AMREX_REAL_ANYD(dx),
                             BL_TO_FORTRAN_ANYD(sold_mf[mfi]), sold_mf[mfi].nCompPtr(),
                             BL_TO_FORTRAN_ANYD(uold_mf[mfi]), uold_mf[mfi].nCompPtr(),
                             BL_TO_FORTRAN_ANYD(vel_force_mf[mfi]), vel_force_mf[mfi].nCompPtr(),
                             BL_TO_FORTRAN_ANYD(S_cc_old_mf[mfi]),
                             BL_TO_FORTRAN_ANYD(gp0_cart_mf[mfi]));
#else
                Abort("FirstDt: Spherical is not valid for DIM < 3");
#endif
            }
        }

        dt_lev = std::min(dt_lev,dt_grid);
        umax_lev = std::max(umax_lev,umax_grid);

        } //end openmp

        // find the smallest dt over all processors
        ParallelDescriptor::ReduceRealMin(dt_lev);

        // find the largest umax over all processors
        ParallelDescriptor::ReduceRealMax(umax_lev);

        // update umax over all levels
        umax = std::max(umax,umax_lev);

        if (maestro_verbose > 0) {
            Print() << "Call to firstdt for level " << lev << " gives dt_lev = " << dt_lev << std::endl;
        }

        // multiply by init_shrink
        dt_lev *= init_shrink;

        if (maestro_verbose > 0) {
            Print() << "Multiplying dt_lev by init_shrink; dt_lev = " << dt_lev << std::endl;
        }

        // update dt over all levels
        dt = std::min(dt,dt_lev);

    }     // end loop over levels

    if (maestro_verbose > 0) {
        Print() << "Minimum firstdt over all levels = " << dt << std::endl;
    }

    if (dt < small_dt) {
        Abort("FirstDt: dt < small_dt");
    }

    if (dt > max_dt) {
        if (maestro_verbose > 0) {
            Print() << "max_dt limits the new dt = " << max_dt << std::endl;
        }
        dt = max_dt;
    }

    if (fixed_dt != -1.0) {
        // fixed dt
        dt = fixed_dt;
        if (maestro_verbose > 0) {
            Print() << "Setting fixed dt = " << dt << std::endl;
        }
    }

    // set rel_eps in fortran module
    umax *= 1.e-8;
    c_rel_eps = umax;
    set_rel_eps(&umax);
}

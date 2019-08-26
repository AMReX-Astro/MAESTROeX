
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void
Maestro::MakeVelForce (Vector<MultiFab>& vel_force,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& uedge,
                       const Vector<MultiFab>& rho,
                       const RealVector& rho0,
                       const RealVector& grav_cell,
                       const RealVector& w0_force,
                       const Vector<MultiFab>& w0_force_cart,
                       int do_add_utilde_force)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeVelForce()",MakeVelForce);

    // For spherical case
    Vector<MultiFab> gradw0_cart(finest_level+1);
    Vector<MultiFab> grav_cart(finest_level+1);
    Vector<MultiFab> rho0_cart(finest_level+1);

    for (int lev=0; lev<=finest_level; ++lev) {
        gradw0_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        gradw0_cart[lev].setVal(0.);

        grav_cart[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        grav_cart[lev].setVal(0.);

        rho0_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        rho0_cart[lev].setVal(0.);

    }

    if (spherical == 1) {
        RealVector gradw0( (max_radial_level+1)*nr_fine );
        gradw0.shrink_to_fit();

        if (use_exact_base_state || average_base_state) {
            std::fill(gradw0.begin(), gradw0.end(), 0.);
        } else {
            compute_grad_phi_rad(w0.dataPtr(), gradw0.dataPtr());
        }

        Put1dArrayOnCart(gradw0,gradw0_cart,0,0,bcs_u,0);
        Put1dArrayOnCart(rho0,rho0_cart,0,0,bcs_f,0);
        Put1dArrayOnCart(grav_cell,grav_cart,0,1,bcs_f,0);
    }

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& vel_force_mf = vel_force[lev];
        const MultiFab& gpi_mf = gpi[lev];
        const MultiFab& rho_mf = rho[lev];
        const MultiFab& uedge_mf = uedge[lev][0];
        const MultiFab& vedge_mf = uedge[lev][1];
#if (AMREX_SPACEDIM == 3)
        const MultiFab& wedge_mf = uedge[lev][2];
        const MultiFab& gradw0_mf = gradw0_cart[lev];
        const MultiFab& normal_mf = normal[lev];
        const MultiFab& w0force_mf = w0_force_cart[lev];
#endif
        const MultiFab& grav_mf = grav_cart[lev];
        const MultiFab& rho0_mf = rho0_cart[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(vel_force_mf, true); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const Real* dx = geom[lev].CellSize();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            if (spherical == 0) {
#pragma gpu box(tileBox)
                make_vel_force(AMREX_INT_ANYD(tileBox.loVect()),
                               AMREX_INT_ANYD(tileBox.hiVect()), lev,
                               BL_TO_FORTRAN_ANYD(vel_force_mf[mfi]),
                               BL_TO_FORTRAN_ANYD(gpi_mf[mfi]),
                               BL_TO_FORTRAN_N_ANYD(rho_mf[mfi],Rho),
                               BL_TO_FORTRAN_ANYD(uedge_mf[mfi]),
                               BL_TO_FORTRAN_ANYD(vedge_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
                               BL_TO_FORTRAN_ANYD(wedge_mf[mfi]),
#endif
                               w0.dataPtr(),
                               w0_force.dataPtr(),
                               rho0.dataPtr(),
                               grav_cell.dataPtr(),
                               do_add_utilde_force);
            } else {

#if (AMREX_SPACEDIM == 3)
#pragma gpu box(tileBox)
                make_vel_force_sphr(AMREX_INT_ANYD(tileBox.loVect()),
                                    AMREX_INT_ANYD(tileBox.hiVect()),
                                    BL_TO_FORTRAN_ANYD(vel_force_mf[mfi]),
                                    BL_TO_FORTRAN_ANYD(gpi_mf[mfi]),
                                    BL_TO_FORTRAN_N_ANYD(rho_mf[mfi],Rho),
                                    BL_TO_FORTRAN_ANYD(uedge_mf[mfi]),
                                    BL_TO_FORTRAN_ANYD(vedge_mf[mfi]),
                                    BL_TO_FORTRAN_ANYD(wedge_mf[mfi]),
                                    BL_TO_FORTRAN_ANYD(normal_mf[mfi]),
                                    BL_TO_FORTRAN_ANYD(gradw0_mf[mfi]),
                                    BL_TO_FORTRAN_ANYD(w0force_mf[mfi]),
                                    BL_TO_FORTRAN_ANYD(rho0_mf[mfi]),
                                    BL_TO_FORTRAN_ANYD(grav_mf[mfi]),
                                    AMREX_REAL_ANYD(dx),
                                    do_add_utilde_force);
#else
                Abort("MakeVelForce: Spherical is not valid for DIM < 3");
#endif
            }
        }
    }

    // average fine data onto coarser cells
    AverageDown(vel_force,0,AMREX_SPACEDIM);

    // note - we need to reconsider the bcs type here
    // it matches fortran MAESTRO but is that correct?
    FillPatch(t_old, vel_force, vel_force, vel_force, 0, 0, AMREX_SPACEDIM, 0,
              bcs_u, 1);

}


void
Maestro::ModifyScalForce(Vector<MultiFab>& scal_force,
                         const Vector<MultiFab>& state,
                         const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
                         const RealVector& s0,
                         const RealVector& s0_edge,
                         const Vector<MultiFab>& s0_cart,
                         int comp,
                         const Vector<BCRec>& bcs,
                         int fullform)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::ModifyScalForce()",ModifyScalForce);

#ifdef AMREX_USE_CUDA
    auto not_launched = Gpu::notInLaunchRegion();
    // turn on GPU
    if (not_launched) Gpu::setLaunchRegion(true);
#endif

    RealVector divu;
    Vector<MultiFab> divu_cart(finest_level+1);

    if (spherical == 1) {
        divu.resize(nr_fine);
        std::fill(divu.begin(), divu.end(), 0.);

        if (!use_exact_base_state) {
            for (int r=0; r<nr_fine-1; ++r) {
                Real dr = r_edge_loc[r+1] - r_edge_loc[r];
                divu[r] = (r_edge_loc[r+1]*r_edge_loc[r+1] * w0[r+1]
                           - r_edge_loc[r]*r_edge_loc[r] * w0[r]) / (dr * r_cc_loc[r]*r_cc_loc[r]);
            }
        }

        for (int lev=0; lev<=finest_level; ++lev) {
            divu_cart[lev].define(grids[lev], dmap[lev], 1, 0);
            divu_cart[lev].setVal(0.);
        }
        Put1dArrayOnCart(divu,divu_cart,0,0,bcs_u,0);
    }

    for (int lev=0; lev<=finest_level; ++lev) {

        // Get the index space and grid spacing of the domain
        const Box& domainBox = geom[lev].Domain();
        const Real* dx = geom[lev].CellSize();

        // get references to the MultiFabs at level lev
        MultiFab& scal_force_mf = scal_force[lev];
        const MultiFab& state_mf = state[lev];
        const MultiFab& umac_mf = umac[lev][0];
        const MultiFab& vmac_mf = umac[lev][1];
#if (AMREX_SPACEDIM == 3)
        const MultiFab& wmac_mf = umac[lev][2];
        const MultiFab& s0cart_mf = s0_cart[lev];
        const MultiFab& divu_mf = divu_cart[lev];
#endif

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(scal_force_mf, true); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            if (spherical == 1) {
#if (AMREX_SPACEDIM == 3)
#pragma gpu box(tileBox)
                modify_scal_force_sphr(AMREX_INT_ANYD(tileBox.loVect()),
                                       AMREX_INT_ANYD(tileBox.hiVect()),
                                       AMREX_INT_ANYD(domainBox.loVect()), AMREX_INT_ANYD(domainBox.hiVect()),
                                       scal_force_mf[mfi].dataPtr(comp),
                                       AMREX_INT_ANYD(scal_force_mf[mfi].loVect()), AMREX_INT_ANYD(scal_force_mf[mfi].hiVect()),
                                       state_mf[mfi].dataPtr(comp),
                                       AMREX_INT_ANYD(state_mf[mfi].loVect()), AMREX_INT_ANYD(state_mf[mfi].hiVect()),
                                       BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
                                       BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
                                       BL_TO_FORTRAN_ANYD(wmac_mf[mfi]),
                                       BL_TO_FORTRAN_ANYD(s0cart_mf[mfi]),
                                       w0.dataPtr(), AMREX_REAL_ANYD(dx), fullform,
                                       BL_TO_FORTRAN_ANYD(divu_mf[mfi]));
#else
                Abort("ModifyScalForce: Spherical is not valid for DIM < 3");
#endif
            } else {
#pragma gpu box(tileBox)
                modify_scal_force(AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),lev,
                                  scal_force_mf[mfi].dataPtr(comp),
                                  AMREX_INT_ANYD(scal_force_mf[mfi].loVect()), AMREX_INT_ANYD(scal_force_mf[mfi].hiVect()),
                                  state_mf[mfi].dataPtr(comp),
                                  AMREX_INT_ANYD(state_mf[mfi].loVect()), AMREX_INT_ANYD(state_mf[mfi].hiVect()),
                                  BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
                                  BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
                                  BL_TO_FORTRAN_ANYD(wmac_mf[mfi]),
#endif
                                  s0.dataPtr(), s0_edge.dataPtr(), w0.dataPtr(),
                                  AMREX_REAL_ANYD(dx), fullform);
            }
        }
    }

    // average fine data onto coarser cells
    AverageDown(scal_force,comp,1);

    // fill ghost cells
    FillPatch(t_old, scal_force, scal_force, scal_force, comp, comp, 1, 0, bcs_f);

#ifdef AMREX_USE_CUDA
    // turn off GPU
    if (not_launched) Gpu::setLaunchRegion(false);
#endif
}

void
Maestro::MakeRhoHForce(Vector<MultiFab>& scal_force,
                       int is_prediction,
                       const Vector<MultiFab>& thermal,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
                       int add_thermal,
                       const int &which_step)

{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeRhoHForce()",MakeRhoHForce);

#ifdef AMREX_USE_CUDA
    auto not_launched = Gpu::notInLaunchRegion();
    // turn on GPU
    if (not_launched) Gpu::setLaunchRegion(true);
#endif

    // if we are doing the prediction, then it only makes sense to be in
    // this routine if the quantity we are predicting is rhoh', h, or rhoh
    if (is_prediction == 1 && !(enthalpy_pred_type == predict_rhohprime ||
                                enthalpy_pred_type == predict_h ||
                                enthalpy_pred_type == predict_rhoh) ) {
        Abort("ERROR: should only call mkrhohforce when predicting rhoh', h, or rhoh");
    }

    RealVector rho0( (max_radial_level+1)*nr_fine );
    RealVector   p0( (max_radial_level+1)*nr_fine );
    RealVector grav( (max_radial_level+1)*nr_fine );
    rho0.shrink_to_fit();
    p0.shrink_to_fit();
    grav.shrink_to_fit();

    if (which_step == 1) {
        rho0 = rho0_old;
        p0 =   p0_old;
    }
    else {
        for(int i=0; i<rho0.size(); ++i) {
            rho0[i] = 0.5*(rho0_old[i]+rho0_new[i]);
            p0[i] = 0.5*(  p0_old[i]+  p0_new[i]);
        }
    }

    Vector<MultiFab> p0_cart(finest_level+1);
    Vector<MultiFab> psi_cart(finest_level+1);
    Vector< std::array< MultiFab, AMREX_SPACEDIM > > p0mac(finest_level+1);
    if (spherical == 1) {
        for (int lev=0; lev<=finest_level; ++lev) {
            p0_cart[lev].define(grids[lev], dmap[lev], 1, 1);
            psi_cart[lev].define(grids[lev], dmap[lev], 1, 1);
            AMREX_D_TERM(p0mac[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1); ,
                         p0mac[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1); ,
                         p0mac[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1); );
            psi_cart[lev].setVal(0.);
        }

        Put1dArrayOnCart(p0, p0_cart, 0, 0, bcs_f, 0);
        MakeS0mac(p0, p0mac);
        Put1dArrayOnCart(psi,psi_cart,0,0,bcs_f,0);
    }

    make_grav_cell(grav.dataPtr(),
                   rho0.dataPtr(),
                   r_cc_loc.dataPtr(),
                   r_edge_loc.dataPtr());

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& scal_force_mf = scal_force[lev];
        const MultiFab& umac_mf = umac[lev][0];
        const MultiFab& vmac_mf = umac[lev][1];
#if (AMREX_SPACEDIM == 3)
        const MultiFab& wmac_mf = umac[lev][2];
        const MultiFab& p0cart_mf = p0_cart[lev];
        const MultiFab& p0macx_mf = p0mac[lev][0];
        const MultiFab& p0macy_mf = p0mac[lev][1];
        const MultiFab& p0macz_mf = p0mac[lev][2];
#endif
        const MultiFab& thermal_mf = thermal[lev];
        const MultiFab& psi_mf = psi_cart[lev];

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(scal_force_mf, true); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const Real* dx = geom[lev].CellSize();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            if (spherical == 1) {
#if (AMREX_SPACEDIM == 3)
                // if use_exact_base_state or average_base_state,
                // psi is set to dpdt in advance subroutine
#pragma gpu box(tileBox)
                mkrhohforce_sphr(AMREX_INT_ANYD(tileBox.loVect()),
                                 AMREX_INT_ANYD(tileBox.hiVect()),
                                 scal_force_mf[mfi].dataPtr(RhoH),
                                 AMREX_INT_ANYD(scal_force_mf[mfi].loVect()),
                                 AMREX_INT_ANYD(scal_force_mf[mfi].hiVect()),
                                 BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
                                 BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
                                 BL_TO_FORTRAN_ANYD(wmac_mf[mfi]),
                                 BL_TO_FORTRAN_ANYD(thermal_mf[mfi]),
                                 BL_TO_FORTRAN_ANYD(p0cart_mf[mfi]),
                                 BL_TO_FORTRAN_ANYD(p0macx_mf[mfi]),
                                 BL_TO_FORTRAN_ANYD(p0macy_mf[mfi]),
                                 BL_TO_FORTRAN_ANYD(p0macz_mf[mfi]),
                                 BL_TO_FORTRAN_ANYD(psi_mf[mfi]),
                                 AMREX_REAL_ANYD(dx),
                                 is_prediction, add_thermal);
#else
                Abort("MakeRhoHForce: Spherical is not valid for DIM < 3");
#endif
            } else {

                // if average_base_state, psi is set to dpdt in advance subroutine
#pragma gpu box(tileBox)
                mkrhohforce(AMREX_INT_ANYD(tileBox.loVect()),
                            AMREX_INT_ANYD(tileBox.hiVect()),
                            lev,
                            scal_force_mf[mfi].dataPtr(RhoH),
                            AMREX_INT_ANYD(scal_force_mf[mfi].loVect()),
                            AMREX_INT_ANYD(scal_force_mf[mfi].hiVect()),
#if (AMREX_SPACEDIM == 2)
                            BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
#elif (AMREX_SPACEDIM == 3)
                            BL_TO_FORTRAN_ANYD(wmac_mf[mfi]),
#endif
                            BL_TO_FORTRAN_ANYD(thermal_mf[mfi]),
                            p0.dataPtr(), rho0.dataPtr(), grav.dataPtr(), psi.dataPtr(),
                            is_prediction, add_thermal);
            }
        }
    }

    // average down and fill ghost cells
    AverageDown(scal_force,RhoH,1);
    FillPatch(t_old,scal_force,scal_force,scal_force,RhoH,RhoH,1,0,bcs_f);

#ifdef AMREX_USE_CUDA
    // turn off GPU
    if (not_launched) Gpu::setLaunchRegion(false);
#endif
}

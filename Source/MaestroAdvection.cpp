
#include <Maestro.H>
#include <Maestro_F.H>
#include <MaestroBCThreads.H>

using namespace amrex;

const int predict_rhoh             = 0;
const int predict_rhohprime        = 1;
const int predict_h                = 2;
const int predict_T_then_rhohprime = 3;
const int predict_T_then_h         = 4;
const int predict_hprime           = 5;
const int predict_Tprime_then_h    = 6;

const int predict_rhoprime_and_X   = 1;
const int predict_rhoX             = 2;
const int predict_rho_and_X        = 3;

// compute unprojected mac velocities
void
Maestro::AdvancePremac (Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
                        const Vector<std::array< MultiFab,AMREX_SPACEDIM > >& w0mac,
                        const RealVector& w0_force,
                        const Vector<MultiFab>& w0_force_cart)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::AdvancePremac()",AdvancePremac);

	// create a uold with filled ghost cells
	Vector<MultiFab> utilde(finest_level+1);
	for (int lev=0; lev<=finest_level; ++lev) {
		utilde[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, ng_adv);
		utilde[lev].setVal(0.);
	}

	FillPatch(t_new, utilde, uold, uold, 0, 0, AMREX_SPACEDIM, 0, bcs_u, 1);

	// create a MultiFab to hold uold + w0
	Vector<MultiFab> ufull(finest_level+1);
	for (int lev=0; lev<=finest_level; ++lev) {
		ufull[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, ng_adv);
                // needed to avoid NaNs in filling corner ghost cells with 2 physical boundaries
                ufull[lev].setVal(0.);
	}

	// create ufull = uold + w0
        for (int lev=0; lev<=finest_level; ++lev) {
            MultiFab::Copy(ufull[lev], w0_cart[lev], 0, 0, AMREX_SPACEDIM, 0);
        }
        // fill ufull ghost cells
        FillPatch(t_old, ufull, ufull, ufull, 0, 0, AMREX_SPACEDIM, 0, bcs_u, 1);
	for (int lev=0; lev<=finest_level; ++lev) {
            MultiFab::Add(ufull[lev],utilde[lev],0,0,AMREX_SPACEDIM,ng_adv);
	}
        
	// create a face-centered MultiFab to hold utrans
	Vector<std::array< MultiFab, AMREX_SPACEDIM > > utrans(finest_level+1);
	for (int lev=0; lev<=finest_level; ++lev) {
		utrans[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);
		utrans[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
#if (AMREX_SPACEDIM == 3)
		utrans[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1);
#endif
		for (int j=0; j < AMREX_SPACEDIM; j++)
			utrans[lev][j].setVal(0.);
	}

	// create utrans
	MakeUtrans(utilde,ufull,utrans,w0mac);

	// create a MultiFab to hold the velocity forcing
	Vector<MultiFab> vel_force(finest_level+1);
	for (int lev=0; lev<=finest_level; ++lev) {
		if (ppm_trace_forces == 0) {
			vel_force[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
		} else {
			// tracing needs more ghost cells
			vel_force[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, ng_s);
		}
		vel_force[lev].setVal(0.);

	}

	int do_add_utilde_force = 1;
	MakeVelForce(vel_force,utrans,sold,rho0_old,grav_cell_old,
	             w0_force_cart,do_add_utilde_force);

	// add w0 to trans velocities
	Addw0 (utrans,w0mac,1.);

	VelPred(utilde,ufull,utrans,umac,w0mac,vel_force);
}


void
Maestro::UpdateScal(const Vector<MultiFab>& stateold,
                    Vector<MultiFab>& statenew,
                    const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& sflux,
                    const Vector<MultiFab>& force,
                    int start_comp, int num_comp,
                    const Vector<MultiFab>& p0_cart)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::UpdateScal()",UpdateScal);

    // Make sure to pass in comp+1 for fortran indexing
    const int startcomp = start_comp + 1;
    const int endcomp = startcomp + num_comp;

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& scalold_mf = stateold[lev];
        MultiFab& scalnew_mf = statenew[lev];
        const MultiFab& sfluxx_mf = sflux[lev][0];
        const MultiFab& sfluxy_mf = sflux[lev][1];
#if (AMREX_SPACEDIM == 3)
        const MultiFab& sfluxz_mf = sflux[lev][2];
#endif
    	const MultiFab& p0cart_mf = p0_cart[lev];
        const MultiFab& force_mf = force[lev];

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(scalold_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const Real* dx = geom[lev].CellSize();

            if (start_comp == RhoH)
            {   // Enthalpy update

                // call fortran subroutine
                // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
                // lo/hi coordinates (including ghost cells), and/or the # of components
                // We will also pass "validBox", which specifies the "valid" region.

#pragma gpu box(tileBox)
                update_rhoh(AMREX_INT_ANYD(tileBox.loVect()), 
                        AMREX_INT_ANYD(tileBox.hiVect()),
                        BL_TO_FORTRAN_ANYD(scalold_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(scalnew_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(sfluxx_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(sfluxy_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
                        BL_TO_FORTRAN_ANYD(sfluxz_mf[mfi]),
#endif
                        BL_TO_FORTRAN_ANYD(force_mf[mfi]),
                        BL_TO_FORTRAN_ANYD(p0cart_mf[mfi]),
                        AMREX_REAL_ANYD(dx), dt,
                        NumSpec);

            } else if (start_comp == FirstSpec) {   // RhoX update

                // call fortran subroutine
                // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
                // lo/hi coordinates (including ghost cells), and/or the # of components
                // We will also pass "validBox", which specifies the "valid" region.
#pragma gpu box(tileBox)
                update_rhoX(AMREX_INT_ANYD(tileBox.loVect()),
                    AMREX_INT_ANYD(tileBox.hiVect()),
                    BL_TO_FORTRAN_ANYD(scalold_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(scalnew_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(sfluxx_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(sfluxy_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
                    BL_TO_FORTRAN_ANYD(sfluxz_mf[mfi]),
#endif
                    BL_TO_FORTRAN_ANYD(force_mf[mfi]),
                    AMREX_REAL_ANYD(dx), dt,
                    startcomp, endcomp);
            } else {
                Abort("Invalid scalar in UpdateScal().");
            } // }
        } // end MFIter loop
    } // end loop over levels

    // synchronize by refluxing and averaging down, starting from the finest_level-1/finest_level pair
    if (reflux_type == 2) {
        for (int lev=finest_level-1; lev>=0; --lev) {
            // update lev based on coarse-fine flux mismatch
            flux_reg_s[lev+1]->Reflux(statenew[lev], 1.0, start_comp, start_comp, num_comp, geom[lev]);
            if (start_comp == FirstSpec) {
                // do the same for density if we updated the species
                flux_reg_s[lev+1]->Reflux(statenew[lev], 1.0, Rho, Rho, 1, geom[lev]);
            }
        }
    }

    // average fine data onto coarser cells
    // fill ghost cells
    AverageDown(statenew,start_comp,num_comp);
    FillPatch(t_old, statenew, statenew, statenew, start_comp, start_comp, 
        num_comp, start_comp, bcs_s);

    // do the same for density if we updated the species
    if (start_comp == FirstSpec) {
        AverageDown(statenew,Rho,1);
        FillPatch(t_old, statenew, statenew, statenew, Rho, Rho, 1, Rho, bcs_s);
    }
}

void
Maestro::UpdateVel (const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
                    const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& uedge,
                    const Vector<MultiFab>& force,
                    const Vector<MultiFab>& sponge,
                    const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& w0mac)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::UpdateVel()",UpdateVel);

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& uold_mf = uold[lev];
        MultiFab& unew_mf = unew[lev];
        const MultiFab& umac_mf   = umac[lev][0];
        const MultiFab& uedgex_mf = uedge[lev][0];
        const MultiFab& vmac_mf   = umac[lev][1];
        const MultiFab& uedgey_mf = uedge[lev][1];
#if (AMREX_SPACEDIM == 3)
        const MultiFab& wmac_mf   = umac[lev][2];
        const MultiFab& uedgez_mf = uedge[lev][2];

        // if spherical == 1
        const MultiFab& w0macx_mf = w0mac[lev][0];
        const MultiFab& w0macy_mf = w0mac[lev][1];
        const MultiFab& w0macz_mf = w0mac[lev][2];
#endif
        const MultiFab& force_mf = force[lev];
        const MultiFab& sponge_mf = sponge[lev];
        const MultiFab& w0_mf = w0_cart[lev];

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(force_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const Real* dx = geom[lev].CellSize();

	    if (spherical == 0) {
#pragma gpu box(tileBox)
		update_velocity(AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()), lev,
				 BL_TO_FORTRAN_ANYD(uold_mf[mfi]),
				 BL_TO_FORTRAN_ANYD(unew_mf[mfi]),
				 BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
				 BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
				 BL_TO_FORTRAN_ANYD(wmac_mf[mfi]),
#endif
				 BL_TO_FORTRAN_ANYD(uedgex_mf[mfi]),
				 BL_TO_FORTRAN_ANYD(uedgey_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
				 BL_TO_FORTRAN_ANYD(uedgez_mf[mfi]),
#endif
				 BL_TO_FORTRAN_ANYD(force_mf[mfi]),
				 BL_TO_FORTRAN_ANYD(sponge_mf[mfi]),
				 BL_TO_FORTRAN_ANYD(w0_mf[mfi]),
				 AMREX_REAL_ANYD(dx), dt);
	    } else {
#if (AMREX_SPACEDIM == 3)
#pragma gpu box(tileBox)
		update_velocity_sphr(AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),
				      BL_TO_FORTRAN_ANYD(uold_mf[mfi]),
				      BL_TO_FORTRAN_ANYD(unew_mf[mfi]),
				      BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
				      BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
				      BL_TO_FORTRAN_ANYD(wmac_mf[mfi]),
				      BL_TO_FORTRAN_ANYD(uedgex_mf[mfi]),
				      BL_TO_FORTRAN_ANYD(uedgey_mf[mfi]),
				      BL_TO_FORTRAN_ANYD(uedgez_mf[mfi]),
				      BL_TO_FORTRAN_ANYD(force_mf[mfi]),
				      BL_TO_FORTRAN_ANYD(sponge_mf[mfi]),
				      BL_TO_FORTRAN_ANYD(w0macx_mf[mfi]),
				      BL_TO_FORTRAN_ANYD(w0macy_mf[mfi]),
				      BL_TO_FORTRAN_ANYD(w0macz_mf[mfi]),
				      AMREX_REAL_ANYD(dx), dt);
#else
		Abort("UpdateVel: Spherical is not valid for DIM < 3");
#endif
	    }
        } // end MFIter loop
    } // end loop over levels

    // average fine data onto coarser cells
    AverageDown(unew,0,AMREX_SPACEDIM);

    // fill ghost cells
    FillPatch(t_old, unew, unew, unew, 0, 0, AMREX_SPACEDIM, 0, bcs_u, 1);

}

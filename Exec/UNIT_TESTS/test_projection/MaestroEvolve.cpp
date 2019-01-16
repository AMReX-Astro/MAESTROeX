
#include <Maestro.H>
#include <Problem_F.H>

using namespace amrex;

// advance solution to final time
void
Maestro::Evolve ()
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::Evolve()",Evolve);

	Print() << "Calling Evolve()" << std::endl;

	int project_type;
	get_project_type(&project_type);

	// need to make umid somewhere

	if (project_type == 2) {
		// MAC projection, need to initialize some more multifabs.
		// init_mac_velocity
	}

	//--------------------------------------------------------------------------
	// 'pollute' the velocity field by adding the gradient of a scalar
	//--------------------------------------------------------------------------

	if (project_type == 1) {
		// add_grad_scalar

		for (int lev=0; lev<=finest_level; ++lev) {
			MultiFab& gphi_mf = gphi[lev];
			MultiFab& umid_mf = umid[lev];
			const Real* dx = geom[lev].CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
			for (MFIter mfi(gphi_mf), true); mfi.isValid(); ++mfi)
			{
				const Box& tilebox = mfi.tilebox();
				const int* lo  = tilebox.loVect();
				const int* hi  = tilebox.hiVect();

				add_grad_scalar(ARLIM_3D(lo), ARLIM_3D(hi),
				                BL_TO_FORTRAN_FAB(gphi_mf[mfi]),
				                BL_TO_FORTRAN_FAB(umid_mf[mfi]),
				                phys_bc.dataPtr(), ZFILL(dx));

			}
		}

		// copy umid to unew
	} else {
		for (int lev=0; lev<=finest_level; ++lev) {
			MultiFab& gphi_mac_mf = gphi_mac[lev];
			MultiFab& umac_mid_mf = umac_mid[lev];
			const Real* dx = geom[lev].CellSize();
#ifdef _OPENMP
#pragma omp parallel
#endif
			for (MFIter mfi(gphi_mac_mf, true); mfi.isValid(); ++mfi)
			{
				const Box& tilebox = mfi.tilebox();
				const int* lo  = tilebox.loVect();
				const int* hi  = tilebox.hiVect();

				// NOTE: I have no idea what 'box_phys_bc' is in MAESTROeX, so am
				// guessing it's bcs_u[0]?

				add_grad_scalar_mac(ARLIM_3D(lo), ARLIM_3D(hi),
				                    BL_TO_FORTRAN_FAB(gphi_mac_mf[mfi]),
				                    BL_TO_FORTRAN_FAB(umac_mid_mf[mfi]),
				                    phys_bc.dataPtr(), bcs_u[0].data(), ZFILL(dx));

			}
		}

		for (int lev=0; lev<=finest_level; ++lev) {
			MultiFab& umac_mid_mf = umac_mid[lev];
			MultiFab& utemp_mf = utemp[lev];
			const Real* dx = geom[lev].CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
			for (MFIter mfi(umac_mid_mf, true); mfi.isValid(); ++mfi)
			{
				const Box& tilebox = mfi.tilebox();
				const int* lo  = tilebox.loVect();
				const int* hi  = tilebox.hiVect();

				convert_MAC_to_cc(ARLIM_3D(lo), ARLIM_3D(hi),
				                  BL_TO_FORTRAN_FAB(umac_mid_mf[mfi]),
				                  BL_TO_FORTRAN_FAB(utemp_mf[mfi]));

			}
		}

        for (int lev=0; lev<=finest_level; ++lev) {
			MultiFab& gphi_mac_mf = gphi_mac[lev];
			MultiFab& utemp_mf = utemp[lev];
			const Real* dx = geom[lev].CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
			for (MFIter mfi(gphi_mac_mf, true); mfi.isValid(); ++mfi)
			{
				const Box& tilebox = mfi.tilebox();
				const int* lo  = tilebox.loVect();
				const int* hi  = tilebox.hiVect();

				convert_MAC_to_cc(ARLIM_3D(lo), ARLIM_3D(hi),
				                  BL_TO_FORTRAN_FAB(gphi_mac_mf[mfi]),
				                  BL_TO_FORTRAN_FAB(utemp_mf[mfi]));

			}
		}

		for (int lev=0; lev<=finest_level; ++lev) {
			for (auto comp=0; comp < AMREX_SPACEDIM; ++comp)
				MultiFab::Copy(umac_new[lev], umac_mid[lev], comp, comp, 1,
				               umac_new[lev].nGrow());
		}
	}

	//--------------------------------------------------------------------------
	// project out the divergent portion of the velocity field
	//--------------------------------------------------------------------------

	if (project_type == 1) {
		// hgprojection -- here pi is nodal and u is cell-centered

		for (int lev=0; lev<=finest_level; ++lev) {
			sold[lev].setVal(0.);
			snew[lev].setVal(0.);
			sold[lev].setVal(1., Rho, 1, 1);
			snew[lev].setVal(1., Rho, 1, 1);
			uold[lev].setVal(0.);
			unew[lev].setVal(0.);
			pi[lev].setVal(0.);
			gpi[lev].setVal(0.);
			rhcc_for_nodalproj[lev].setVal(0.);
		}

		std::fill(beta0_old.begin(), beta0_old.end(), 1.);
		std::fill(beta0_nph.begin(), beta0_nph.end(), 1.);
		std::fill(beta0_new.begin(), beta0_new.end(), 1.);

		t_new = t_old + 1.;

		// hgproject
		NodalProj(initial_projection_comp, rhcc_for_nodalproj);

	} else {
		// mac projection -- here pi is cell-centered and u is MAC

		Vector<MultiFab> rhohalf(finest_level+1);
		Vector<MultiFab> macpi(finest_level+1);
		Vector<MultiFab> macrhs(finest_level+1);

		std::fill(beta0_old.begin(), beta0_old.end(), 1.);

		for (int lev=0; lev<=finest_level; ++lev) {
			// cell-centered MultiFabs
			rhohalf[lev].define(grids[lev], dmap[lev], Nscal,    1);
			rhohalf[lev].setVal(1.);
			macpi[lev].define(grids[lev], dmap[lev],       1,    1);
			macpi[lev].setVal(0.);
			macrhs[lev].define(grids[lev], dmap[lev],       1,    1);
			macrhs[lev].setVal(0.);

			sold[lev].setVal(0.);
			snew[lev].setVal(0.);
			sold[lev].setVal(1., Rho, 1, 1);
			snew[lev].setVal(1., Rho, 1, 1);
		}

		// macproject
		auto is_predictor = 0;
		MacProj(umac_new,macpi,macrhs,beta0_old,is_predictor);

        for (int lev=0; lev<=finest_level; ++lev) {
			MultiFab& umac_new_mf = umac_new[lev];
			MultiFab& utemp_mf = utemp[lev];
			const Real* dx = geom[lev].CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
			for (MFIter mfi(umac_new_mf, true); mfi.isValid(); ++mfi)
			{
				const Box& tilebox = mfi.tilebox();
				const int* lo  = tilebox.loVect();
				const int* hi  = tilebox.hiVect();

				convert_MAC_to_cc(ARLIM_3D(lo), ARLIM_3D(hi),
				                  BL_TO_FORTRAN_FAB(umac_new_mf[mfi]),
				                  BL_TO_FORTRAN_FAB(utemp_mf[mfi]));

			}
		}

	}

}

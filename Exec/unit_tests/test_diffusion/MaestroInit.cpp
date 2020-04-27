
#include <Maestro.H>
#include <Maestro_F.H>
#include <AMReX_VisMF.H>
using namespace amrex;


void
Maestro::Init ()
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::Init()",Init);

	Print() << "Calling Init()" << std::endl;

	start_step = 1;

	// fill in multifab and base state data
	InitData();

	// set finest_radial_level in fortran
	// compute numdisjointchunks, r_start_coord, r_end_coord
	init_multilevel(tag_array.dataPtr(),&finest_level);

	// compute initial time step
	FirstDt();

	if (stop_time >= 0. && t_old+dt > stop_time) {
		dt = std::min(dt,stop_time-t_old);
		Print() << "Stop time limits dt = " << dt << std::endl;
	}

	dtold = dt;
	t_new = t_old + dt;

}

// fill in multifab and base state data
void
Maestro::InitData ()
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::InitData()",InitData);

	Print() << "Calling InitData()" << std::endl;

	// read in model file and fill in s0_init and p0_init for all levels
	for (auto lev = 0; lev <= base_geom.max_radial_level; ++lev) {
	        InitBaseState(rho0_old, rhoh0_old,
			      p0_old, lev);
	}
	
	// calls AmrCore::InitFromScratch(), which calls a MakeNewGrids() function
	// that repeatedly calls Maestro::MakeNewLevelFromScratch() to build and initialize
	InitFromScratch(t_old);

	// reset tagging array to include buffer zones
	TagArray();

	// set finest_radial_level in fortran
	// compute numdisjointchunks, r_start_coord, r_end_coord
	init_multilevel(tag_array.dataPtr(),&finest_level);

	// average down data and fill ghost cells
	AverageDown(sold,0,Nscal);
	FillPatch(t_old,sold,sold,sold,0,0,Nscal,0,bcs_s);

	// free memory in s0_init and p0_init by swapping it
	// with an empty vector that will go out of scope
	Vector<Real> s0_swap, p0_swap;
	std::swap(s0_swap,s0_init);
	std::swap(p0_swap,p0_init);

    p0_new = p0_old;
}

// During initialization of a simulation, Maestro::InitData() calls
// AmrCore::InitFromScratch(), which calls
// a MakeNewGrids() function that repeatedly calls this function to build
// and initialize finer levels.  This function creates a new fine
// level that did not exist before by interpolating from the coarser level
// overrides the pure virtual function in AmrCore
void Maestro::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
                                       const DistributionMapping& dm)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::MakeNewLevelFromScratch()",MakeNewLevelFromScratch);

	sold              [lev].define(ba, dm,          Nscal, ng_s);
	snew              [lev].define(ba, dm,          Nscal, ng_s);
	S_cc_old          [lev].define(ba, dm,              1,    0);

	sold              [lev].setVal(0.);
	snew              [lev].setVal(0.);
	S_cc_old          [lev].setVal(0.);

	const auto prob_lo = geom[lev].ProbLoArray();
	const auto dx = geom[lev].CellSizeArray();

	const auto center_p = center;
	const Real * AMREX_RESTRICT s0_p = s0_init.dataPtr();
	auto& p0_p = p0_old;

	const int max_lev = base_geom.max_radial_level + 1;
	const auto nrf = base_geom.nr_fine;

	const auto peak_h_loc = peak_h;
	const auto ambient_h_loc = ambient_h;
	const auto diff_coeff = diffusion_coefficient;
	const auto t0_loc = t0;

	// convergence parameters
	const auto max_iter = 50;
	const auto tol = 1.e-12;

	// Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(sold[lev], true); mfi.isValid(); ++mfi)
	{
		const Box& tileBox = mfi.tilebox();

		const Array4<Real> scal = sold[lev].array(mfi);

		// set scalars to zero 
		AMREX_PARALLEL_FOR_4D(tileBox, Nscal, i, j, k, n, {
			scal(i,j,k,n) = 0.0;
		});

		AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
			const auto r = j;

			const auto x = prob_lo[0] + (Real(i) + 0.5) * dx[0] - center_p[0];
			const auto y = prob_lo[1] + (Real(j) + 0.5) * dx[1] - center_p[1];

			// apply the guassian enthalpy pulse at constant density
			const auto dist2 = x*x + y*y;

			const auto h_zone = (peak_h_loc - ambient_h_loc) * std::exp(-dist2 / (4.0 * diff_coeff * t0_loc)) + ambient_h_loc;

			auto temp_zone = s0_p[lev+max_lev*(r+nrf*Temp)];

			eos_t eos_state;

			for (auto comp = 0; comp < NumSpec; ++comp) {
				eos_state.xn[comp] = s0_p[lev+max_lev*(r+nrf*(FirstSpec+comp))] / s0_p[lev+max_lev*(r+nrf*Rho)];
			}

			eos_state.rho = s0_p[lev+max_lev*(r+nrf*Rho)];

			auto converged = false;

			for (auto iter = 0; iter < max_iter; ++iter) {
				eos_state.T = temp_zone;

				eos(eos_input_rt, eos_state);

				const auto dhdt = eos_state.cv + eos_state.dpdT / eos_state.rho;

				const auto del_temp = -(eos_state.h - h_zone) / dhdt;

				temp_zone += del_temp;

				if (fabs(del_temp) < tol*temp_zone) {
					converged = true;
					break;
				}
			}

			if (!converged) {
				Abort("Iters did not converge in InitLevelData.");
			}

			// call eos one last time
			eos_state.T = temp_zone;

			eos(eos_input_rt, eos_state);
			
			scal(i,j,k,Rho) = eos_state.rho;
			scal(i,j,k,RhoH) = eos_state.rho * eos_state.h;
			scal(i,j,k,Temp) = eos_state.T;
			for (auto comp = 0; comp < NumSpec; ++comp) {
				scal(i,j,k,FirstSpec+comp) = eos_state.xn[comp] * eos_state.rho;
			}

			p0_p(lev,r) = eos_state.p;
		}); 
	}
}

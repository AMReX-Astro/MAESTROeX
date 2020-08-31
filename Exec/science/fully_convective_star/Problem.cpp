#include "Maestro.H"

using namespace amrex;

#ifdef DO_PROBLEM_POST_INIT

void Maestro::ProblemPostInit() {
    BL_PROFILE_VAR("Maestro::ProblemPostInit()", ProblemPostInit);
    
    // initialize radial diagnostic
    // 0: convective velocity
    // 1: meridional circulation (radial)
    postfile_data.resize(base_geom.max_radial_level+1, base_geom.nr_fine, 2);
    postfile_data.setVal(0.);

    // compute diagnostics for initial state
    // make radial and circumferential velocity
    Vector<MultiFab> rad_vel(finest_level + 1);
    Vector<MultiFab> circ_vel(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
	rad_vel[lev] .define(grids[lev], dmap[lev], 1, 0);
	circ_vel[lev].define(grids[lev], dmap[lev], 1, 0);
    }
    MakeVelrc(uold, w0_cart, rad_vel, circ_vel);

    // make convection velocity
    BaseState<Real> convect_vel(base_geom.max_radial_level + 1,
				base_geom.nr_fine);
    MakeConvectionVel(rad_vel, convect_vel, postfile_data, 0, t_new);

    // wrote a plotfile
    if (plot_int > 0 || plot_deltat > 0) {
	std::string plotfilename = plot_base_name + "0000000";

	VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
    
	// write out the cell-centered diagnostics
	if (ParallelDescriptor::IOProcessor()) {
	    for (int lev = 0; lev <= base_geom.max_radial_level; ++lev) {
		std::ofstream PostFile;
		PostFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(),
					    io_buffer.size());
		std::string PostFileName(plotfilename + "/PostTimestep_");
		std::string levStr = std::to_string(lev);
		PostFileName.append(levStr);
		PostFile.open(PostFileName.c_str(), std::ofstream::out |
                                                    std::ofstream::trunc |
                                                    std::ofstream::binary);
		if (!PostFile.good()) {
		    amrex::FileOpenFailed(PostFileName);
		}

		PostFile.precision(17);

		PostFile << "r_cc  convection_vel   \n";

		for (int i = 0; i < base_geom.nr(lev); ++i) {
		    PostFile << base_geom.r_cc_loc(lev, i) << " "
			     << convect_vel.array()(lev, i) << "\n";
		}
	    }
        }
    }
}

#endif

#ifdef DO_PROBLEM_POST_TIMESTEP

void Maestro::ProblemPostTimestep() {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::ProblemPostInit()", ProblemPostInit);

    // compute diagnostics after every time step
    // make radial and circumferential velocity
    Vector<MultiFab> rad_vel(finest_level + 1);
    Vector<MultiFab> circ_vel(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
	rad_vel[lev] .define(grids[lev], dmap[lev], 1, 0);
	circ_vel[lev].define(grids[lev], dmap[lev], 1, 0);
    }
    MakeVelrc(unew, w0_cart, rad_vel, circ_vel);

    // make convection velocity
    BaseState<Real> convect_vel(base_geom.max_radial_level + 1,
				base_geom.nr_fine);
    MakeConvectionVel(rad_vel, convect_vel, postfile_data, 0, t_new);
    
    // wrote a plotfile       
    if ((plot_int > 0 && istep % plot_int == 0) ||
        (plot_deltat > 0 && std::fmod(t_new, plot_deltat) < dt) ||
	((plot_int > 0 || plot_deltat > 0) &&
	 (istep == max_step || t_old >= stop_time))) {

	std::string plotfilename = plot_base_name;
        PlotFileName(istep, &plotfilename);
	
	VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
    
	// write out the cell-centered diagnostics
	if (ParallelDescriptor::IOProcessor()) {
	    for (int lev = 0; lev <= base_geom.max_radial_level; ++lev) {
		std::ofstream PostFile;
		PostFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(),
					    io_buffer.size());
		std::string PostFileName(plotfilename + "/PostTimestep_");
		std::string levStr = std::to_string(lev);
		PostFileName.append(levStr);
		PostFile.open(PostFileName.c_str(), std::ofstream::out |
                                                    std::ofstream::trunc |
                                                    std::ofstream::binary);
		if (!PostFile.good()) {
		    amrex::FileOpenFailed(PostFileName);
		}

		PostFile.precision(17);
		
		PostFile << "r_cc  convection_vel   \n";

		for (int i = 0; i < base_geom.nr(lev); ++i) {
		    PostFile << base_geom.r_cc_loc(lev, i) << " "
			     << convect_vel.array()(lev, i) << "\n";
		}
	    }
	}
    }
}

#endif

void Maestro::MakeConvectionVel(const Vector<MultiFab>& velr,
				BaseState<Real>& vel_conv,
				BaseState<Real>& totsum, const int comp,
				const Real t_interval) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeConvectionVel()", MakeConvectionVel);
    
    Vector<MultiFab> vel2(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        vel2[lev].define(velr[lev].boxArray(), velr[lev].DistributionMap(), 1,
                         0);
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(velr[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const Array4<const Real> velr_arr = velr[lev].array(mfi);
            const Array4<Real> vel2_arr = vel2[lev].array(mfi);

            AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                vel2_arr(i, j, k) = 0.0;

                // square radial vel
                vel2_arr(i, j, k) = velr_arr(i, j, k) * velr_arr(i, j, k);
            });
        }
    }

    // average down and fill ghost cells
    AverageDown(vel2, 0, 1);
    // FillPatch(t0, vel2, vel2, vel2, 0, 0, 1, 0, bcs_f);

    // radial average of square of radial vel
    Average(vel2, vel_conv, 0);

    // add vel^2 to total sum and compute convection vel
    auto vel_conv_arr = vel_conv.array();
    auto totsum_arr = totsum.array();

    for (auto r = 0; r < base_geom.nr_fine; ++r) {
	totsum_arr(0, r, comp) += dt*vel_conv_arr(0, r);
	
	// root-mean-squared radial velocity
        vel_conv_arr(0, r) = std::sqrt(totsum_arr(0, r) / t_interval);
    }
}

#include "Maestro.H"

using namespace amrex;

#ifdef DO_PROBLEM_POST_INIT

void Maestro::ProblemPostInit() {
    BL_PROFILE_VAR("Maestro::ProblemPostInit()", ProblemPostInit);
    
    // initialize radial diagnostic
    // 0: convective velocity
    // 1: latitudinal shear
    postfile_data.resize(base_geom.max_radial_level+1, base_geom.nr_fine, 2);

    // compute diagnostics for initial state


    // wrote a plotfile
    std::string plotfilename = plot_base_name;

    if (plotfilename.back() == '_') {
	plotfilename += "InitData";
    } else {
	plotfilename += "_InitData";
    }

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

            PostFile << "r_cc  rho0   \n";

            for (int i = 0; i < base_geom.nr(lev); ++i) {
                PostFile << base_geom.r_cc_loc(lev, i) << " "
			 << rho0_old.array()(lev, i) << "\n";
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
    
    
    // wrote a plotfile       
    if ((plot_int > 0 && istep % plot_int == 0) ||
        (plot_deltat > 0 && std::fmod(t_new, plot_deltat) < dt) ||
        (istep == max_step) || (t_old >= stop_time)) {

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
		
		PostFile << "r_cc  rho0   \n";

		for (int i = 0; i < base_geom.nr(lev); ++i) {
		    PostFile << base_geom.r_cc_loc(lev, i) << " "
			     << rho0_new.array()(lev, i) << "\n";
		}
	    }
	}
    }
}

#endif

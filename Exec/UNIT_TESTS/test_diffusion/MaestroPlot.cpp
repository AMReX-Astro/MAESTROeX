
#include <Maestro.H>
#include <MaestroPlot.H>
#include <AMReX_buildInfo.H>

using namespace amrex;

// write plotfile to disk
void
Maestro::WritePlotFile (const int step,
                        const Real t_in,
                        const Real dt_in,
                        const Vector<Real>& a,
                        const Vector<Real>& b,
                        const Vector<Real>& c,
                        const Vector<Real>& d,
                        const Vector<MultiFab>& state,
                        Vector<MultiFab>& analytic,
                        const Vector<MultiFab>& err)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::WritePlotFile()",WritePlotFile);

	// wallclock time
	const Real strt_total = ParallelDescriptor::second();

	std::string plotfilename = PlotFileName(step);

	int nPlot = 0;
	const auto& varnames = PlotFileVarNames(&nPlot);
	const Vector<std::string> analytic_varnames = {"h"};
	const Vector<std::string> error_varnames = {"error"};

    // make plot mfs
    const Vector<MultiFab> dummy;
	const auto& state_mf = PlotFileMF(nPlot,0,dt_in,state,dummy,dummy,
        dummy,dummy,analytic,a,a,dummy);
	const auto& analytic_mf = PlotFileMF(nPlot,1,dt_in,dummy,dummy,
        dummy,dummy,dummy,analytic,a,a,dummy);
	const auto& error_mf = PlotFileMF(nPlot,2,dt_in,err,dummy,dummy,
        dummy,dummy,analytic,a,a,dummy);

	// WriteMultiLevelPlotfile expects an array of step numbers
	Vector<int> step_array;
	step_array.resize(maxLevel()+1, step);

	WriteMultiLevelPlotfile(plotfilename, finest_level+1, state_mf, varnames,
	                        Geom(), t_in, step_array, refRatio());
	WriteMultiLevelPlotfile("analytic_" + plotfilename, finest_level+1, analytic_mf,
	                        analytic_varnames,
	                        Geom(), t_in, step_array, refRatio());
	WriteMultiLevelPlotfile("error_" + plotfilename, finest_level+1, error_mf,
	                        error_varnames,
	                        Geom(), t_in, step_array, refRatio());

}


// get plotfile name
std::string
Maestro::PlotFileName (int lev) const
{
	return Concatenate(plot_base_name, lev, 7);
}

// put together a vector of multifabs for writing
Vector<const MultiFab*>
Maestro::PlotFileMF (const int nPlot,
                     const Real t_in,
                     const Real dt_in,
                     const Vector<MultiFab>& s_in,
                     const Vector<MultiFab>& h,
                     const Vector<MultiFab>& j,
                     const Vector<MultiFab>& a,
                     const Vector<MultiFab>& b,
                     Vector<MultiFab>& analytic,
                     const Vector<Real>& d,
                     const Vector<Real>& e,
                     const Vector<MultiFab>& f)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::PlotFileMF()",PlotFileMF);

	// MultiFab to hold plotfile data
	Vector<const MultiFab*> plot_mf;

	// temporary MultiFab to hold plotfile data
	Vector<MultiFab*> plot_mf_data(finest_level+1);

	int dest_comp = 0;

    if (t_in == 0) { // this is the multifab containing the full state

    	// build temporary MultiFab to hold plotfile data
    	for (int i = 0; i <= finest_level; ++i) {
    		plot_mf_data[i] = new MultiFab((s_in[i]).boxArray(),(s_in[i]).DistributionMap(),nPlot,0);
    	}

    	// rho
    	for (int i = 0; i <= finest_level; ++i) {
    		plot_mf_data[i]->copy((s_in[i]),Rho,dest_comp,1);
    	}
    	++dest_comp;

        // h
    	for (int i = 0; i <= finest_level; ++i) {
    		plot_mf_data[i]->copy((s_in[i]),RhoH,dest_comp,1);
    		MultiFab::Divide(*plot_mf_data[i], s_in[i], Rho, dest_comp, 1, 0);
    	}
    	++dest_comp;

    	// rhoX
    	for (int i = 0; i <= finest_level; ++i) {
    		plot_mf_data[i]->copy((s_in[i]),FirstSpec,dest_comp,NumSpec);
    	}
    	dest_comp += NumSpec;

    	// compute tfromh
    	for (int i = 0; i <= finest_level; ++i) {
    		// tfromh
    		plot_mf_data[i]->copy((s_in[i]),Temp,dest_comp,1);
    	}
    	++dest_comp;

    } else if (t_in == 1) {
        // this is the multifab containing the analytic solution

    	// build temporary MultiFab to hold plotfile data
    	for (int i = 0; i <= finest_level; ++i) {
    		plot_mf_data[i] = new MultiFab((analytic[i]).boxArray(),(analytic[i]).DistributionMap(),1,0);
    	}

    	for (int i = 0; i <= finest_level; ++i) {
    		plot_mf_data[i]->copy((analytic[i]),0,dest_comp,1);
    	}
    	++dest_comp;
    } else {
        // this is the multifab containing the error

    	// build temporary MultiFab to hold plotfile data
    	for (int i = 0; i <= finest_level; ++i) {
    		plot_mf_data[i] = new MultiFab((s_in[i]).boxArray(),(s_in[i]).DistributionMap(),1,0);
    	}

        for (int i = 0; i <= finest_level; ++i) {
    		plot_mf_data[i]->copy((s_in[i]),0,dest_comp,1);
    	}
    	++dest_comp;

    }

	// add plot_mf_data[i] to plot_mf
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf.push_back(plot_mf_data[i]);
		// delete [] plot_mf_data[i];
	}

	return plot_mf;

}

// set plotfile variable names
Vector<std::string>
Maestro::PlotFileVarNames (int * nPlot) const
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::PlotFileVarNames()",PlotFileVarNames);

	(*nPlot) = 3 + NumSpec;

	Vector<std::string> names(*nPlot);

	int cnt = 0;

	// density and enthalpy
	names[cnt++] = "rho";
	names[cnt++] = "h";

	for (int i = 0; i < NumSpec; i++) {
		int len = 20;
		Vector<int> int_spec_names(len);
		//
		// This call return the actual length of each string in "len"
		//
		get_spec_names(int_spec_names.dataPtr(),&i,&len);
		char* spec_name = new char[len+1];
		for (int j = 0; j < len; j++)
			spec_name[j] = int_spec_names[j];
		spec_name[len] = '\0';
		std::string spec_string = "rhoX(";
		spec_string += spec_name;
		spec_string += ')';

		names[cnt++] = spec_string;
		delete [] spec_name;
	}

	names[cnt++] = "tfromp";

	return names;

}


void
Maestro::WriteJobInfo (const std::string& dir) const
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::WriteJobInfo()",WriteJobInfo);

	if (ParallelDescriptor::IOProcessor())
	{
		// job_info file with details about the run
		std::ofstream jobInfoFile;
		std::string FullPathJobInfoFile = dir;

		std::string PrettyLine = std::string(78, '=') + "\n";
		std::string OtherLine = std::string(78, '-') + "\n";
		std::string SkipSpace = std::string(8, ' ') + "\n";

		FullPathJobInfoFile += "/job_info";
		jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);

		// job information
		jobInfoFile << PrettyLine;
		jobInfoFile << " MAESTROeX Job Information\n";
		jobInfoFile << PrettyLine;

		jobInfoFile << "job name: " << job_name << "\n\n";
		jobInfoFile << "inputs file: " << inputs_name << "\n\n";

		jobInfoFile << "number of MPI processes: " << ParallelDescriptor::NProcs() << "\n";
#ifdef _OPENMP
		jobInfoFile << "number of threads:       " << omp_get_max_threads() << "\n";

		jobInfoFile << "tile size: ";
		for (int d=0; d<AMREX_SPACEDIM; ++d) {
			jobInfoFile << FabArrayBase::mfiter_tile_size[d] << " ";
		}
		jobInfoFile << "\n";
#endif
		jobInfoFile << "\n";
		jobInfoFile << "CPU time used since start of simulation (CPU-hours): " <<
		        getCPUTime()/3600.0;

		jobInfoFile << "\n\n";

		// plotfile information
		jobInfoFile << PrettyLine;
		jobInfoFile << " Plotfile Information\n";
		jobInfoFile << PrettyLine;

		time_t now = time(0);

		// Convert now to tm struct for local timezone
		tm* localtm = localtime(&now);
		jobInfoFile << "output data / time: " << asctime(localtm);

		char currentDir[FILENAME_MAX];
		if (getcwd(currentDir, FILENAME_MAX)) {
			jobInfoFile << "output dir:         " << currentDir << "\n";
		}

		jobInfoFile << "\n\n";

		// build information
		jobInfoFile << PrettyLine;
		jobInfoFile << " Build Information\n";
		jobInfoFile << PrettyLine;

		jobInfoFile << "build date:    " << buildInfoGetBuildDate() << "\n";
		jobInfoFile << "build machine: " << buildInfoGetBuildMachine() << "\n";
		jobInfoFile << "build dir:     " << buildInfoGetBuildDir() << "\n";
		jobInfoFile << "AMReX dir:     " << buildInfoGetAMReXDir() << "\n";

		jobInfoFile << "\n";

		jobInfoFile << "COMP:          " << buildInfoGetComp() << "\n";
		jobInfoFile << "COMP version:  " << buildInfoGetCompVersion() << "\n";

		jobInfoFile << "\n";

		jobInfoFile << "C++ compiler:  " << buildInfoGetCXXName() << "\n";
		jobInfoFile << "C++ flags:     " << buildInfoGetCXXFlags() << "\n";

		jobInfoFile << "\n";

		jobInfoFile << "Fortran comp:  " << buildInfoGetFName() << "\n";
		jobInfoFile << "Fortran flags: " << buildInfoGetFFlags() << "\n";

		jobInfoFile << "\n";

		jobInfoFile << "Link flags:    " << buildInfoGetLinkFlags() << "\n";
		jobInfoFile << "Libraries:     " << buildInfoGetLibraries() << "\n";

		jobInfoFile << "\n";

		for (int n = 1; n <= buildInfoGetNumModules(); n++) {
			jobInfoFile << buildInfoGetModuleName(n) << ": " << buildInfoGetModuleVal(n) << "\n";
		}

		const char* githash1 = buildInfoGetGitHash(1);
		const char* githash2 = buildInfoGetGitHash(2);
		const char* githash3 = buildInfoGetGitHash(3);
		if (strlen(githash1) > 0) {
			jobInfoFile << "MAESTROeX git describe: " << githash1 << "\n";
		}
		if (strlen(githash2) > 0) {
			jobInfoFile << "AMReX git describe: " << githash2 << "\n";
		}
		if (strlen(githash3) > 0) {
			jobInfoFile << "Microphysics git describe: " << githash3 << "\n";
		}

		const char* buildgithash = buildInfoGetBuildGitHash();
		const char* buildgitname = buildInfoGetBuildGitName();
		if (strlen(buildgithash) > 0) {
			jobInfoFile << buildgitname << " git describe: " << buildgithash << "\n";
		}

		jobInfoFile << "\n\n";

		// grid information
		jobInfoFile << PrettyLine;
		jobInfoFile << " Grid Information\n";
		jobInfoFile << PrettyLine;

		for (int i = 0; i <= finest_level; i++)
		{
			jobInfoFile << " level: " << i << "\n";
			jobInfoFile << "   number of boxes = " << grids[i].size() << "\n";
			jobInfoFile << "   maximum zones   = ";
			for (int n = 0; n < BL_SPACEDIM; n++)
			{
				jobInfoFile << geom[i].Domain().length(n) << " ";
			}
			jobInfoFile << "\n\n";
		}

		jobInfoFile << " Boundary conditions\n";
		Vector<int> lo_bc_out(BL_SPACEDIM), hi_bc_out(BL_SPACEDIM);
		ParmParse pp("maestro");
		pp.getarr("lo_bc",lo_bc_out,0,BL_SPACEDIM);
		pp.getarr("hi_bc",hi_bc_out,0,BL_SPACEDIM);


// these names correspond to the integer flags setup in the
// Castro_setup.cpp
		const char* names_bc[] =
		{ "interior", "inflow", "outflow",
		  "symmetry", "slipwall", "noslipwall" };


		jobInfoFile << "   -x: " << names_bc[lo_bc_out[0]] << "\n";
		jobInfoFile << "   +x: " << names_bc[hi_bc_out[0]] << "\n";
		if (BL_SPACEDIM >= 2) {
			jobInfoFile << "   -y: " << names_bc[lo_bc_out[1]] << "\n";
			jobInfoFile << "   +y: " << names_bc[hi_bc_out[1]] << "\n";
		}
		if (BL_SPACEDIM == 3) {
			jobInfoFile << "   -z: " << names_bc[lo_bc_out[2]] << "\n";
			jobInfoFile << "   +z: " << names_bc[hi_bc_out[2]] << "\n";
		}

		jobInfoFile << "\n\n";

// species info
		Real Aion = 0.0;
		Real Zion = 0.0;

		int mlen = 20;

		jobInfoFile << PrettyLine;
		jobInfoFile << " Species Information\n";
		jobInfoFile << PrettyLine;

		jobInfoFile <<
		        std::setw(6) << "index" << SkipSpace <<
		        std::setw(mlen+1) << "name" << SkipSpace <<
		        std::setw(7) << "A" << SkipSpace <<
		        std::setw(7) << "Z" << "\n";
		jobInfoFile << OtherLine;

		for (int i = 0; i < NumSpec; i++)
		{

			int len = mlen;
			Vector<int> int_spec_names(len);
			//
			// This call return the actual length of each string in "len"
			//
			get_spec_names(int_spec_names.dataPtr(),&i,&len);
			char* spec_name = new char[len+1];
			for (int j = 0; j < len; j++)
				spec_name[j] = int_spec_names[j];
			spec_name[len] = '\0';

			// get A and Z
			get_spec_az(&i, &Aion, &Zion);

			jobInfoFile <<
			        std::setw(6) << i << SkipSpace <<
			        std::setw(mlen+1) << std::setfill(' ') << spec_name << SkipSpace <<
			        std::setw(7) << Aion << SkipSpace <<
			        std::setw(7) << Zion << "\n";
			delete [] spec_name;
		}
		jobInfoFile << "\n\n";

		// runtime parameters
		jobInfoFile << PrettyLine;
		jobInfoFile << " Inputs File Parameters\n";
		jobInfoFile << PrettyLine;

		ParmParse::dumpTable(jobInfoFile, true);

		jobInfoFile.close();

		// now the external parameters
		const int jobinfo_file_length = FullPathJobInfoFile.length();
		Vector<int> jobinfo_file_name(jobinfo_file_length);

		for (int i = 0; i < jobinfo_file_length; i++)
			jobinfo_file_name[i] = FullPathJobInfoFile[i];

	}
}

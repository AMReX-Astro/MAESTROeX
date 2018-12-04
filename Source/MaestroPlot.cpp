
#include <Maestro.H>
#include <MaestroPlot.H>
#include <AMReX_buildInfo.H>

using namespace amrex;

// write plotfile to disk
void
Maestro::WritePlotFile (const int step,
                        const Real t_in,
                        const Vector<Real>& rho0_in,
                        const Vector<Real>& rhoh0_in,
                        const Vector<Real>& p0_in,
                        const Vector<Real>& gamma1bar_in,
                        const Vector<MultiFab>& u_in,
                        Vector<MultiFab>& s_in)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::WritePlotFile()",WritePlotFile);

	// wallclock time
	const Real strt_total = ParallelDescriptor::second();

	std::string plotfilename;

	if (step == 9999999) {
		plotfilename = "plt_InitData";
	}
	else if (step == 9999998) {
		plotfilename = "plt_after_InitProj";
	}
	else if (step == 9999997) {
		plotfilename = "plt_after_DivuIter";
	}
	else {
		plotfilename = PlotFileName(step);
	}

	// convert rho0 to multi-D MultiFab
	Vector<MultiFab> rho0_cart(finest_level+1);
	for (int lev=0; lev<=finest_level; ++lev) {
		rho0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
	}
	Put1dArrayOnCart(rho0_in,rho0_cart,0,0);

	// convert rhoh0 to multi-D MultiFab
	Vector<MultiFab> rhoh0_cart(finest_level+1);
	for (int lev=0; lev<=finest_level; ++lev) {
		rhoh0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
	}
	Put1dArrayOnCart(rhoh0_in,rhoh0_cart,0,0);

	// convert p0 to multi-D MultiFab
	Vector<MultiFab> p0_cart(finest_level+1);
	for (int lev=0; lev<=finest_level; ++lev) {
		p0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
	}
	Put1dArrayOnCart(p0_in,p0_cart,0,0);

	// convert gamma1bar to multi-D MultiFab
	Vector<MultiFab> gamma1bar_cart(finest_level+1);
	for (int lev=0; lev<=finest_level; ++lev) {
		gamma1bar_cart[lev].define(grids[lev], dmap[lev], 1, 0);
	}
	Put1dArrayOnCart(gamma1bar_in,gamma1bar_cart,0,0);

	const auto& mf = PlotFileMF(rho0_cart,rhoh0_cart,p0_cart,gamma1bar_cart,u_in,s_in,p0_in,gamma1bar_in);
	const auto& varnames = PlotFileVarNames();

	// WriteMultiLevelPlotfile expects an array of step numbers
	Vector<int> step_array;
	step_array.resize(maxLevel()+1, step);

	WriteMultiLevelPlotfile(plotfilename, finest_level+1, mf, varnames,
	                        Geom(), t_in, step_array, refRatio());

	WriteJobInfo(plotfilename);

	// wallclock time
	Real end_total = ParallelDescriptor::second() - strt_total;

	// print wallclock time
	ParallelDescriptor::ReduceRealMax(end_total,ParallelDescriptor::IOProcessorNumber());
	if (maestro_verbose > 0) {
		Print() << "Time to write plotfile: " << end_total << '\n';
	}

	for (int i = 0; i <= finest_level; ++i) {
		delete mf[i];
	}

}


// get plotfile name
std::string
Maestro::PlotFileName (int lev) const
{
	return Concatenate(plot_base_name, lev, 7);
}

// put together a vector of multifabs for writing
Vector<const MultiFab*>
Maestro::PlotFileMF (const Vector<MultiFab>& rho0_cart,
                     const Vector<MultiFab>& rhoh0_cart,
                     const Vector<MultiFab>& p0_cart,
                     const Vector<MultiFab>& gamma1bar_cart,
                     const Vector<MultiFab>& u_in,
                     Vector<MultiFab>& s_in,
                     const Vector<Real>& p0_in,
                     const Vector<Real>& gamma1bar_in)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::PlotFileMF()",PlotFileMF);

	// velocities (AMREX_SPACEDIM)
	// magvel
	// rho, rhoh, rhoX, tfromp, tfromh, deltaT Pi (Nscal+2 -- the extra 2 are tfromh and deltaT)
	// X (NumSpec), omegadot(NumSpec)
	// rho' and rhoh' (2)
	// rho0, rhoh0, p0, w0 (3+AMREX_SPACEDIM)
	// MachNumber, deltagamma
	int nPlot = 2*AMREX_SPACEDIM + Nscal + 2*NumSpec + 11;

	// MultiFab to hold plotfile data
	Vector<const MultiFab*> plot_mf;

	// temporary MultiFab to hold plotfile data
	Vector<MultiFab*> plot_mf_data(finest_level+1);

	// temporary MultiFab for calculations
	Vector<MultiFab> tempmf(finest_level+1);
	Vector<MultiFab> tempmf_state(finest_level+1);
    Vector<MultiFab> tempmf_scalar1(finest_level+1);
    Vector<MultiFab> tempmf_scalar2(finest_level+1);


	int dest_comp = 0;

	// build temporary MultiFab to hold plotfile data
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i] = new MultiFab((s_in[i]).boxArray(),(s_in[i]).DistributionMap(),nPlot,0);
		tempmf[i].define(grids[i],dmap[i],AMREX_SPACEDIM,0);
        tempmf_state[i].define(grids[i],dmap[i],Nscal,0);
        tempmf_scalar1[i].define(grids[i],dmap[i],1,0);
        tempmf_scalar2[i].define(grids[i],dmap[i],1,0);
	}

	// velocity
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy((u_in[i]),0,dest_comp,AMREX_SPACEDIM);
	}
	dest_comp += AMREX_SPACEDIM;

	// magvel
	MakeMagvel(u_in, tempmf);
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy((tempmf[i]),0,dest_comp,1);
	}
	++dest_comp;

	// vorticity
	MakeVorticity(u_in, tempmf);
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy((tempmf[i]),0,dest_comp,1);
	}
	++dest_comp;

	// rho
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy((s_in[i]),Rho,dest_comp,1);
	}
	++dest_comp;

	// rhoh
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy((s_in[i]),RhoH,dest_comp,1);
	}
	++dest_comp;

	// rhoX
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy((s_in[i]),FirstSpec,dest_comp,NumSpec);
	}
	dest_comp += NumSpec;

	// X
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy((s_in[i]),FirstSpec,dest_comp,NumSpec);
		for (int comp=0; comp<NumSpec; ++comp) {
			MultiFab::Divide(*plot_mf_data[i],s_in[i],Rho,dest_comp+comp,1,0);
		}
	}
	dest_comp += NumSpec;

    // omegadot
    React(s_in, tempmf_state, tempmf_scalar1, tempmf, tempmf_scalar2, p0_in, dt);

    for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy((tempmf[i]),0,dest_comp,NumSpec);
		for (int comp=0; comp<NumSpec; ++comp) {
			MultiFab::Divide(*plot_mf_data[i],s_in[i],Rho,dest_comp+comp,1,0);
		}
	}
	dest_comp += NumSpec;

	// compute tfromp
	TfromRhoP(s_in,p0_in);

	// tfromp
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy((s_in[i]),Temp,dest_comp,1);
	}
	++dest_comp;

	// compute tfromh
	TfromRhoH(s_in,p0_in);

	for (int i = 0; i <= finest_level; ++i) {
		// tfromh
		plot_mf_data[i]->copy((s_in[i]),Temp,dest_comp,1);
	}
	++dest_comp;

	// deltaT
	// compute & copy tfromp
	TfromRhoP(s_in,p0_in);
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy((s_in[i]),Temp,dest_comp,1);
	}

	// compute tfromh
	TfromRhoH(s_in,p0_in);
	// compute deltaT = (tfromp - tfromh) / tfromh
	for (int i = 0; i <= finest_level; ++i) {
		MultiFab::Subtract(*plot_mf_data[i],s_in[i],Temp,dest_comp,1,0);
		MultiFab::Divide(*plot_mf_data[i],s_in[i],Temp,dest_comp,1,0);
	}
	++dest_comp;

	// restore tfromp if necessary
	if (use_tfromp) {
		TfromRhoP(s_in,p0_in);
	}

	// pi
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy((s_in[i]),Pi,dest_comp,1);
	}
	++dest_comp;

	// rhopert
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy((s_in[i]),Rho,dest_comp,1);
		MultiFab::Subtract(*plot_mf_data[i],rho0_cart[i],0,dest_comp,1,0);
	}
	++dest_comp;

	// rhohpert
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy((s_in[i]),RhoH,dest_comp,1);
		MultiFab::Subtract(*plot_mf_data[i],rhoh0_cart[i],0,dest_comp,1,0);
	}
	++dest_comp;

	// rho0, rhoh0, and p0
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy(( rho0_cart[i]),0,dest_comp,1);
		plot_mf_data[i]->copy((rhoh0_cart[i]),0,dest_comp+1,1);
		plot_mf_data[i]->copy((   p0_cart[i]),0,dest_comp+2,1);
	}
	dest_comp += 3;

	if (spherical == 1) {
		Vector<std::array< MultiFab, AMREX_SPACEDIM > > w0mac(finest_level+1);
		Vector<MultiFab> w0r_cart(finest_level+1);

		for (int lev=0; lev<=finest_level; ++lev) {
			// w0mac will contain an edge-centered w0 on a Cartesian grid,
			// for use in computing divergences.
			AMREX_D_TERM(w0mac[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1); ,
			             w0mac[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1); ,
			             w0mac[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1); );
			for (int idim=0; idim<AMREX_SPACEDIM; ++idim) {
				w0mac[lev][idim].setVal(0.);
			}

			// w0r_cart is w0 but onto a Cartesian grid in cell-centered as
			// a scalar.  Since w0 is the radial expansion velocity, w0r_cart
			// is the radial w0 in a zone
			w0r_cart[lev].define(grids[lev], dmap[lev], 1, 0);
			w0r_cart[lev].setVal(0.);
		}

		if (evolve_base_state == 1 && use_exact_base_state == 0) {
			MakeW0mac(w0mac);
			Put1dArrayOnCart(w0,w0r_cart,1,0,bcs_f,0);
		}

		// Mach number
		MachfromRhoHSphr(s_in,u_in,p0_in,w0r_cart,tempmf);
	} else {
		// Mach number
		MachfromRhoH(s_in,u_in,p0_in,tempmf);
	}

	// MachNumber
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy((tempmf[i]),0,dest_comp,1);
	}
	++dest_comp;

	// deltagamma
	MakeDeltaGamma(s_in, p0_in, p0_cart, gamma1bar_in, gamma1bar_cart, tempmf);
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy((tempmf[i]),0,dest_comp,1);
	}
	++dest_comp;

	// w0
	Put1dArrayOnCart(w0,tempmf,1,1,bcs_u,0);
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy((tempmf[i]),0,dest_comp,AMREX_SPACEDIM);
	}
	dest_comp += AMREX_SPACEDIM;

	// add plot_mf_data[i] to plot_mf
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf.push_back(plot_mf_data[i]);
		// delete [] plot_mf_data[i];
	}

	return plot_mf;


}

// set plotfile variable names
Vector<std::string>
Maestro::PlotFileVarNames () const
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::PlotFileVarNames()",PlotFileVarNames);

	// velocities (AMREX_SPACEDIM)
	// magvel
	// rho, rhoh, rhoX, tfromp, tfromh, deltaT Pi (Nscal+2 -- the extra 2 are tfromh and deltaT)
	// X (NumSpec), omegadot(NumSpec)
	// rho' and rhoh' (2)
	// rho0, rhoh0, p0, w0 (3+AMREX_SPACEDIM)
	// MachNumber, deltagamma
	int nPlot = 2*AMREX_SPACEDIM + Nscal + 2*NumSpec + 11;
	Vector<std::string> names(nPlot);

	int cnt = 0;

	// add velocities
	for (int i=0; i<AMREX_SPACEDIM; ++i) {
		std::string x = "vel";
		x += (120+i);
		names[cnt++] = x;
	}

	names[cnt++] = "magvel";
	names[cnt++] = "vort";

	// density and enthalpy
	names[cnt++] = "rho";
	names[cnt++] = "rhoh";

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

	for (int i = 0; i < NumSpec; i++) {
		int len = 20;
		Vector<int> int_spec_names(len);
		//
		// This call return the actual length of each string in "len"
		//
		get_spec_names(int_spec_names.dataPtr(),&i,&len);
		char* spec_name = new char[len+1];
		for (int j = 0; j < len; j++) {
			spec_name[j] = int_spec_names[j];
		}
		spec_name[len] = '\0';
		std::string spec_string = "X(";
		spec_string += spec_name;
		spec_string += ')';

		names[cnt++] = spec_string;

		delete [] spec_name;
	}

    for (int i = 0; i < NumSpec; i++) {
		int len = 20;
		Vector<int> int_spec_names(len);
		//
		// This call return the actual length of each string in "len"
		//
		get_spec_names(int_spec_names.dataPtr(),&i,&len);
		char* spec_name = new char[len+1];
		for (int j = 0; j < len; j++) {
			spec_name[j] = int_spec_names[j];
		}
		spec_name[len] = '\0';
		std::string spec_string = "omegadot(";
		spec_string += spec_name;
		spec_string += ')';

		names[cnt++] = spec_string;

		delete [] spec_name;
	}


	names[cnt++] = "tfromp";
	names[cnt++] = "tfromh";
	names[cnt++] = "deltaT";
	names[cnt++] = "Pi";

	names[cnt++] = "rhopert";
	names[cnt++] = "rhohpert";

	names[cnt++] = "rho0";
	names[cnt++] = "rhoh0";
	names[cnt++] = "p0";
	names[cnt++] = "MachNumber";
	names[cnt++] = "deltagamma";

	// add w0
	for (int i=0; i<AMREX_SPACEDIM; ++i) {
		std::string x = "w0";
		x += (120+i);
		names[cnt++] = x;
	}

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

		// runtime_pretty_print(jobinfo_file_name.dataPtr(), &jobinfo_file_length);
	}
}

void
Maestro::MakeMagvel (const Vector<MultiFab>& vel,
                     Vector<MultiFab>& magvel)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::MakeMagvel()",MakeMagvel);

	Vector<MultiFab> w0_cart(finest_level+1);
	if (spherical == 1) {
		for (int lev=0; lev<=finest_level; ++lev) {
			w0_cart[lev].define(grids[lev], dmap[lev], 1, 1);
		}
	}

	for (int lev=0; lev<=finest_level; ++lev) {

		// get references to the MultiFabs at level lev
		const MultiFab& vel_mf = vel[lev];
		MultiFab& magvel_mf = magvel[lev];

		// Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
		if (spherical == 0) {
#ifdef _OPENMP
#pragma omp parallel
#endif
			for ( MFIter mfi(vel_mf, true); mfi.isValid(); ++mfi ) {

				// Get the index space of the valid region
				const Box& tileBox = mfi.tilebox();

				// call fortran subroutine
				// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
				// lo/hi coordinates (including ghost cells), and/or the # of components
				// We will also pass "validBox", which specifies the "valid" region.
				make_magvel(&lev,ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
				            BL_TO_FORTRAN_3D(vel_mf[mfi]),
				            w0.dataPtr(),
				            BL_TO_FORTRAN_3D(magvel_mf[mfi]));
			}

		} else {

			const MultiFab& w0cart_mf = w0_cart[lev];

#ifdef _OPENMP
#pragma omp parallel
#endif
			for ( MFIter mfi(vel_mf, true); mfi.isValid(); ++mfi ) {

				// Get the index space of the valid region
				const Box& tileBox = mfi.tilebox();

				// call fortran subroutine
				// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
				// lo/hi coordinates (including ghost cells), and/or the # of components
				// We will also pass "validBox", which specifies the "valid" region.

				make_magvel_sphr(ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
				                 BL_TO_FORTRAN_3D(vel_mf[mfi]),
				                 BL_TO_FORTRAN_3D(w0cart_mf[mfi]),
				                 BL_TO_FORTRAN_3D(magvel_mf[mfi]));
			}
		}
	}

	// average down and fill ghost cells
	AverageDown(magvel,0,1);
	FillPatch(t_old,magvel,magvel,magvel,0,0,1,0,bcs_f);
}


void
Maestro::MakeVelrc (const Vector<MultiFab>& vel,
                    Vector<MultiFab>& rad_vel,
                    Vector<MultiFab>& circ_vel)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::MakeVelrc()",MakeVelrc);

	Vector<MultiFab> w0r_cart(finest_level+1);
	for (int lev=0; lev<=finest_level; ++lev) {
		w0r_cart[lev].define(grids[lev], dmap[lev], 1, 1);
	}

	Put1dArrayOnCart(w0,w0r_cart,1,0,bcs_u,0);

	for (int lev=0; lev<=finest_level; ++lev) {

		// get references to the MultiFabs at level lev
		const MultiFab& vel_mf = vel[lev];
		MultiFab& radvel_mf = rad_vel[lev];
		MultiFab& circvel_mf = circ_vel[lev];
		const MultiFab& w0rcart_mf = w0r_cart[lev];
		const MultiFab& normal_mf = normal[lev];

#ifdef _OPENMP
#pragma omp parallel
#endif
		for ( MFIter mfi(vel_mf, true); mfi.isValid(); ++mfi ) {

			// Get the index space of the valid region
			const Box& tileBox = mfi.tilebox();

			// call fortran subroutine
			// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
			// lo/hi coordinates (including ghost cells), and/or the # of components
			// We will also pass "validBox", which specifies the "valid" region.

			make_velrc(ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
			           BL_TO_FORTRAN_3D(vel_mf[mfi]),
			           BL_TO_FORTRAN_3D(w0rcart_mf[mfi]),
			           BL_TO_FORTRAN_3D(normal_mf[mfi]),
			           BL_TO_FORTRAN_3D(radvel_mf[mfi]),
			           BL_TO_FORTRAN_3D(circvel_mf[mfi]));
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
			for ( MFIter mfi(state_mf, true); mfi.isValid(); ++mfi ) {

				// Get the index space of the valid region
				const Box& tileBox = mfi.tilebox();

				// call fortran subroutine
				// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
				// lo/hi coordinates (including ghost cells), and/or the # of components
				// We will also pass "validBox", which specifies the "valid" region.
				make_ad_excess(ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
				               BL_TO_FORTRAN_FAB(state_mf[mfi]),
				               BL_TO_FORTRAN_3D(ad_excess_mf[mfi]));
			}

		} else {

			const MultiFab& normal_mf = normal[lev];

#ifdef _OPENMP
#pragma omp parallel
#endif
			for ( MFIter mfi(state_mf, true); mfi.isValid(); ++mfi ) {

				// Get the index space of the valid region
				const Box& tileBox = mfi.tilebox();

				// call fortran subroutine
				// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
				// lo/hi coordinates (including ghost cells), and/or the # of components
				// We will also pass "validBox", which specifies the "valid" region.

				make_ad_excess_sphr(ARLIM_3D(tileBox.loVect()),
				                    ARLIM_3D(tileBox.hiVect()),
				                    BL_TO_FORTRAN_FAB(state_mf[mfi]),
				                    BL_TO_FORTRAN_3D(normal_mf[mfi]),
				                    BL_TO_FORTRAN_3D(ad_excess_mf[mfi]));
			}
		}
	}

	// average down and fill ghost cells
	AverageDown(ad_excess,0,1);
	FillPatch(t_old,ad_excess,ad_excess,ad_excess,0,0,1,0,bcs_f);
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
		for ( MFIter mfi(vel_mf, true); mfi.isValid(); ++mfi ) {

			// Get the index space of the valid region
			const Box& tileBox = mfi.tilebox();
			const Real* dx = geom[lev].CellSize();

			// call fortran subroutine
			// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
			// lo/hi coordinates (including ghost cells), and/or the # of components
			// We will also pass "validBox", which specifies the "valid" region.
			make_vorticity(ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
			               BL_TO_FORTRAN_3D(vel_mf[mfi]), dx,
			               BL_TO_FORTRAN_3D(vorticity_mf[mfi]), phys_bc.dataPtr());
		}


	}

	// average down and fill ghost cells
	AverageDown(vorticity,0,1);
	FillPatch(t_old,vorticity,vorticity,vorticity,0,0,1,0,bcs_f);
}

void
Maestro::MakeDeltaGamma (const Vector<MultiFab>& state,
                         const Vector<Real>& p0,
                         const Vector<MultiFab>& p0_cart,
                         const Vector<Real>& gamma1bar,
                         const Vector<MultiFab>& gamma1bar_cart,
                         Vector<MultiFab>& deltagamma)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::MakeDeltaGamma()",MakeDeltaGamma);

	for (int lev=0; lev<=finest_level; ++lev) {

		// get references to the MultiFabs at level lev
		const MultiFab& state_mf = state[lev];
		MultiFab& deltagamma_mf = deltagamma[lev];

        if (spherical == 0) {

    		// Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
    		for ( MFIter mfi(state_mf, true); mfi.isValid(); ++mfi ) {

    			// Get the index space of the valid region
    			const Box& tileBox = mfi.tilebox();

    			// call fortran subroutine
    			// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
    			// lo/hi coordinates (including ghost cells), and/or the # of components
    			// We will also pass "validBox", which specifies the "valid" region.
    			make_deltagamma(&lev,ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
    			                BL_TO_FORTRAN_FAB(state_mf[mfi]),
    			                p0.dataPtr(), gamma1bar.dataPtr(),
    			                BL_TO_FORTRAN_3D(deltagamma_mf[mfi]));
    		}

        } else {

            const MultiFab& p0cart_mf = p0_cart[lev];
            const MultiFab& gamma1barcart_mf = gamma1bar_cart[lev];

            // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
            for ( MFIter mfi(state_mf, true); mfi.isValid(); ++mfi ) {

                // Get the index space of the valid region
                const Box& tileBox = mfi.tilebox();

                // call fortran subroutine
                // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
                // lo/hi coordinates (including ghost cells), and/or the # of components
                // We will also pass "validBox", which specifies the "valid" region.
                make_deltagamma_sphr(ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
                                BL_TO_FORTRAN_FAB(state_mf[mfi]),
                                BL_TO_FORTRAN_3D(p0cart_mf[mfi]), BL_TO_FORTRAN_3D(gamma1barcart_mf[mfi]),
                                BL_TO_FORTRAN_3D(deltagamma_mf[mfi]));
            }

        }


	}

	// average down and fill ghost cells
	AverageDown(deltagamma,0,1);
	FillPatch(t_old,deltagamma,deltagamma,deltagamma,0,0,1,0,bcs_f);
}

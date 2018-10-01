
#include <Maestro.H>
#include <MaestroPlot.H>
#include <AMReX_buildInfo.H>

using namespace amrex;

// write plotfile to disk
void
Maestro::WritePlotFile (const int step,
                        const Real t_in,
                        const Vector<Real>& rho0_in,
                        const Vector<Real>& p0_in,
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

		// convert p0 to multi-D MultiFab
		Vector<MultiFab> p0_cart(finest_level+1);
		for (int lev=0; lev<=finest_level; ++lev) {
				p0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
		}
		Put1dArrayOnCart(p0_in,p0_cart,0,0);

		// convert rho0 to multi-D MultiFab
		Vector<MultiFab> rho0_cart(finest_level+1);
		for (int lev=0; lev<=finest_level; ++lev) {
				rho0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
		}
		Put1dArrayOnCart(rho0_in,rho0_cart,0,0);

		const auto& mf = PlotFileMF(p0_cart,rho0_cart,u_in,s_in,p0_in);
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

}


// get plotfile name
std::string
Maestro::PlotFileName (int lev) const
{
		return Concatenate(plot_base_name, lev, 7);
}

// put together a vector of multifabs for writing
Vector<const MultiFab*>
Maestro::PlotFileMF (const Vector<MultiFab>& p0_cart,
                     const Vector<MultiFab>& rho0_cart,
                     const Vector<MultiFab>& u_in,
                     Vector<MultiFab>& s_in,
                     const Vector<Real>& p0_in)
{
		// timer for profiling
		BL_PROFILE_VAR("Maestro::PlotFileMF()",PlotFileMF);

		// velocities (AMREX_SPACEDIM)
		// rho, rhoh, rhoX, tfromp, tfromh, Pi (Nscal+1)
		// X (NumSpec)
		// rho0, p0 (2)
		// deltaT (1)
		// w0 (AMREX_SPACEDIM)
		int nPlot = 2*AMREX_SPACEDIM + Nscal + NumSpec + 4;

		// MultiFab to hold plotfile data
		Vector<const MultiFab*> plot_mf;

		// temporary MultiFab to hold plotfile data
		Vector<MultiFab*> plot_mf_data(finest_level+1);

		// temporary MultiFab for calculations
		Vector<MultiFab> tempmf(finest_level+1);

		int dest_comp = 0;

		// build temporary MultiFab to hold plotfile data
		for (int i = 0; i <= finest_level; ++i) {
				plot_mf_data[i] = new MultiFab((s_in[i]).boxArray(),(s_in[i]).DistributionMap(),nPlot,0);
				tempmf[i].define(grids[i],dmap[i],AMREX_SPACEDIM,0);
		}

		// velocity
		for (int i = 0; i <= finest_level; ++i) {
				plot_mf_data[i]->copy((u_in[i]),0,dest_comp,AMREX_SPACEDIM);
		}
		dest_comp += AMREX_SPACEDIM;

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

		// rho0 and p0
		for (int i = 0; i <= finest_level; ++i) {
				plot_mf_data[i]->copy((rho0_cart[i]),0,dest_comp,1);
				plot_mf_data[i]->copy((  p0_cart[i]),0,dest_comp+1,1);
		}
		dest_comp += 2;

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
		} // spherical

		// w0
		Put1dArrayOnCart(w0,tempmf,1,1,bcs_u,0);
		for (int i = 0; i <= finest_level; ++i) {
				plot_mf_data[i]->copy((tempmf[i]),0,dest_comp,AMREX_SPACEDIM);
		}
		dest_comp += AMREX_SPACEDIM;

		// add plot_mf_data[i] to plot_mf
		for (int i = 0; i <= finest_level; ++i) {
				plot_mf.push_back(plot_mf_data[i]);
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
		// rho, rhoh, rhoX, tfromp, tfromh, Pi (Nscal+1)
		// X (NumSpec)
		// rho0, p0 (2)
		// deltaT (1)
		// w0 (AMREX_SPACEDIM)
		int nPlot = 2*AMREX_SPACEDIM + Nscal + NumSpec + 4;
		Vector<std::string> names(nPlot);

		int cnt = 0;

		// add velocities
		for (int i=0; i<AMREX_SPACEDIM; ++i) {
				std::string x = "vel";
				x += (120+i);
				names[cnt++] = x;
		}

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
		}

		names[cnt++] = "tfromp";
		names[cnt++] = "tfromh";
		names[cnt++] = "deltaT";
		names[cnt++] = "Pi";

		names[cnt++] = "rho0";
		names[cnt++] = "p0";

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

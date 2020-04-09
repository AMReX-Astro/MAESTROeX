
#include <Maestro.H>
#include <Maestro_F.H>
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
                        const Vector<MultiFab>& u_in,
                        Vector<MultiFab>& e,
                        const Vector<MultiFab>& f,
			const bool is_small)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::WritePlotFile()",WritePlotFile);

	// wallclock time
	const Real strt_total = ParallelDescriptor::second();

    std::string plotfilename = plot_base_name;

    if (step == 0) {
    	plotfilename += "_uold";
    } else if (step == 1) {
        plotfilename += "_umid";
    } else if (step == 2) {
        plotfilename += "_gphi";
    } else {
        plotfilename += "_unew";
    }

	int nPlot = 0;

    // make plot mf
    const Vector<MultiFab> dummy;

	// WriteMultiLevelPlotfile expects an array of step numbers
	Vector<int> step_array;
	step_array.resize(maxLevel()+1, step);

    if (step == 2) {
        const Vector<std::string> varnames = {"gphix", "gphiy", "gphiz"};
        const auto& mf = PlotFileMF(nPlot,t_in,dt_in,dummy,dummy,dummy,dummy,u_in,e,a,a, dummy);
        WriteMultiLevelPlotfile(plotfilename, finest_level+1, mf, varnames,
    	                        Geom(), t_in, step_array, refRatio());

    	for (int i = 0; i <= finest_level; ++i) {
    		delete mf[i];
    	}

    } else {
    	const auto& varnames = PlotFileVarNames(&nPlot);
        const auto& mf = PlotFileMF(nPlot,t_in,dt_in,dummy,dummy,dummy,dummy,u_in,e,a,a, dummy);

        WriteMultiLevelPlotfile(plotfilename, finest_level+1, mf, varnames,
    	                        Geom(), t_in, step_array, refRatio());

    	for (int i = 0; i <= finest_level; ++i) {
    		delete mf[i];
    	}
    }
}


// get plotfile name
void
Maestro::PlotFileName (const int lev, std::string* plotfilename)
{
	*plotfilename = Concatenate(*plotfilename, lev, 7);
}

// put together a vector of multifabs for writing
Vector<const MultiFab*>
Maestro::PlotFileMF (const int nPlot,
                     const Real t_in,
                     const Real dt_in,
                     const Vector<MultiFab>& a,
                     const Vector<MultiFab>& b,
                     const Vector<MultiFab>& c,
                     const Vector<MultiFab>& d,
                     const Vector<MultiFab>& u_in,
                     Vector<MultiFab>& e,
                     const Vector<Real>& f,
                     const Vector<Real>& g,
                     const Vector<MultiFab>& h)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::PlotFileMF()",PlotFileMF);

	// MultiFab to hold plotfile data
	Vector<const MultiFab*> plot_mf;

	// temporary MultiFab to hold plotfile data
	Vector<MultiFab*> plot_mf_data(finest_level+1);

	int dest_comp = 0;

	// build temporary MultiFab to hold plotfile data
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i] = new MultiFab((u_in[i]).boxArray(),(u_in[i]).DistributionMap(),AMREX_SPACEDIM,0);
	}

	// velocity
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy((u_in[i]),0,dest_comp,AMREX_SPACEDIM);
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
Maestro::PlotFileVarNames (int * nPlot) const
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::PlotFileVarNames()",PlotFileVarNames);

	(*nPlot) = AMREX_SPACEDIM;

	Vector<std::string> names(*nPlot);

	int cnt = 0;

	// add velocities
	for (int i=0; i<AMREX_SPACEDIM; ++i) {
		std::string x = "vel";
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
Maestro::WriteBuildInfo ()
{
        std::string PrettyLine = std::string(78, '=') + "\n";
	std::string OtherLine = std::string(78, '-') + "\n";
	std::string SkipSpace = std::string(8, ' ');

	// build information
	std::cout << PrettyLine;
	std::cout << " MAESTROeX Build Information\n";
	std::cout << PrettyLine;
	
	std::cout << "build date:    " << buildInfoGetBuildDate() << "\n";
	std::cout << "build machine: " << buildInfoGetBuildMachine() << "\n";
	std::cout << "build dir:     " << buildInfoGetBuildDir() << "\n";
	std::cout << "AMReX dir:     " << buildInfoGetAMReXDir() << "\n";
	
	std::cout << "\n";
	
	std::cout << "COMP:          " << buildInfoGetComp() << "\n";
	std::cout << "COMP version:  " << buildInfoGetCompVersion() << "\n";
	
	std::cout << "\n";
	
	std::cout << "C++ compiler:  " << buildInfoGetCXXName() << "\n";
	std::cout << "C++ flags:     " << buildInfoGetCXXFlags() << "\n";
	
	std::cout << "\n";
	
	std::cout << "Fortran comp:  " << buildInfoGetFName() << "\n";
	std::cout << "Fortran flags: " << buildInfoGetFFlags() << "\n";
	
	std::cout << "\n";
	
	std::cout << "Link flags:    " << buildInfoGetLinkFlags() << "\n";
	std::cout << "Libraries:     " << buildInfoGetLibraries() << "\n";
	
	std::cout << "\n";
	
	for (int n = 1; n <= buildInfoGetNumModules(); n++) {
	        std::cout << buildInfoGetModuleName(n) << ": " << buildInfoGetModuleVal(n) << "\n";
	}
	
	const char* githash1 = buildInfoGetGitHash(1);
	const char* githash2 = buildInfoGetGitHash(2);
	const char* githash3 = buildInfoGetGitHash(3);
	if (strlen(githash1) > 0) {
	        std::cout << "MAESTROeX git describe: " << githash1 << "\n";
	}
	if (strlen(githash2) > 0) {
	        std::cout << "AMReX git describe: " << githash2 << "\n";
	}
	if (strlen(githash3) > 0) {
	        std::cout << "Microphysics git describe: " << githash3 << "\n";
	}

	const char* buildgithash = buildInfoGetBuildGitHash();
	const char* buildgitname = buildInfoGetBuildGitName();
	if (strlen(buildgithash) > 0) {
	        std::cout << buildgitname << " git describe: " << buildgithash << "\n";
	}
	
	std::cout << "\n\n";
}

void
Maestro::MakeMagvel (const Vector<MultiFab>& vel,
                     Vector<MultiFab>& magvel)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::MakeMagvel()",MakeMagvel);

        Vector<std::array< MultiFab, AMREX_SPACEDIM > > w0mac(finest_level+1);

#if (AMREX_SPACEDIM == 3)
        if (spherical == 1) {
            for (int lev=0; lev<=finest_level; ++lev) {
                w0mac[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);
                w0mac[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
                w0mac[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1);
            }
            MakeW0mac(w0mac);
        }
#endif

	for (int lev=0; lev<=finest_level; ++lev) {

		// get references to the MultiFabs at level lev
		const MultiFab& vel_mf = vel[lev];

		// Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
		if (spherical == 0) {
#ifdef _OPENMP
#pragma omp parallel
#endif
		        for ( MFIter mfi(vel_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

				// Get the index space of the valid region
				const Box& tileBox = mfi.tilebox();

				const Array4<const Real> vel_arr = vel[lev].array(mfi);
				const Array4<const Real> w0_arr = w0_cart[lev].array(mfi);
				const Array4<Real> magvel_arr = magvel[lev].array(mfi);

				AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
#if (AMREX_SPACEDIM == 2)
					Real v_total = vel_arr(i,j,k,1) + 0.5 * (w0_arr(i,j,k,1) + w0_arr(i,j+1,k,1));
					magvel_arr(i,j,k) = sqrt(vel_arr(i,j,k,0)*vel_arr(i,j,k,0) +
								 v_total*v_total);
#else
					Real w_total = vel_arr(i,j,k,2) + 0.5 * (w0_arr(i,j,k,2) + w0_arr(i,j,k+1,2));
					magvel_arr(i,j,k) = sqrt(vel_arr(i,j,k,0)*vel_arr(i,j,k,0) +
								 vel_arr(i,j,k,1)*vel_arr(i,j,k,1) +
								 w_total*w_total);
#endif
				});
			}

		} else {

#ifdef _OPENMP
#pragma omp parallel
#endif
			for ( MFIter mfi(vel_mf, true); mfi.isValid(); ++mfi ) {

				// Get the index space of the valid region
				const Box& tileBox = mfi.tilebox();

				const Array4<const Real> vel_arr = vel[lev].array(mfi);
				const Array4<const Real> w0macx = w0mac[lev][0].array(mfi);
				const Array4<const Real> w0macy = w0mac[lev][1].array(mfi);
				const Array4<const Real> w0macz = w0mac[lev][2].array(mfi);
				const Array4<Real> magvel_arr = magvel[lev].array(mfi);

				AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
					Real u_total = vel_arr(i,j,k,0) + 0.5 * (w0macx(i,j,k) + w0macx(i+1,j,k));
					Real v_total = vel_arr(i,j,k,1) + 0.5 * (w0macy(i,j,k) + w0macy(i,j+1,k));
					Real w_total = vel_arr(i,j,k,2) + 0.5 * (w0macz(i,j,k) + w0macz(i,j,k+1));
					magvel_arr(i,j,k) = sqrt(u_total*u_total +
								 v_total*v_total + w_total*w_total);
				});
			}
		}
	}

	// average down and fill ghost cells
	AverageDown(magvel,0,1);
	FillPatch(t_old,magvel,magvel,magvel,0,0,1,0,bcs_f);
}


void
Maestro::MakeVelrc (const Vector<MultiFab>& vel,
                    const Vector<MultiFab>& w0rcart,
                    Vector<MultiFab>& rad_vel,
                    Vector<MultiFab>& circ_vel)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::MakeVelrc()",MakeVelrc);

	for (int lev=0; lev<=finest_level; ++lev) {

#ifdef _OPENMP
#pragma omp parallel
#endif
		for ( MFIter mfi(vel[lev], true); mfi.isValid(); ++mfi ) {

			// Get the index space of the valid region
			const Box& tileBox = mfi.tilebox();

			const Array4<const Real> vel_arr = vel[lev].array(mfi);
			const Array4<Real> radvel_arr = rad_vel[lev].array(mfi);
			const Array4<Real> circvel_arr = circ_vel[lev].array(mfi);
			const Array4<const Real> w0rcart_arr = w0rcart[lev].array(mfi);
			const Array4<const Real> normal_arr = normal[lev].array(mfi);

			AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
				circvel_arr(i,j,k) = 0.0;
				radvel_arr(i,j,k) = 0.0;

				for (auto n = 0; n < AMREX_SPACEDIM; ++n) {
				    radvel_arr(i,j,k) += vel_arr(i,j,k,n) * normal_arr(i,j,k,n);
				}

				for (auto n = 0; n < AMREX_SPACEDIM; ++n) {
				    Real circ_comp = vel_arr(i,j,k,n) - radvel_arr(i,j,k) * normal_arr(i,j,k,n);
				    circvel_arr(i,j,k) += circ_comp * circ_comp;
				}

				circvel_arr(i,j,k) = sqrt(circvel_arr(i,j,k));

				// add base state vel to get full radial velocity
				radvel_arr(i,j,k) += w0rcart_arr(i,j,k);
			});
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
    BL_PROFILE_VAR("Maestro::MakeAdExcess()", MakeAdExcess);

    const auto base_cutoff_density_loc = base_cutoff_density;

    for (int lev=0; lev<=finest_level; ++lev) {

        // create MultiFabs to hold pressure and gradient
        MultiFab pres_mf(grids[lev], dmap[lev], 1, 0);
        MultiFab nabla_ad_mf(grids[lev], dmap[lev], 1, 0);

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(state[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const Array4<const Real> state_arr = state[lev].array(mfi);
            const Array4<Real> ad_excess_arr = ad_excess[lev].array(mfi);
            const Array4<Real> pres = pres_mf.array(mfi);
            const Array4<Real> nabla_ad = nabla_ad_mf.array(mfi);
#if (AMREX_SPACEDIM == 3)
            const Array4<const Real> normal_arr = normal[lev].array(mfi);
#endif

            AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                eos_t eos_state;

                eos_state.rho   = state_arr(i,j,k,Rho);
                eos_state.T     = state_arr(i,j,k,Temp);
                for (auto comp = 0; comp < NumSpec; ++comp) {
                    eos_state.xn[comp] = state_arr(i,j,k,FirstSpec+comp)/eos_state.rho;
                }

                eos(eos_input_rt, eos_state);

                pres(i,j,k) = eos_state.p;
                // Print() << "pres = " << pres(i,j,k) << std::endl;

                Real chi_rho = eos_state.rho * eos_state.dpdr / eos_state.p;
                Real chi_t = eos_state.T * eos_state.dpdT / eos_state.p;
                nabla_ad(i,j,k) = (eos_state.gam1 - chi_rho) / (chi_t * eos_state.gam1);
            });

            const auto lo = tileBox.loVect3d();
            const auto hi = tileBox.hiVect3d();

            if (spherical == 0) {
                AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                    Real nabla = 0.0;

                    if (state_arr(i,j,k,Rho) > base_cutoff_density_loc) {
                        Real dtemp = 0.0;
                        Real dp = 0.0;
#if (AMREX_SPACEDIM == 2)
                        // forward difference
                        if (j == lo[1]) {
                            dtemp = state_arr(i,j+1,k,Temp) - state_arr(i,j,k,Temp);
                            dp = pres(i,j+1,k) - pres(i,j,k);
                        // backward difference
                        } else if (j == hi[1]) {
                            dtemp = state_arr(i,j,k,Temp) - state_arr(i,j-1,k,Temp);
                            dp = pres(i,j,k) - pres(i,j-1,k);
                        // centered difference
                        } else {
                            dtemp = state_arr(i,j+1,k,Temp) - state_arr(i,j-1,k,Temp);
                            dp = pres(i,j+1,k) - pres(i,j-1,k);
                        }
#else 
                        // forward difference
                        if (k == lo[2]) {
                            dtemp = state_arr(i,j,k+1,Temp) - state_arr(i,j,k,Temp);
                            dp = pres(i,j,k+1) - pres(i,j,k);
                        // backward difference
                        } else if (k == hi[2]) {
                            dtemp = state_arr(i,j,k,Temp) - state_arr(i,j,k-1,Temp);
                            dp = pres(i,j,k) - pres(i,j,k-1);
                        // centered difference
                        } else {
                            dtemp = state_arr(i,j,k+1,Temp) - state_arr(i,j,k-1,Temp);
                            dp = pres(i,j,k+1) - pres(i,j,k-1);
                        }
#endif
                        // prevent Inf
                        if (dp * state_arr(i,j,k,Temp) == 0.0) {
                            nabla = std::numeric_limits<Real>::min();
                        } else {
                            nabla = pres(i,j,k) * dtemp / (dp * state_arr(i,j,k,Temp));
                        }
                    }

                    ad_excess_arr(i,j,k) = nabla - nabla_ad(i,j,k);
                });
            } else {
#if (AMREX_SPACEDIM == 3)
                RealVector dtemp_vec(AMREX_SPACEDIM, 0.0);
                RealVector dp_vec(AMREX_SPACEDIM, 0.0);

                Real * AMREX_RESTRICT dtemp = dtemp_vec.dataPtr();
                Real * AMREX_RESTRICT dp = dp_vec.dataPtr();

                AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                    Real nabla = 0.0;

                    if (state_arr(i,j,k,Rho) > base_cutoff_density_loc) {
                        // compute gradient
                        // forward difference
                        if (i == lo[0]) {
                            dtemp[0] = state_arr(i+1,j,k,Temp) - state_arr(i,j,k,Temp);
                            dp[0] = pres(i+1,j,k) - pres(i,j,k);
                        // backward difference
                        } else if (i == hi[0]) {
                            dtemp[0] = state_arr(i,j,k,Temp) - state_arr(i-1,j,k,Temp);
                            dp[0] = pres(i,j,k) - pres(i-1,j,k);
                        // centered difference
                        } else {
                            dtemp[0] = state_arr(i+1,j,k,Temp) - state_arr(i-1,j,k,Temp);
                            dp[0] = pres(i+1,j,k) - pres(i-1,j,k);
                        }
                        // forward difference
                        if (j == lo[1]) {
                            dtemp[1] = state_arr(i,j+1,k,Temp) - state_arr(i,j,k,Temp);
                            dp[1] = pres(i,j+1,k) - pres(i,j,k);
                        // backward difference
                        } else if (j == hi[1]) {
                            dtemp[1] = state_arr(i,j,k,Temp) - state_arr(i,j-1,k,Temp);
                            dp[1] = pres(i,j,k) - pres(i,j-1,k);
                        // centered difference
                        } else {
                            dtemp[1] = state_arr(i,j+1,k,Temp) - state_arr(i,j-1,k,Temp);
                            dp[1] = pres(i,j+1,k) - pres(i,j-1,k);
                        }
                        // forward difference
                        if (k == lo[2]) {
                            dtemp[2] = state_arr(i,j,k+1,Temp) - state_arr(i,j,k,Temp);
                            dp[2] = pres(i,j,k+1) - pres(i,j,k);
                        // backward difference
                        } else if (k == hi[2]) {
                            dtemp[2] = state_arr(i,j,k,Temp) - state_arr(i,j,k-1,Temp);
                            dp[2] = pres(i,j,k) - pres(i,j,k-1);
                        // centered difference
                        } else {
                            dtemp[2] = state_arr(i,j,k+1,Temp) - state_arr(i,j,k-1,Temp);
                            dp[2] = pres(i,j,k+1) - pres(i,j,k-1);
                        }

                        Real dp_dot = 0.0;
                        Real dtemp_dot = 0.0;
                        for (auto c = 0; c < AMREX_SPACEDIM; ++c) {
                            dp_dot += dp[c] * normal_arr(i,j,k,c);
                            dtemp_dot += dtemp[c] * normal_arr(i,j,k,c);
                        }

                        // prevent Inf
                        if (dp_dot * state_arr(i,j,k,Temp) == 0.0) {
                            nabla = std::numeric_limits<Real>::min();
                        } else {
                            nabla = pres(i,j,k) * dtemp_dot / (dp_dot * state_arr(i,j,k,Temp));
                        }
                    }

                    ad_excess_arr(i,j,k) = nabla - nabla_ad(i,j,k);
                });
#endif
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

        const Real* dx = geom[lev].CellSize();
        const Box& domainBox = geom[lev].Domain();

        const Real hx = dx[0];
        const Real hy = dx[1];
#if (AMREX_SPACEDIM == 3)
        const Real hz = dx[2];
#endif
        const int ilo = domainBox.loVect()[0];
        const int ihi = domainBox.hiVect()[0];
        const int jlo = domainBox.loVect()[1];
        const int jhi = domainBox.hiVect()[1];
#if (AMREX_SPACEDIM == 3)
        const int klo = domainBox.loVect()[2];
        const int khi = domainBox.hiVect()[2];
#endif

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(vel_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            Array4<const Real> const u = vel[lev].array(mfi);
            Array4<Real> const vort = vorticity[lev].array(mfi);
            GpuArray<int,AMREX_SPACEDIM*2> physbc;
            for (int n = 0; n < AMREX_SPACEDIM*2; ++n) {
                physbc[n] = phys_bc[n];
            } 

#if (AMREX_SPACEDIM == 2)

            AMREX_PARALLEL_FOR_3D(tileBox, i, j, k,
            {
                Real vx = 0.5*(u(i+1,j,k,1)-u(i-1,j,k,1))/hx;
                Real uy = 0.5*(u(i,j+1,k,0)-u(i,j-1,k,0))/hy;

                if (i == ilo && 
                    (physbc[0] == Inflow || 
                     physbc[0] == SlipWall || 
                     physbc[0] == NoSlipWall)) 
                {
                    vx = (u(i+1,j,k,1) + 3.0*u(i,j,k,1) - 
                          4.0*u(i-1,j,k,1)) / hx;
                    uy = 0.5 * (u(i,j+1,k,0) - u(i,j-1,k,0)) / hy;

                } else if (i == ihi+1 &&
                        (physbc[AMREX_SPACEDIM] == Inflow || 
                        physbc[AMREX_SPACEDIM] == SlipWall || 
                        physbc[AMREX_SPACEDIM] == NoSlipWall))
                {
                    vx = -(u(i-1,j,k,1) + 3.0*u(i,j,k,1) - 
                         4.0*u(i+1,j,k,1)) / hx;
                    uy = 0.5 * (u(i,j+1,k,0) - u(i,j-1,k,0)) / hy;
                }

                if (j == jlo &&
                    (physbc[1] == Inflow || 
                     physbc[1] == SlipWall || 
                     physbc[1] == NoSlipWall))
                {
                    vx = 0.5 * (u(i+1,j,k,1) - u(i-1,j,k,0)) / hx;
                    uy = (u(i,j+1,k,0) + 3.0*u(i,j,k,0) - 
                         4.0*u(i,j-1,k,0)) / hy;

                } else if (j == jhi+1 && 
                           (physbc[AMREX_SPACEDIM+1] == Inflow || 
                            physbc[AMREX_SPACEDIM+1] == SlipWall || 
                            physbc[AMREX_SPACEDIM+1] == NoSlipWall))
                {
                    vx = 0.5 * (u(i+1,j,k,1) - u(i-1,j,k,1)) / hx;
                    uy = -(u(i,j-1,k,0) + 3.0*u(i,j,k,0) - 
                         4.0*u(i,j+1,k,0)) / hy;
                }

                vort(i,j,k) = vx - uy;
            });

#else 
            AMREX_PARALLEL_FOR_3D(tileBox, i, j, k,
            {
                Real uy = 0.5*(u(i,j+1,k,0)-u(i,j-1,k,0))/hy;
                Real uz = 0.5*(u(i,j,k+1,0)-u(i,j,k-1,0))/hz;
                Real vx = 0.5*(u(i+1,j,k,1)-u(i-1,j,k,1))/hx;
                Real vz = 0.5*(u(i,j,k+1,1)-u(i,j,k-1,1))/hz;
                Real wx = 0.5*(u(i+1,j,k,2)-u(i-1,j,k,2))/hx;
                Real wy = 0.5*(u(i,j+1,k,2)-u(i,j-1,k,2))/hy;

                bool fix_lo_x = (physbc[0] == Inflow || 
                                 physbc[0] == NoSlipWall);
                bool fix_hi_x = (physbc[AMREX_SPACEDIM] == Inflow || 
                                 physbc[AMREX_SPACEDIM] == NoSlipWall);

                bool fix_lo_y = (physbc[1] == Inflow || 
                                 physbc[1] == NoSlipWall);
                bool fix_hi_y = (physbc[AMREX_SPACEDIM+1] == Inflow ||
                                 physbc[AMREX_SPACEDIM+1] == NoSlipWall);

                bool fix_lo_z = (physbc[2] == Inflow || 
                                 physbc[2] == NoSlipWall);
                bool fix_hi_z = (physbc[AMREX_SPACEDIM+2] == Inflow ||
                                 physbc[AMREX_SPACEDIM+2] == NoSlipWall);

                // First do all the faces
                if (fix_lo_x && i == ilo) {
                    vx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    wx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                } else if (fix_hi_x && i == ihi+1) {
                    vx = -(u(i-1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i+1,j,k,1))/(3.0*hx);
                    wx = -(u(i-1,j,k,2)+3.0*u(i,j,k,2)-4.0*u(i+1,j,k,2))/(3.0*hx);
                }

                if (fix_lo_y && j == jlo) {
                    uy = (u(i,j+1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j-1,k,0))/(3.0*hy);
                    wy = (u(i,j+1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j-1,k,2))/(3.0*hy);
                } else if (fix_hi_y && j == jhi+1) {
                    uy = -(u(i,j-1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j+1,k,0))/(3.0*hy);
                    wy = -(u(i,j-1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j+1,k,2))/(3.0*hy);
                }

                if (fix_lo_z && k == klo) {
                    uz = (u(i,j,k+1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k-1,0))/(3.0*hz);
                    vz = (u(i,j,k+1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k-1,1))/(3.0*hz);
                }

                if (fix_hi_z && k == khi+1) {
                    uz = -(u(i,j,k-1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k+1,0))/(3.0*hz);
                    vz = -(u(i,j,k-1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k+1,1))/(3.0*hz);
                }

                // Next do all the edges
                if (fix_lo_x && fix_lo_y && i == ilo && j == jlo) {
                    vx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    wx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    uy = (u(i,j+1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j-1,k,0))/(3.0*hy);
                    wy = (u(i,j+1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j-1,k,2))/(3.0*hy);
                }

                if (fix_hi_x && fix_lo_y && i == ihi+1 && j == jlo) {
                    vx = -(u(i-1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i+1,j,k,1))/(3.0*hx);
                    wx = -(u(i-1,j,k,2)+3.0*u(i,j,k,2)-4.0*u(i+1,j,k,2))/(3.0*hx);
                    uy = (u(i,j+1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j-1,k,0))/(3.0*hy);
                    wy = (u(i,j+1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j-1,k,2))/(3.0*hy);
                }

                if (fix_lo_x && fix_hi_y && i == ilo && j == jhi+1) {
                    vx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    wx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    uy = -(u(i,j-1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j+1,k,0))/(3.0*hy);
                    wy = -(u(i,j-1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j+1,k,2))/(3.0*hy);
                }

                if (fix_lo_x && fix_lo_z && i == ilo && k == klo) {
                    vx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    wx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    uz = (u(i,j,k+1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k-1,0))/(3.0*hz);
                    vz = (u(i,j,k+1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k-1,1))/(3.0*hz);
                }

                if (fix_hi_x && fix_lo_z && i == ihi+1 && k == klo) {
                    vx = -(u(i-1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i+1,j,k,1))/(3.0*hx);
                    wx = -(u(i-1,j,k,2)+3.0*u(i,j,k,2)-4.0*u(i+1,j,k,2))/(3.0*hx);
                    uz = (u(i,j,k+1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k-1,0))/(3.0*hz);
                    vz = (u(i,j,k+1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k-1,1))/(3.0*hz);
                }

                if (fix_lo_x && fix_hi_z && i == ilo && k == khi+1) {
                    vx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    wx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    uz = -(u(i,j,k-1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k+1,0))/(3.0*hz);
                    vz = -(u(i,j,k-1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k+1,1))/(3.0*hz);
                }

                if (fix_hi_x && fix_hi_z && i == ihi+1 && k == khi+1) {
                    vx = -(u(i-1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i+1,j,k,1))/(3.0*hx);
                    wx = -(u(i-1,j,k,2)+3.0*u(i,j,k,2)-4.0*u(i+1,j,k,2))/(3.0*hx);
                    uz = -(u(i,j,k-1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k+1,0))/(3.0*hz);
                    vz = -(u(i,j,k-1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k+1,1))/(3.0*hz);
                }

                if (fix_lo_y && fix_lo_z && j == jlo && k == klo) {
                    uy = (u(i,j+1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j-1,k,0))/(3.0*hy);
                    wy = (u(i,j+1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j-1,k,2))/(3.0*hy);
                    uz = (u(i,j,k+1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k-1,0))/(3.0*hz);
                    vz = (u(i,j,k+1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k-1,1))/(3.0*hz);
                }

                if (fix_hi_y && fix_lo_z && j == jhi+1 && k == klo) {
                    uy = -(u(i,j-1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j+1,k,0))/(3.0*hy);
                    wy = -(u(i,j-1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j+1,k,2))/(3.0*hy);
                    uz = (u(i,j,k+1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k-1,0))/(3.0*hz);
                    vz = (u(i,j,k+1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k-1,1))/(3.0*hz);
                }

                if (fix_lo_y && fix_hi_z && j == jlo && k == khi+1) {
                    uy = (u(i,j+1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j-1,k,0))/(3.0*hy);
                    wy = (u(i,j+1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j-1,k,2))/(3.0*hy);
                    uz = -(u(i,j,k-1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k+1,0))/(3.0*hz);
                    vz = -(u(i,j,k-1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k+1,1))/(3.0*hz);
                }

                if (fix_hi_y && fix_hi_z && j == jhi+1 && k == khi+1) {
                    uy = -(u(i,j-1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j+1,k,0))/(3.0*hy);
                    wy = -(u(i,j-1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j+1,k,2))/(3.0*hy);
                    uz = -(u(i,j,k-1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k+1,0))/(3.0*hz);
                    vz = -(u(i,j,k-1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k+1,1))/(3.0*hz);
                }
                
                // Finally do all the corners
                if (fix_lo_x && fix_lo_y && fix_lo_z && 
                    i == ilo && j == jlo && k == klo) {
                    vx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    wx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    uy = (u(i,j+1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j-1,k,0))/(3.0*hy);
                    wy = (u(i,j+1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j-1,k,2))/(3.0*hy);
                    uz = (u(i,j,k+1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k-1,0))/(3.0*hz);
                    vz = (u(i,j,k+1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k-1,1))/(3.0*hz);
                }

                if (fix_hi_x && fix_lo_y && fix_lo_z &&
                    i == ihi+1 && j == jlo && k == klo) {
                    vx = -(u(i-1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i+1,j,k,1))/(3.0*hx);
                    wx = -(u(i-1,j,k,2)+3.0*u(i,j,k,2)-4.0*u(i+1,j,k,2))/(3.0*hx);
                    uy = (u(i,j+1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j-1,k,0))/(3.0*hy);
                    wy = (u(i,j+1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j-1,k,2))/(3.0*hy);
                    uz = (u(i,j,k+1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k-1,0))/(3.0*hz);
                    vz = (u(i,j,k+1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k-1,1))/(3.0*hz);
                }

                if (fix_lo_x && fix_hi_y && fix_lo_z &&
                    i == ilo && j == jhi+1 && k == klo) {
                    vx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    wx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    uy = -(u(i,j-1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j+1,k,0))/(3.0*hy);
                    wy = -(u(i,j-1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j+1,k,2))/(3.0*hy);
                    uz = (u(i,j,k+1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k-1,0))/(3.0*hz);
                    vz = (u(i,j,k+1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k-1,1))/(3.0*hz);
                }

                if (fix_hi_x && fix_hi_y && fix_lo_z &&
                    i == ihi+1 && j == jhi+1 && k == klo) {
                    vx = -(u(i-1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i+1,j,k,1))/(3.0*hx);
                    wx = -(u(i-1,j,k,2)+3.0*u(i,j,k,2)-4.0*u(i+1,j,k,2))/(3.0*hx);
                    uy = -(u(i,j-1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j+1,k,0))/(3.0*hy);
                    wy = -(u(i,j-1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j+1,k,2))/(3.0*hy);
                    uz = (u(i,j,k+1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k-1,0))/(3.0*hz);
                    vz = (u(i,j,k+1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k-1,1))/(3.0*hz);
                }

                if (fix_lo_x && fix_lo_y && fix_hi_z &&
                    i == ilo && j == jlo && k == khi+1) {
                    vx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    wx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    uy = (u(i,j+1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j-1,k,0))/(3.0*hy);
                    wy = (u(i,j+1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j-1,k,2))/(3.0*hy);
                    uz = -(u(i,j,k-1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k+1,0))/(3.0*hz);
                    vz = -(u(i,j,k-1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k+1,1))/(3.0*hz);
                }

                if (fix_hi_x && fix_lo_y && fix_hi_z &&
                    i == ihi+1 && j == jlo && k == khi+1) {
                    vx = -(u(i-1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i+1,j,k,1))/(3.0*hx);
                    wx = -(u(i-1,j,k,2)+3.0*u(i,j,k,2)-4.0*u(i+1,j,k,2))/(3.0*hx);
                    uy = (u(i,j+1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j-1,k,0))/(3.0*hy);
                    wy = (u(i,j+1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j-1,k,2))/(3.0*hy);
                    uz = -(u(i,j,k-1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k+1,0))/(3.0*hz);
                    vz = -(u(i,j,k-1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k+1,1))/(3.0*hz);
                }

                if (fix_lo_x && fix_hi_y && fix_hi_z &&
                    i == ilo && j == jhi+1 && k == khi+1) {
                    vx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    wx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    uy = -(u(i,j-1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j+1,k,0))/(3.0*hy);
                    wy = -(u(i,j-1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j+1,k,2))/(3.0*hy);
                    uz = -(u(i,j,k-1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k+1,0))/(3.0*hz);
                    vz = -(u(i,j,k-1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k+1,1))/(3.0*hz);
                }

                if (fix_hi_x && fix_hi_y && fix_hi_z &&
                    i == ihi+1 && j == jhi+1 && k == khi+1) {
                    vx = -(u(i-1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i+1,j,k,1))/(3.0*hx);
                    wx = -(u(i-1,j,k,2)+3.0*u(i,j,k,2)-4.0*u(i+1,j,k,2))/(3.0*hx);
                    uy = -(u(i,j-1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j+1,k,0))/(3.0*hy);
                    wy = -(u(i,j-1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j+1,k,2))/(3.0*hy);
                    uz = -(u(i,j,k-1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k+1,0))/(3.0*hz);
                    vz = -(u(i,j,k-1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k+1,1))/(3.0*hz);
                }

                vort(i,j,k) = sqrt((wy-vz)*(wy-vz)+
                    (uz-wx)*(uz-wx)+(vx-uy)*(vx-uy));
            });
#endif
        }
    }

    // average down and fill ghost cells
    AverageDown(vorticity,0,1);
    FillPatch(t_old,vorticity,vorticity,vorticity,0,0,1,0,bcs_f);
}

void
Maestro::MakeDeltaGamma (const Vector<MultiFab>& state,
                         const RealVector& p0,
                         const Vector<MultiFab>& p0_cart,
                         const RealVector& gamma1bar,
                         const Vector<MultiFab>& gamma1bar_cart,
                         Vector<MultiFab>& deltagamma)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeDeltaGamma()", MakeDeltaGamma);

    const auto use_pprime_in_tfromp_loc = use_pprime_in_tfromp;

    for (int lev=0; lev<=finest_level; ++lev) {

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(state[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const Array4<const Real> state_arr = state[lev].array(mfi);
            const Array4<const Real> p0_arr = p0_cart[lev].array(mfi);
            const Array4<const Real> gamma1bar_arr = gamma1bar_cart[lev].array(mfi);
            const Array4<Real> deltagamma_arr = deltagamma[lev].array(mfi);

            AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                eos_t eos_state;

                eos_state.rho   = state_arr(i,j,k,Rho);
                eos_state.T     = state_arr(i,j,k,Temp);
                if (use_pprime_in_tfromp_loc) {
                    eos_state.p     = p0_arr(i,j,k) + state_arr(i,j,k,Pi);
                } else {
                    eos_state.p     = p0_arr(i,j,k);
                }

                for (auto comp = 0; comp < NumSpec; ++comp) {
                    eos_state.xn[comp] = state_arr(i,j,k,FirstSpec+comp)/eos_state.rho;
                }

                eos(eos_input_rp, eos_state);

                deltagamma_arr(i,j,k) = eos_state.gam1 - gamma1bar_arr(i,j,k);
            });
        }
    }

    // average down and fill ghost cells
    AverageDown(deltagamma, 0, 1);
    FillPatch(t_old, deltagamma, deltagamma, deltagamma, 0, 0, 1, 0, bcs_f);
}


void
Maestro::MakeDivw0 (const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& w0mac,
                    Vector<MultiFab>& divw0)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeDivw0()", MakeDivw0);

    for (int lev=0; lev<=finest_level; ++lev) {

        if (spherical == 0) {

            // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
            for ( MFIter mfi(divw0[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

                // Get the index space of the valid region
                const Box& tileBox = mfi.tilebox();
                const auto dx = geom[lev].CellSizeArray();
                
                const Array4<const Real> w0_arr = w0_cart[lev].array(mfi);
                const Array4<Real> divw0_arr = divw0[lev].array(mfi);

                AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
#if (AMREX_SPACEDIM == 2)
                    divw0_arr(i,j,k) = (w0_arr(i,j+1,k,1) - w0_arr(i,j,k,1)) / dx[1];
#else
                    divw0_arr(i,j,k) = (w0_arr(i,j,k+1,2) - w0_arr(i,j,k,2)) / dx[2];
#endif
                });
            }

        } else {

            // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
            for ( MFIter mfi(divw0[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

                // Get the index space of the valid region
                const Box& tileBox = mfi.tilebox();
                const auto dx = geom[lev].CellSizeArray();
                
                const Array4<const Real> w0macx = w0mac[lev][0].array(mfi);
                const Array4<const Real> w0macy = w0mac[lev][1].array(mfi);
                const Array4<const Real> w0macz = w0mac[lev][2].array(mfi);
                const Array4<Real> divw0_arr = divw0[lev].array(mfi);

                AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                    divw0_arr(i,j,k) = (w0macx(i+1,j,k) - w0macx(i,j,k)) / dx[0] + 
                        (w0macy(i,j+1,k) - w0macy(i,j,k)) / dx[1] + 
                        (w0macz(i,j,k+1) - w0macz(i,j,k)) / dx[2];
                });
            }
        }
    }

    // average down and fill ghost cells
    AverageDown(divw0,0,1);
    FillPatch(t_old,divw0,divw0,divw0,0,0,1,0,bcs_f);
}

void
Maestro::MakePiDivu (const Vector<MultiFab>& vel,
                     const Vector<MultiFab>& state,
                     Vector<MultiFab>& pidivu)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakePiDivu()",MakePiDivu);

    for (int lev=0; lev<=finest_level; ++lev) {

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(pidivu[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const auto dx = geom[lev].CellSizeArray();

            const Array4<const Real> vel_arr = vel[lev].array(mfi);
            const Array4<const Real> pi_cc = state[lev].array(mfi, Pi);
            const Array4<Real> pidivu_arr = pidivu[lev].array(mfi);

            AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                pidivu_arr(i,j,k) = pi_cc(i,j,k) * 0.5 * (
                    (vel_arr(i+1,j,k,0) - vel_arr(i-1,j,k,0))/dx[0]
                  + (vel_arr(i,j+1,k,1) - vel_arr(i,j-1,k,1))/dx[1]
#if (AMREX_SPACEDIM == 3)
                  + (vel_arr(i,j,k+1,2) - vel_arr(i,j,k-1,2))/dx[2]
#endif
                    );
            });
        }
    }

    // average down and fill ghost cells
    AverageDown(pidivu, 0, 1);
    FillPatch(t_old, pidivu, pidivu, pidivu, 0, 0, 1, 0, bcs_f);
}

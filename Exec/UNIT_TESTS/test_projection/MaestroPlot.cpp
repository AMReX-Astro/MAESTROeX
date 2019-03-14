
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
                        const Vector<MultiFab>& u_in,
                        Vector<MultiFab>& e,
                        const Vector<MultiFab>& f)
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

#ifdef _OPENMP
#pragma omp parallel
#endif
			for ( MFIter mfi(vel_mf, true); mfi.isValid(); ++mfi ) {

				// Get the index space of the valid region
				const Box& tileBox = mfi.tilebox();

                                MultiFab& w0macx_mf = w0mac[lev][0];
                                MultiFab& w0macy_mf = w0mac[lev][1];
                                MultiFab& w0macz_mf = w0mac[lev][2];

				// call fortran subroutine
				// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
				// lo/hi coordinates (including ghost cells), and/or the # of components
				// We will also pass "validBox", which specifies the "valid" region.
				make_magvel_sphr(ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
				                 BL_TO_FORTRAN_3D(vel_mf[mfi]),
                                                 BL_TO_FORTRAN_3D(w0macx_mf[mfi]),
                                                 BL_TO_FORTRAN_3D(w0macy_mf[mfi]),
                                                 BL_TO_FORTRAN_3D(w0macz_mf[mfi]),
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
                    const Vector<MultiFab>& w0rcart,
                    Vector<MultiFab>& rad_vel,
                    Vector<MultiFab>& circ_vel)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::MakeVelrc()",MakeVelrc);

	for (int lev=0; lev<=finest_level; ++lev) {

		// get references to the MultiFabs at level lev
		const MultiFab& vel_mf = vel[lev];
		MultiFab& radvel_mf = rad_vel[lev];
		MultiFab& circvel_mf = circ_vel[lev];
		const MultiFab& w0rcart_mf = w0rcart[lev];
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

void
Maestro::MakeDivw0 (const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& w0mac,
                    Vector<MultiFab>& divw0)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::MakeDivw0()",MakeDivw0);

	for (int lev=0; lev<=finest_level; ++lev) {

		// get references to the MultiFabs at level lev
		MultiFab& divw0_mf = divw0[lev];

		if (spherical == 0) {

			// Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
			for ( MFIter mfi(divw0_mf, true); mfi.isValid(); ++mfi ) {

				// Get the index space of the valid region
				const Box& tileBox = mfi.tilebox();
				const Real* dx = geom[lev].CellSize();

				// call fortran subroutine
				// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
				// lo/hi coordinates (including ghost cells), and/or the # of components
				// We will also pass "validBox", which specifies the "valid" region.
				make_divw0(&lev,ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
				           w0.dataPtr(), dx,
				           BL_TO_FORTRAN_3D(divw0_mf[mfi]));
			}

		} else {

			const MultiFab& w0macx_mf = w0mac[lev][0];
			const MultiFab& w0macy_mf = w0mac[lev][1];
			const MultiFab& w0macz_mf = w0mac[lev][2];

			// Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
			for ( MFIter mfi(divw0_mf, true); mfi.isValid(); ++mfi ) {

				// Get the index space of the valid region
				const Box& tileBox = mfi.tilebox();
				const Real* dx = geom[lev].CellSize();

				// call fortran subroutine
				// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
				// lo/hi coordinates (including ghost cells), and/or the # of components
				// We will also pass "validBox", which specifies the "valid" region.
				make_divw0_sphr(ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
				                BL_TO_FORTRAN_3D(w0macx_mf[mfi]),
				                BL_TO_FORTRAN_3D(w0macy_mf[mfi]),
				                BL_TO_FORTRAN_3D(w0macz_mf[mfi]), dx,
				                BL_TO_FORTRAN_3D(divw0_mf[mfi]));
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
                     // const Vector<MultiFab>& pi_cc,
                     Vector<MultiFab>& pidivu)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::MakePiDivu()",MakePiDivu);

	for (int lev=0; lev<=finest_level; ++lev) {

		// get references to the MultiFabs at level lev
		const MultiFab& vel_mf = vel[lev];
		const MultiFab& state_mf = state[lev];
		// const MultiFab& pi_mf = pi_cc[lev];
		MultiFab& pidivu_mf = pidivu[lev];

		// Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
		for ( MFIter mfi(pidivu_mf, true); mfi.isValid(); ++mfi ) {

			// Get the index space of the valid region
			const Box& tileBox = mfi.tilebox();
			const Real* dx = geom[lev].CellSize();

			// call fortran subroutine
			// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
			// lo/hi coordinates (including ghost cells), and/or the # of components
			// We will also pass "validBox", which specifies the "valid" region.
			make_pidivu(ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
			            BL_TO_FORTRAN_3D(vel_mf[mfi]), dx,
			            // BL_TO_FORTRAN_3D(pi_mf[mfi]),
			            BL_TO_FORTRAN_FAB(state_mf[mfi]),
			            BL_TO_FORTRAN_3D(pidivu_mf[mfi]));
		}

	}

	// average down and fill ghost cells
	AverageDown(pidivu,0,1);
	FillPatch(t_old,pidivu,pidivu,pidivu,0,0,1,0,bcs_f);
}

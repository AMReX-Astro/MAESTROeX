
#include <AMReX_buildInfo.H>
#include <Maestro.H>
#include <MaestroPlot.H>
#include <unistd.h>  // getcwd

using namespace amrex;

// write plotfile to disk
void Maestro::WritePlotFile(const int step, const Real t_in, const Real dt_in,
                            const BaseState<Real>& a, const BaseState<Real>& b,
                            const BaseState<Real>& c, const BaseState<Real>& d,
                            const Vector<MultiFab>& u_in, Vector<MultiFab>& e,
                            const Vector<MultiFab>& f, const bool is_small) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::WritePlotFile()", WritePlotFile);

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
    step_array.resize(maxLevel() + 1, step);

    if (step == 2) {
        const Vector<std::string> varnames = {"gphix", "gphiy", "gphiz"};
        const auto& mf = PlotFileMF(nPlot, t_in, dt_in, dummy, dummy, dummy,
                                    dummy, u_in, e, b, b, dummy);
        WriteMultiLevelPlotfile(plotfilename, finest_level + 1, mf, varnames,
                                Geom(), t_in, step_array, refRatio());

        for (int i = 0; i <= finest_level; ++i) {
            delete mf[i];
        }

    } else {
        const auto& varnames = PlotFileVarNames(&nPlot);
        const auto& mf = PlotFileMF(nPlot, t_in, dt_in, dummy, dummy, dummy,
                                    dummy, u_in, e, b, b, dummy);

        WriteMultiLevelPlotfile(plotfilename, finest_level + 1, mf, varnames,
                                Geom(), t_in, step_array, refRatio());

        for (int i = 0; i <= finest_level; ++i) {
            delete mf[i];
        }
    }
}

// get plotfile name
void Maestro::PlotFileName(const int lev, std::string* plotfilename) {
    *plotfilename = Concatenate(*plotfilename, lev, 7);
}

// put together a vector of multifabs for writing
Vector<const MultiFab*> Maestro::PlotFileMF(
    const int nPlot, const Real t_in, const Real dt_in,
    const Vector<MultiFab>& a, const Vector<MultiFab>& b,
    const Vector<MultiFab>& c, const Vector<MultiFab>& d,
    const Vector<MultiFab>& u_in, Vector<MultiFab>& e, const BaseState<Real>& f,
    const BaseState<Real>& g, const Vector<MultiFab>& h) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::PlotFileMF()", PlotFileMF);

    // MultiFab to hold plotfile data
    Vector<const MultiFab*> plot_mf;

    // temporary MultiFab to hold plotfile data
    Vector<MultiFab*> plot_mf_data(finest_level + 1);

    int dest_comp = 0;

    // build temporary MultiFab to hold plotfile data
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i] =
            new MultiFab((u_in[i]).boxArray(), (u_in[i]).DistributionMap(),
                         AMREX_SPACEDIM, 0);
    }

    // velocity
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i]->copy((u_in[i]), 0, dest_comp, AMREX_SPACEDIM);
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
Vector<std::string> Maestro::PlotFileVarNames(int* nPlot) const {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::PlotFileVarNames()", PlotFileVarNames);

    (*nPlot) = AMREX_SPACEDIM;

    Vector<std::string> names(*nPlot);

    int cnt = 0;

    // add velocities
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        std::string x = "vel";
        x += (120 + i);
        names[cnt++] = x;
    }

    return names;
}

void Maestro::WriteJobInfo(const std::string& dir) const {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::WriteJobInfo()", WriteJobInfo);

    if (ParallelDescriptor::IOProcessor()) {
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

        jobInfoFile << "number of MPI processes: "
                    << ParallelDescriptor::NProcs() << "\n";
#ifdef _OPENMP
        jobInfoFile << "number of threads:       " << omp_get_max_threads()
                    << "\n";

        jobInfoFile << "tile size: ";
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            jobInfoFile << FabArrayBase::mfiter_tile_size[d] << " ";
        }
        jobInfoFile << "\n";
#endif
        jobInfoFile << "\n";
        jobInfoFile << "CPU time used since start of simulation (CPU-hours): "
                    << getCPUTime() / 3600.0;

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
            jobInfoFile << buildInfoGetModuleName(n) << ": "
                        << buildInfoGetModuleVal(n) << "\n";
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
            jobInfoFile << buildgitname << " git describe: " << buildgithash
                        << "\n";
        }

        jobInfoFile << "\n\n";

        // grid information
        jobInfoFile << PrettyLine;
        jobInfoFile << " Grid Information\n";
        jobInfoFile << PrettyLine;

        for (int i = 0; i <= finest_level; i++) {
            jobInfoFile << " level: " << i << "\n";
            jobInfoFile << "   number of boxes = " << grids[i].size() << "\n";
            jobInfoFile << "   maximum zones   = ";
            for (int n = 0; n < BL_SPACEDIM; n++) {
                jobInfoFile << geom[i].Domain().length(n) << " ";
            }
            jobInfoFile << "\n\n";
        }

        jobInfoFile << " Boundary conditions\n";
        Vector<int> lo_bc_out(BL_SPACEDIM), hi_bc_out(BL_SPACEDIM);
        ParmParse pp("maestro");
        pp.getarr("lo_bc", lo_bc_out, 0, BL_SPACEDIM);
        pp.getarr("hi_bc", hi_bc_out, 0, BL_SPACEDIM);

        // these names correspond to the integer flags setup in the
        // Castro_setup.cpp
        const char* names_bc[] = {"interior", "inflow",   "outflow",
                                  "symmetry", "slipwall", "noslipwall"};

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

        jobInfoFile << std::setw(6) << "index" << SkipSpace
                    << std::setw(mlen + 1) << "name" << SkipSpace
                    << std::setw(7) << "A" << SkipSpace << std::setw(7) << "Z"
                    << "\n";
        jobInfoFile << OtherLine;

        for (int i = 0; i < NumSpec; i++) {
            auto spec_name = short_spec_names_cxx[i];
            jobInfoFile << std::setw(6) << i << SkipSpace << std::setw(mlen + 1)
                        << std::setfill(' ') << spec_name << SkipSpace
                        << std::setw(7) << aion[i] << SkipSpace << std::setw(7)
                        << zion[i] << "\n";
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

void Maestro::WriteBuildInfo() {
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
        std::cout << buildInfoGetModuleName(n) << ": "
                  << buildInfoGetModuleVal(n) << "\n";
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

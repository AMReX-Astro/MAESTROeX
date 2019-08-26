
#include <Maestro.H>
#include <AMReX_VisMF.H>
#include <Maestro_F.H>

using namespace amrex;

namespace
{
const std::string level_prefix {"Level_"};
}

// compute S at cell-centers
void
Maestro::WriteCheckPoint (int step) {

		// chk00010            write a checkpoint file with this root directory
		// chk00010/Header     this contains information you need to save (e.g., finest_level, t_new, etc.) and also
		//                     the BoxArrays at each level
		// chk00010/Level_0/
		// chk00010/Level_1/
		// etc.                these subdirectories will hold the MultiFab data at each level of refinement

		// timer for profiling
		BL_PROFILE_VAR("Maestro::WriteCheckPoint()",WriteCheckPoint);

		// checkpoint file name, e.g., chk00010
		const std::string& checkpointname = amrex::Concatenate(check_base_name,step,7);

		amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

		const int nlevels = finest_level+1;

		// ---- prebuild a hierarchy of directories
		// ---- dirName is built first.  if dirName exists, it is renamed.  then build
		// ---- dirName/subDirPrefix_0 .. dirName/subDirPrefix_nlevels-1
		// ---- if callBarrier is true, call ParallelDescriptor::Barrier()
		// ---- after all directories are built
		// ---- ParallelDescriptor::IOProcessor() creates the directories
		amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, true);

		VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

		// write Header file
		if (ParallelDescriptor::IOProcessor()) {

				std::ofstream HeaderFile;
				HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
				std::string HeaderFileName(checkpointname + "/Header");
				HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
				                std::ofstream::trunc |
				                std::ofstream::binary);

				if( !HeaderFile.good()) {
						amrex::FileOpenFailed(HeaderFileName);
				}

				HeaderFile.precision(17);

				// write out title line
				HeaderFile << "Checkpoint file for MAESTRO\n";

				// write out the time step number
				HeaderFile << step << "\n";

				// write out finest_level
				HeaderFile << finest_level << "\n";

				// write out step
				HeaderFile << step << "\n";

				// write out dt
				HeaderFile << dt << "\n";

				// write out t_new
				HeaderFile << t_new << "\n";

				// write out rel_eps
				Real rel_eps;
				get_rel_eps(&rel_eps);
				HeaderFile << rel_eps << "\n";

				// write the BoxArray at each level
				for (int lev = 0; lev <= finest_level; ++lev) {
						boxArray(lev).writeOn(HeaderFile);
						HeaderFile << '\n';
				}

				{
						// store elapsed CPU time
						std::ofstream CPUFile;
						std::string FullPathCPUFile(checkpointname + "/CPUtime");
						CPUFile.open(FullPathCPUFile.c_str(), std::ios::out);

						CPUFile << std::setprecision(15) << getCPUTime();
						CPUFile.close();
				}


		}

		// write the MultiFab data to, e.g., chk00010/Level_0/
		for (int lev = 0; lev <= finest_level; ++lev) {
				VisMF::Write(snew[lev],
				             amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "snew"));
				VisMF::Write(unew[lev],
				             amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "unew"));
				VisMF::Write(gpi[lev],
				             amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "gpi"));
				VisMF::Write(dSdt[lev],
				             amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "dSdt"));
				VisMF::Write(S_cc_new[lev],
				             amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "S_cc_new"));
		}

		// write out the cell-centered base state
		if (ParallelDescriptor::IOProcessor()) {

				std::ofstream BaseCCFile;
				BaseCCFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
				std::string BaseCCFileName(checkpointname + "/BaseCC");
				BaseCCFile.open(BaseCCFileName.c_str(), std::ofstream::out   |
				                std::ofstream::trunc |
				                std::ofstream::binary);
				if( !BaseCCFile.good()) {
						amrex::FileOpenFailed(BaseCCFileName);
				}

				BaseCCFile.precision(17);

				for (int i=0; i<rho0_new.size(); ++i) {
						BaseCCFile << rho0_new[i] << " "
						           << p0_new[i] << " "
						           << gamma1bar_new[i] << " "
						           << rhoh0_new[i] << " "
						           << beta0_new[i] << " "
						           << psi[i] << " "
						           << tempbar[i] << " "
						           << etarho_cc[i] << " "
						           << tempbar_init[i] << " "
							   << p0_old[i] << " "
							   << beta0_nm1[i] << "\n";
				}
		}

		// write out the face-centered base state
		if (ParallelDescriptor::IOProcessor()) {

				std::ofstream BaseFCFile;
				BaseFCFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
				std::string BaseFCFileName(checkpointname + "/BaseFC");
				BaseFCFile.open(BaseFCFileName.c_str(), std::ofstream::out   |
				                std::ofstream::trunc |
				                std::ofstream::binary);
				if( !BaseFCFile.good()) {
						amrex::FileOpenFailed(BaseFCFileName);
				}

				BaseFCFile.precision(17);

				for (int i=0; i<w0.size(); ++i) {
						BaseFCFile << w0[i] << " "
						           << etarho_ec[i] << "\n";
				}
		}

		WriteJobInfo(checkpointname);
}

int
Maestro::ReadCheckPoint ()
{
		// timer for profiling
		BL_PROFILE_VAR("Maestro::ReadCheckPoint()",ReadCheckPoint);

		amrex::Print() << "Restart from checkpoint " << restart_file << "\n";

		VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

		std::string line, word;
		int step;

		// Header
		{
				std::string File(restart_file + "/Header");
				Vector<char> fileCharPtr;
				ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
				std::string fileCharPtrString(fileCharPtr.dataPtr());
				std::istringstream is(fileCharPtrString, std::istringstream::in);

				// read in title line
				std::getline(is, line);

				// read in time step number
				is >> start_step;
				GotoNextLine(is);
				++start_step;

				// read in finest_level
				is >> finest_level;
				GotoNextLine(is);

				// read in step
				is >> step;
				GotoNextLine(is);

				// read in dt
				is >> dt;
				GotoNextLine(is);

				// read in time
				is >> t_old;
				GotoNextLine(is);

				// read in rel_eps
				Real rel_eps;
				is >> rel_eps;
				GotoNextLine(is);
				set_rel_eps(&rel_eps);

				for (int lev = 0; lev <= finest_level; ++lev) {

						// read in level 'lev' BoxArray from Header
						BoxArray ba;
						ba.readFrom(is);
						GotoNextLine(is);

						// create a distribution mapping
						DistributionMapping dm { ba, ParallelDescriptor::NProcs() };

						// set BoxArray grids and DistributionMapping dmap in AMReX_AmrMesh.H class
						SetBoxArray(lev, ba);
						SetDistributionMap(lev, dm);

						// build MultiFab data
						sold              [lev].define(ba, dm,          Nscal, ng_s);
						snew              [lev].define(ba, dm,          Nscal, ng_s);
						uold              [lev].define(ba, dm, AMREX_SPACEDIM, ng_s);
						unew              [lev].define(ba, dm, AMREX_SPACEDIM, ng_s);
						S_cc_old          [lev].define(ba, dm,              1,    0);
						S_cc_new          [lev].define(ba, dm,              1,    0);
						gpi               [lev].define(ba, dm, AMREX_SPACEDIM,    0);
						dSdt              [lev].define(ba, dm,              1,    0);
						rhcc_for_nodalproj[lev].define(ba, dm,              1,    1);

						pi[lev].define(convert(ba,nodal_flag), dm, 1, 0); // nodal

						if (spherical == 1) {
								normal[lev].define(ba, dm, 1, 1);
								cell_cc_to_r[lev].define(ba, dm, 1, 0);
						}

						// build FluxRegister data
						if (lev > 0 && reflux_type == 2) {
								flux_reg_s[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, Nscal));
						}
				}
		}

		// read in the MultiFab data - put it in the "old" MultiFabs
		for (int lev = 0; lev <= finest_level; ++lev) {
				VisMF::Read(sold[lev],
				            amrex::MultiFabFileFullPrefix(lev, restart_file, "Level_", "snew"));
				VisMF::Read(uold[lev],
				            amrex::MultiFabFileFullPrefix(lev, restart_file, "Level_", "unew"));
				VisMF::Read(gpi[lev],
				            amrex::MultiFabFileFullPrefix(lev, restart_file, "Level_", "gpi"));
				VisMF::Read(dSdt[lev],
				            amrex::MultiFabFileFullPrefix(lev, restart_file, "Level_", "dSdt"));
				VisMF::Read(S_cc_old[lev],
				            amrex::MultiFabFileFullPrefix(lev, restart_file, "Level_", "S_cc_new"));
		}

		// get the elapsed CPU time to now;
		{
				// get ellapsed CPU time
				std::string CPUFile(restart_file + "/CPUtime");
				Vector<char> fileCharPtr;
				ParallelDescriptor::ReadAndBcastFile(CPUFile, fileCharPtr);
				std::string fileCharPtrString(fileCharPtr.dataPtr());
				std::istringstream is(fileCharPtrString, std::istringstream::in);

				// read in line
				std::getline(is, line);

				is >> previousCPUTimeUsed;

				std::cout << "read CPU time: " << previousCPUTimeUsed << "\n";

		}


		// BaseCC
		{
				std::string File(restart_file + "/BaseCC");
				Vector<char> fileCharPtr;
				ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
				std::string fileCharPtrString(fileCharPtr.dataPtr());
				std::istringstream is(fileCharPtrString, std::istringstream::in);

				// read in cell-centered base state
				for (int i=0; i<(max_radial_level+1)*nr_fine; ++i) {
						std::getline(is, line);
						std::istringstream lis(line);
						lis >> word;
						rho0_old[i] = std::stod(word);
						lis >> word;
						p0_old[i] = std::stod(word);
						lis >> word;
						gamma1bar_old[i] = std::stod(word);
						lis >> word;
						rhoh0_old[i] = std::stod(word);
						lis >> word;
						beta0_old[i] = std::stod(word);
						lis >> word;
						psi[i] = std::stod(word);
						lis >> word;
						tempbar[i] = std::stod(word);
						lis >> word;
						etarho_cc[i] = std::stod(word);
						lis >> word;
						tempbar_init[i] = std::stod(word);
						lis >> word;
						p0_nm1[i] = std::stod(word);
						lis >> word;
						beta0_nm1[i] = std::stod(word);
				}
		}


		if (do_smallscale) {
		    Average(sold,rho0_old,Rho);
		    compute_cutoff_coords(rho0_old.dataPtr());
		    std::fill(rho0_old.begin(),  rho0_old.end(),  0.);
		}

		// BaseFC
		{
				std::string File(restart_file + "/BaseFC");
				Vector<char> fileCharPtr;
				ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
				std::string fileCharPtrString(fileCharPtr.dataPtr());
				std::istringstream is(fileCharPtrString, std::istringstream::in);


				// read in face-centered base state
				for (int i=0; i<(max_radial_level+1)*nr_fine+1; ++i) {
						std::getline(is, line);
						std::istringstream lis(line);
						lis >> word;
						w0[i] = std::stod(word);
						lis >> word;
						etarho_ec[i] = std::stod(word);
				}
		}

		return step;
}


// utility to skip to next line in Header
void
Maestro::GotoNextLine (std::istream& is)
{
		// timer for profiling
		BL_PROFILE_VAR("Maestro::GotoNextLine()",GotoNextLine);

		constexpr std::streamsize bl_ignore_max { 100000 };
		is.ignore(bl_ignore_max, '\n');
}

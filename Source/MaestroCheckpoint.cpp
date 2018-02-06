
#include <Maestro.H>
#include <AMReX_VisMF.H>

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

    // checkpoint file name, e.g., chk00010
    const std::string& checkpointname = amrex::Concatenate(check_base_name,step);

    amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

    const int nlevels = finest_level+1;

    // ---- prebuild a hierarchy of directories
    // ---- dirName is built first.  if dirName exists, it is renamed.  then build
    // ---- dirName/subDirPrefix_0 .. dirName/subDirPrefix_nlevels-1
    // ---- if callBarrier is true, call ParallelDescriptor::Barrier()
    // ---- after all directories are built
    // ---- ParallelDescriptor::IOProcessor() creates the directories
    amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, true);

    // write Header file
   if (ParallelDescriptor::IOProcessor()) {

       std::string HeaderFileName(checkpointname + "/Header");
       std::ofstream HeaderFile(HeaderFileName.c_str(), std::ofstream::out   |
				                        std::ofstream::trunc |
				                        std::ofstream::binary);
       if( ! HeaderFile.good()) {
           amrex::FileOpenFailed(HeaderFileName);
       }

       HeaderFile.precision(17);

       VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
       HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

       // write out title line
       HeaderFile << "Checkpoint file for MAESTRO\n";

       // write out finest_level
       HeaderFile << finest_level << "\n";

       // write out step
       HeaderFile << step << "\n";

       // write out dt
       HeaderFile << dt << "\n";

       // write out t_new
       HeaderFile << t_new << "\n";

       // write the BoxArray at each level
       for (int lev = 0; lev <= finest_level; ++lev) {
           boxArray(lev).writeOn(HeaderFile);
           HeaderFile << '\n';
       }
   }

   // write the MultiFab data to, e.g., chk00010/Level_0/
   for (int lev = 0; lev <= finest_level; ++lev) {
       VisMF::Write(snew[lev],
                    amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "snew"));
   }

}

int
Maestro::ReadCheckPoint ()
{

    amrex::Print() << "Restart from checkpoint " << restart_file << "\n";

    // Header
    std::string File(restart_file + "/Header");

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line, word;

    // read in title line
    std::getline(is, line);

    // read in finest_level
    is >> finest_level;
    GotoNextLine(is);

    // read in step
    int step;
    is >> step;
    GotoNextLine(is);

    // read in dt
    is >> dt;
    GotoNextLine(is);

    // read in t_new
    is >> t_new;
    GotoNextLine(is);

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

        // build MultiFab and FluxRegister data
        sold    [lev].define(ba, dm,          Nscal, 3);
        snew    [lev].define(ba, dm,          Nscal, 3);
        uold    [lev].define(ba, dm, AMREX_SPACEDIM, 3);
        unew    [lev].define(ba, dm, AMREX_SPACEDIM, 3);
        S_cc_old[lev].define(ba, dm,              1, 0);
        S_cc_new[lev].define(ba, dm,              1, 0);
        gpi     [lev].define(ba, dm, AMREX_SPACEDIM, 0);
        dSdt    [lev].define(ba, dm,              1, 0);
        pi      [lev].define(convert(ba,nodal_flag), dm, 1, 0); // nodal

        rhcc_for_nodalproj[lev].define(ba, dm, 1, 1);

        if (spherical == 1) {
            normal[lev].define(ba, dm, 1, 1);
        }

        if (lev > 0 && do_reflux) {
            flux_reg_s[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, Nscal));
            flux_reg_u[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, AMREX_SPACEDIM));
        }
    }

    // read in the MultiFab data
    for (int lev = 0; lev <= finest_level; ++lev) {
        VisMF::Read(snew[lev],
                    amrex::MultiFabFileFullPrefix(lev, restart_file, "Level_", "phi"));
    }

    return step;
}


// utility to skip to next line in Header
void
Maestro::GotoNextLine (std::istream& is)
{
    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore(bl_ignore_max, '\n');
}

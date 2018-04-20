
#include <Maestro.H>

using namespace amrex;


// print out the contents of a Vector of MultiFabs
void
Maestro::PrintMF (Vector<MultiFab>& MF)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::PrintMF()",PrintMF);

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& MF_mf = MF[lev];

        const BoxArray& ba = MF_mf.boxArray();
        const DistributionMapping& dm = MF_mf.DistributionMap();
        const int myProc = ParallelDescriptor::MyProc();

        for (int i=0; i<ba.size(); ++i) {
            if (dm[i] == myProc) {

                // we want all processors to write, not just the IOProcessor
                std::cout << "Grid #" << i << std::endl;
                std::cout << "Processor #" << myProc << std::endl;

                const Box& validBox = ba[i];

                // call fortran subroutine
                // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
                // lo/hi coordinates (including ghost cells), and/or the # of components
                // We will also pass "validBox", which specifies the "valid" region.
                print_mf(&lev, ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                         BL_TO_FORTRAN_FAB(MF_mf[i]));

            }
            // add this barrier so only one grid gets printed out at a time
            ParallelDescriptor::Barrier();
        }
    }
}

void
Maestro::PrintEdge (Vector<std::array< MultiFab, AMREX_SPACEDIM > >& EDGE,
                    int dir)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::PrintEdge()",PrintEdge);

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& EDGE_mf = EDGE[lev][dir];

        const BoxArray& ba = EDGE_mf.boxArray();
        const DistributionMapping& dm = EDGE_mf.DistributionMap();
        const int myProc = ParallelDescriptor::MyProc();

        for (int i=0; i<ba.size(); ++i) {
            if (dm[i] == myProc) {

                // we want all processors to write, not just the IOProcessor
                std::cout << "Grid #" << i << std::endl;
                std::cout << "Processor #" << myProc << std::endl;

                // EDGE BASED
                const Box& validBox = ba[i];

                // call fortran subroutine
                // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
                // lo/hi coordinates (including ghost cells), and/or the # of components
                // We will also pass "validBox", which specifies the "valid" region.
                print_edge(&lev, ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                           BL_TO_FORTRAN_FAB(EDGE_mf[i]));

            }
            // add this barrier so only one grid gets printed out at a time
            ParallelDescriptor::Barrier();
        }
    }
}

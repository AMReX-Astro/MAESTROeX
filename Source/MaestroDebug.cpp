
#include <Maestro.H>

using namespace amrex;


// print out the contents of a Vector of MultiFabs
void
Maestro::PrintMF (const Vector<MultiFab>& MF)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::PrintMF()",PrintMF);

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& MF_mf = MF[lev];

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
Maestro::PrintEdge (const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& EDGE,
                    int dir)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::PrintEdge()",PrintEdge);

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& EDGE_mf = EDGE[lev][dir];

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

// utility to write out a multilevel multifab to a plotfile
void Maestro::WriteMF (const Vector<MultiFab>& mf,
                       std::string name)
{
    int nComp = mf[0].nComp();

    Vector<std::string> varnames;
    varnames.resize(nComp);
    for (int i=0; i<nComp; ++i) {
        varnames[i] = "X";
    }

    // temporary MultiFab to hold plotfile data
    Vector<MultiFab*> plot_mf_data(finest_level+1);

    // build temporary MultiFab to hold plotfile data
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i] = new MultiFab((mf[i]).boxArray(),(mf[i]).DistributionMap(),nComp,0);
    }
    
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i]->copy((mf[i]),0,0,nComp);
    }
                    
    // MultiFab to hold plotfile data
    Vector<const MultiFab*> plot_mf;
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf.push_back(plot_mf_data[i]);
    }    
    
    Vector<int> step_array;
    step_array.resize(maxLevel()+1, 0);

    WriteMultiLevelPlotfile(name, finest_level+1, plot_mf, varnames,
                            Geom(), 0., step_array, refRatio());
    

}

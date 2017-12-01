
#include <Maestro.H>

using namespace amrex;

// get plotfile name
std::string
Maestro::PlotFileName (int lev) const
{
    return Concatenate(plot_base_name, lev, 5);
}

// put together a vector of multifabs for writing
Vector<const MultiFab*>
Maestro::PlotFileMF (Vector<MultiFab>& p0_cart,
                     Vector<MultiFab>& rho0_cart) const
{

    // velocities (AMREX_SPACEDIM)
    // rho, rhoh, rhoX, Temp, Pi (Nscal)
    // rho0, p0 (2)
    int nPlot = AMREX_SPACEDIM + Nscal + 2;

    // temporary MultiFab to hold plotfile data
    Vector<const MultiFab*> plot_mf;

    // plotfile data gets copied into here
    Vector<MultiFab*> plot_mf_data(finest_level+1);

    for (int i = 0; i <= finest_level; ++i)
    {
        plot_mf_data[i] = new MultiFab((snew[i]).boxArray(),(snew[i]).DistributionMap(),nPlot,0);

        // copy velocity and scalars into plot_mf_data[i]
        int dest_comp = 0;
        plot_mf_data[i]->copy((unew[i]),     0,dest_comp,AMREX_SPACEDIM);
        dest_comp += AMREX_SPACEDIM;
        plot_mf_data[i]->copy((snew[i]),     0,dest_comp,Nscal);
        dest_comp += Nscal;
        plot_mf_data[i]->copy((rho0_cart[i]),0,dest_comp,1);
        dest_comp += 1;
        plot_mf_data[i]->copy((p0_cart[i]),  0,dest_comp,1);
        dest_comp += 1;

        // add plot_mf_data[i] to plot_mf
        plot_mf.push_back(plot_mf_data[i]);
    }

    return plot_mf;

}

// set plotfile variable names
Vector<std::string>
Maestro::PlotFileVarNames () const
{

    // velocities (AMREX_SPACEDIM)
    // rho, rhoh, rhoX, Temp, Pi (Nscal)
    // rho0, p0 (2)
    int nPlot = Nscal + AMREX_SPACEDIM + 2;
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

    for (int i = 0; i < NumSpec; i++)
    {
        int len = 20;
        Array<int> int_spec_names(len);
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

    names[cnt++] = "Temp";
    names[cnt++] = "Pi";

    names[cnt++] = "rho0";
    names[cnt++] = "p0";

    return names;
}

// write plotfile to disk
void
Maestro::WritePlotFile (int step)
{
    // wallclock time
    const Real strt_total = ParallelDescriptor::second();

    const std::string& plotfilename = PlotFileName(step);

    // convert p0 to multi-D MultiFab
    Vector<MultiFab> p0_cart(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        p0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
    }
    Put1dArrayOnCart(p0_new,p0_cart,bcs_f,0,0);

    // convert rho0 to multi-D MultiFab
    Vector<MultiFab> rho0_cart(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        rho0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
    }
    Put1dArrayOnCart(rho0_new,rho0_cart,bcs_f,0,0);

    const auto& mf = PlotFileMF(p0_cart,rho0_cart);
    const auto& varnames = PlotFileVarNames();

    // WriteMultiLevelPlotfile expects an array of step numbers
    Vector<int> step_array;
    step_array.resize(maxLevel()+1, step);

    WriteMultiLevelPlotfile(plotfilename, finest_level+1, mf, varnames,
                            Geom(), t_new, step_array, refRatio());

    // wallclock time
    Real end_total = ParallelDescriptor::second() - strt_total;
	
    // print wallclock time
    ParallelDescriptor::ReduceRealMax(end_total ,ParallelDescriptor::IOProcessorNumber());
    if (verbose > 0) {
        Print() << "Time to write plotfile: " << end_total << '\n';
    }

}

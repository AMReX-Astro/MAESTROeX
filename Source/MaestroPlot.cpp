
#include <Maestro.H>

// get plotfile name
std::string
Maestro::PlotFileName (int lev) const
{
    return amrex::Concatenate(plot_file, lev, 5);
}

// put together an array of multifabs for writing
Array<const MultiFab*>
Maestro::PlotFileMF () const
{
    Array<const MultiFab*> r;
    for (int i = 0; i <= finest_level; ++i) {
        r.push_back(snew[i].get());
    }
    return r;
}

// set plotfile variable names
Array<std::string>
Maestro::PlotFileVarNames () const
{
    return {"phi1", "phi2", "phi3", "phi4", "phi5", "phi6"};
}

// write plotfile to disk
void
Maestro::WritePlotFile (int step) const
{
    const std::string& plotfilename = PlotFileName(step);
    const auto& mf = PlotFileMF();
    const auto& varnames = PlotFileVarNames();

    // WriteMultiLevelPlotfile expects an array of step numbers
    amrex::Array<int> step_array;
    step_array.resize(maxLevel()+1, step);

    amrex::WriteMultiLevelPlotfile(plotfilename, finest_level+1, mf, varnames,
                                   Geom(), t_new, step_array, refRatio());
}

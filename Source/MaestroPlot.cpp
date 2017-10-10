
#include <Maestro.H>

using namespace amrex;

// get plotfile name
std::string
Maestro::PlotFileName (int lev) const
{
    return Concatenate(plot_file, lev, 5);
}

// put together a vector of multifabs for writing
Vector<const MultiFab*>
Maestro::PlotFileMF () const
{
    Vector<const MultiFab*> r;
    for (int i = 0; i <= finest_level; ++i) {
        r.push_back(snew[i].get());
    }
    return r;
}

// set plotfile variable names
Vector<std::string>
Maestro::PlotFileVarNames () const
{
    Vector<std::string> names(NSCAL);

    int cnt = 0;

    names[cnt++] = "rho";
    names[cnt++] = "rhoh";

    for (int i = 0; i < NumSpec; i++)
    {
        int len = 20;
        Array<int> int_spec_names(len);
        //
        // This call return the actual length of each string in "len"
        //
        ca_get_spec_names(int_spec_names.dataPtr(),&i,&len);
        char* spec_name = new char[len+1];
        for (int j = 0; j < len; j++)
            spec_name[j] = int_spec_names[j];
        spec_name[len] = '\0';
        std::string spec_string = "X(";
        spec_string += spec_name;
        spec_string += ')';

        names[cnt++] = spec_string;
    }

    names[cnt++] = "Pi";
    names[cnt++] = "Temp";

/*
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        std::string x = "vel";
        x += (120+i);
        Print() << x << endl;
    }
*/
    return names;
}

// write plotfile to disk
void
Maestro::WritePlotFile (int step) const
{
    const std::string& plotfilename = PlotFileName(step);
    const auto& mf = PlotFileMF();
    const auto& varnames = PlotFileVarNames();

    // WriteMultiLevelPlotfile expects an array of step numbers
    Vector<int> step_array;
    step_array.resize(maxLevel()+1, step);

    WriteMultiLevelPlotfile(plotfilename, finest_level+1, mf, varnames,
                            Geom(), t_new, step_array, refRatio());
}

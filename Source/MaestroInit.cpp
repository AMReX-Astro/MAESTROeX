

#include <Maestro.H>

// initializes multilevel data
void
Maestro::InitData ()
{
    const Real time = 0.0;
    InitFromScratch(time);
    AverageDown();

    if (plot_int > 0) {
        WritePlotFile(0);
    }
}


#include <Maestro.H>
#include <extern_parameters.H>

using namespace amrex;

// advance solution to final time
void Maestro::Evolve() {
    Vector<MultiFab> rho_omegadot(finest_level + 1);
    Vector<MultiFab> rho_Hnuc(finest_level + 1);
    Vector<MultiFab> rho_Hext(finest_level + 1);

    for (int lev = 0; lev <= finest_level; ++lev) {
        // cell-centered MultiFabs
        rho_omegadot[lev].define(grids[lev], dmap[lev], NumSpec, 0);
        rho_Hnuc[lev].define(grids[lev], dmap[lev], 1, 0);
        rho_Hext[lev].define(grids[lev], dmap[lev], 1, 0);
        rho_Hext[lev].setVal(0.);
    }

    auto dbo = do_burning;
    auto dho = do_heating;

    // This problem uses a custom WritePlotFile, but as it's a member function of
    // the Maestro class, it has to use the same prototype as the original.
    // We shall therefore create a dummy variable to fill up all the variables
    // passed into the function that won't be used.
    BaseState<Real> dummy;

    // Model 1: No burning, no heating
    Print() << "\nModel 1: No burning, no heating\n";
    do_burning = false;
    do_heating = false;
    React(sold, snew, rho_Hext, rho_omegadot, rho_Hnuc, p0_old, dt, t_old);
    WritePlotFile(-1, t_new, dt, dummy, dummy, dummy, dummy, rho_omegadot,
                  rho_Hnuc, rho_Hext);

    // Model 2: Burning without heating
    Print() << "\nModel 2: Burning without heating\n";
    do_burning = true;
    do_heating = false;
    React(sold, snew, rho_Hext, rho_omegadot, rho_Hnuc, p0_old, dt, t_old);
    WritePlotFile(-2, t_new, dt, dummy, dummy, dummy, dummy, rho_omegadot,
                  rho_Hnuc, rho_Hext);

    // Model 3: Heating without burning
    Print() << "\nModel 3: Heating without burning\n";
    do_burning = false;
    do_heating = true;
    React(sold, snew, rho_Hext, rho_omegadot, rho_Hnuc, p0_old, dt, t_old);
    WritePlotFile(-3, t_new, dt, dummy, dummy, dummy, dummy, rho_omegadot,
                  rho_Hnuc, rho_Hext);

    // Model 4: Burning and heating
    Print() << "\nModel 4: Burning and heating\n";
    do_burning = true;
    do_heating = true;
    React(sold, snew, rho_Hext, rho_omegadot, rho_Hnuc, p0_old, dt, t_old);
    WritePlotFile(-4, t_new, dt, dummy, dummy, dummy, dummy, rho_omegadot,
                  rho_Hnuc, rho_Hext);

    // Explore ten orders of magnitude of the time domain using user inputs.
    do_burning = dbo;
    do_heating = dho;

    for (auto i = 0; i < react_its; ++i) {
        React(sold, snew, rho_Hext, rho_omegadot, rho_Hnuc, p0_old, dt, t_old);
        WritePlotFile(i, t_new, dt, dummy, dummy, dummy, dummy, rho_omegadot,
                      rho_Hnuc, rho_Hext);
    }
}

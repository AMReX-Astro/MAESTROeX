

#include <Maestro.H>

using namespace amrex;

// define variable mappings (Rho, RhoH, ..., NUM_STATE, etc.)
void Maestro::VariableSetup ()
{

    int cnt = 0;
    Rho = cnt++;
    RhoH = cnt++;

    FirstSpec = cnt;
    ca_get_num_spec(&NumSpec);
    cnt += NumSpec;

    Temp = cnt++;
    Pi = cnt++;

    NUM_STATE = cnt;

}

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



// Make a new level from scratch using provided BoxArray and DistributionMapping.
// Only used during initialization.
// overrides the pure virtual function in AmrCore
void Maestro::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
				       const DistributionMapping& dm)
{
    const int nghost = 0;

    snew[lev].reset(new MultiFab(ba, dm, NUM_STATE, nghost));
    sold[lev].reset(new MultiFab(ba, dm, NUM_STATE, nghost));

    t_new = time;
    t_old = time - 1.e200;

    if (lev > 0 && do_reflux) {
        flux_reg_s[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, NUM_STATE));
        flux_reg_u[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, AMREX_SPACEDIM));
    }

    const Real* dx = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();
    Real cur_time = t_new;

    MultiFab& state = *snew[lev];

    for (MFIter mfi(state); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        const int* lo  = box.loVect();
        const int* hi  = box.hiVect();

        initdata(lev, cur_time, ARLIM_3D(lo), ARLIM_3D(hi),
                 BL_TO_FORTRAN_3D(state[mfi]), ZFILL(dx),
                 ZFILL(prob_lo), NUM_STATE);
    }
}

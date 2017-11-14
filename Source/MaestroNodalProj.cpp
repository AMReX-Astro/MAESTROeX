
#include <Maestro.H>

using namespace amrex;


void
Maestro::NodalProj (Vector<MultiFab>& phi,
                    Vector<MultiFab>& vel,
                    Vector<MultiFab>& rhcc,
                    Vector<MultiFab>& rhnd,
                    Vector<MultiFab>& beta0,
                    int* mg_bcs,
                    Real rel_tol,
                    Real abs_tol)

{
    const Vector<Geometry>& mg_geom = Geom();
    const Vector<BoxArray>& mg_ba = boxArray();
    const Vector<DistributionMapping>& mg_dm = DistributionMap();

    const bool nodal = true;
    const int hg_stencil = ND_CROSS_STENCIL;
    const bool have_rhcc = false;
    const int nc = 0;
    const int ncomp = 1;
    const int verbose = 0;

    MGT_Solver mgt_solver(mg_geom, mg_bcs, mg_ba, mg_dm, nodal, hg_stencil, have_rhcc,
                          nc, ncomp, verbose);

    mgt_solver.set_nodal_coefficients(GetVecOfPtrs(beta0));

    mgt_solver.nodal_project(GetVecOfPtrs(phi),
                             GetVecOfPtrs(vel),
                             GetVecOfPtrs(rhcc),
                             GetVecOfPtrs(rhnd),
                             rel_tol,
                             abs_tol,
                             &lo_inflow[0],
                             &hi_inflow[0]);
}

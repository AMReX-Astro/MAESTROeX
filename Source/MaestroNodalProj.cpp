
#include <Maestro.H>

using namespace amrex;


void
Maestro::NodalProj (const Vector<std::unique_ptr<MultiFab> >& phi,
                    const Vector<std::unique_ptr<MultiFab> >& vel,
                    const Vector<std::unique_ptr<MultiFab> >& rhcc,
                    const Vector<std::unique_ptr<MultiFab> >& rhnd,
                    const Vector<std::unique_ptr<MultiFab> >& beta0,
                    int* mg_bcs,
                    const Real& rel_tol,
                    const Real& abs_tol)

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

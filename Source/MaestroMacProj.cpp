
#include <Maestro.H>

using namespace amrex;

void
Maestro::MacProj (Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
                  Vector<MultiFab>& macphi,
                  const Vector<MultiFab>& macrhs,
                  const Vector<Real>& beta0)
{
    // umac enters with face-centered, time-centered Utilde^* and should leave with Utilde
    // macphi is the solution to the elliptic solve and 
    //   enters as either zero, or the solution to the predictor MAC projection
    // macrhs enters as beta0*(S-Sbar)
    // beta0 is a 1d cell-centered array    

    // this will hold solver RHS = macrhs - div(beta0*umac)
    Vector<MultiFab> solverrhs(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        solverrhs[lev].define(grids[lev], dmap[lev], 1, 0);
    }

    // we also need beta0 at edges
    // allocate AND compute it here
    Vector<Real> beta0_edge( (max_radial_level+1)*(nr_fine+1) );
    cell_to_edge(beta0.dataPtr(),beta0_edge.dataPtr());

    // convert Utilde^* to beta0*Utilde^*
    MultUmacByBeta0(umac,beta0,beta0_edge,1);

    // compute the RHS for the solve, RHS = macrhs - div(beta0*umac)
    ComputeMACSolverRHS(solverrhs,macrhs,umac);

    // set face-centered B coefficients to beta/rho
    //
    //

    // set cell-centered A coefficient to zero
    //
    //

    // solve -div B grad phi = RHS
    //
    //

    // update velocity, beta0 * Utilde = beta0 * Utilde^* - B grad phi
    //
    //

    // convert beta0*Utilde to Utilde
    MultUmacByBeta0(umac,beta0,beta0_edge,0);

    // fill peroidic ghost cells
    for (int lev=0; lev<=finest_level; ++lev) {
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            umac[lev][d].FillBoundary(geom[lev].periodicity());
        }
    }

    // fill ghost cells behind physical boundaries
    FillUmacGhost(umac);

}

// multiply (or divide) umac by beta0
void Maestro::MultUmacByBeta0 (Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
                               const Vector<Real>& beta0,
                               const Vector<Real>& beta0_edge,
                               const int& mult_or_div)
{





}

// compute the RHS for the solve, RHS = macrhs - div(beta0*umac)
void Maestro::ComputeMACSolverRHS (Vector<MultiFab>& solverrhs,
                                   const Vector<MultiFab>& macrhs,
                                   const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac)
{



}

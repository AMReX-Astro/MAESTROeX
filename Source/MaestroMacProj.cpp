
#include <Maestro.H>

using namespace amrex;

// umac enters with face-centered, time-centered Utilde^* and should leave with Utilde
// macphi is the solution to the elliptic solve and 
//   enters as either zero, or the solution to the predictor MAC projection
// macrhs enters as beta0*(S-Sbar)
// beta0 is a 1d cell-centered array    
void
Maestro::MacProj (Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
                  Vector<MultiFab>& macphi,
                  const Vector<MultiFab>& macrhs,
                  const Vector<Real>& beta0,
                  const bool& is_predictor)
{
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
    int mult_or_div = 1;
    MultFacesByBeta0(umac,beta0,beta0_edge,mult_or_div);

    // compute the RHS for the solve, RHS = macrhs - div(beta0*umac)
    ComputeMACSolverRHS(solverrhs,macrhs,umac);

    // create a MultiFab filled with rho and 1 ghost cell.
    // if this is the predictor mac projection, use rho^n
    // if this is the corrector mac projection, use (1/2)(rho^n + rho^{n+1,*})
    Vector<MultiFab> rho(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        rho[lev].define(grids[lev], dmap[lev], 1, 1);
    }
    Real rho_time = is_predictor ? t_old : 0.5*(t_old+t_new);
    for (int lev=0; lev<=finest_level; ++lev) {
        FillPatch(lev, rho_time, rho[lev], sold, snew, Rho, Rho, 1, bcs_s);
    }

    // set face-centered B coefficients to 1/rho by averaging neighboring
    // cell-centered values of rho
//    MakeMacCoeffs(face_bcoeff,rho);
    //

    // multiply face-centered B coefficients by beta0 so they contain beta0/rho
    // MultFacesByBeta0();
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
    mult_or_div = 0;
    MultFacesByBeta0(umac,beta0,beta0_edge,mult_or_div);

    // fill peroidic ghost cells
    for (int lev=0; lev<=finest_level; ++lev) {
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            umac[lev][d].FillBoundary(geom[lev].periodicity());
        }
    }

    // fill ghost cells behind physical boundaries
    FillUmacGhost(umac);

}

// multiply (or divide) face-data by beta0
void Maestro::MultFacesByBeta0 (Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
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

void Maestro::MakeMacCoeffs(Vector<std::array< MultiFab, AMREX_SPACEDIM > >& bcoeff,
                            const Vector<MultiFab>& rho)
{

}

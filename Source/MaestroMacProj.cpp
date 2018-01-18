
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

    // coefficients for solver
    Vector<MultiFab> acoef(finest_level+1);
    Vector<MultiFab> bcoef(finest_level+1);
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > face_bcoef(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        acoef[lev].define(grids[lev], dmap[lev], 1, 0);
        bcoef[lev].define(grids[lev], dmap[lev], 1, 1);
        AMREX_D_TERM(face_bcoef[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 0);,
                     face_bcoef[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 0);,
                     face_bcoef[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 0););

    }

    // set cell-centered A coefficient to zero
    for (int lev=0; lev<=finest_level; ++lev) {
        acoef[lev].setVal(0.);
    }

    // set face-centered B coefficients to 1/rho
    // first set the cell-centered B coefficients to 1/rho
    for (int lev=0; lev<=finest_level; ++lev) {
        bcoef[lev].setVal(1.);
        bcoef[lev].divide(rho[lev],0,0,1);
    }

    // average bcoef to faces
    for (int lev=0; lev<=finest_level; ++lev) {
        amrex::average_cellcenter_to_face({AMREX_D_DECL(&face_bcoef[lev][0],
                                                        &face_bcoef[lev][1],
                                                        &face_bcoef[lev][2])},
                                            bcoef[lev], geom[lev]);
    }

    // multiply face-centered B coefficients by beta0 so they contain beta0/rho
    mult_or_div = 1;
    MultFacesByBeta0(face_bcoef,beta0,beta0_edge,mult_or_div);

    // set up the solver
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

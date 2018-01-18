
#include <Maestro.H>

using namespace amrex;

void
Maestro::MacProj (Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
                  Vector<MultiFab>& macphi,
                  const Vector<MultiFab>& macrhs,
                  const Vector<Real>& beta0)
{

    Vector<Real> beta0_edge( (max_radial_level+1)*(nr_fine+1) );

    // compute beta0 on edges
    cell_to_edge(beta0.dataPtr(),beta0_edge.dataPtr());

    // convert umac to beta0*umac

    // compute div(beta0*umac)

    // set beta coefficients to beta/rho

    // solve

    // update velocity

    // convert beta0*umac to umac

    // fill ghost cells

}

// this contains interfaces to fortran routines that deal with the base state

#include <Maestro.H>

using namespace amrex;

void
Maestro::MakeGravCC(Vector<Real>& grav,
                    const Vector<Real>& rho0)
{

    make_grav_cell(grav.dataPtr(),rho0.dataPtr(),r_cc_loc.dataPtr(),r_edge_loc.dataPtr());

}

void
Maestro::MakeDivCoeff(Vector<Real>& div_coeff,
                      const Vector<Real>& rho0,
                      const Vector<Real>& p0,
                      const Vector<Real>& gamma1bar_in,
                      const Vector<Real>& grav_cell_in)
{
    make_div_coeff(div_coeff.dataPtr(),rho0.dataPtr(),p0.dataPtr(),
                   gamma1bar_in.dataPtr(),grav_cell_in.dataPtr());
}

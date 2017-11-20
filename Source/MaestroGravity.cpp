
#include <Maestro.H>

using namespace amrex;

void
Maestro::MakeGravCC(Vector<Real>& grav,
                    const Vector<Real>& rho0)
{

    make_grav_cell(grav.dataPtr(),rho0.dataPtr(),r_cc_loc.dataPtr(),r_edge_loc.dataPtr());

}

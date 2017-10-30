
#include <Maestro.H>
#include <AMReX_FMultiGrid.H>

using namespace amrex;


void
Maestro::NodalProj ()
{

    const Vector<Geometry>& mg_geom = Geom();

    // order is xlo, xhi, ylo, yhi, zlo, zhi
    // MGT_BC_INT (Interior)
    // MGT_BC_DIR (Dirichlet)
    // MGT_BC_NEU (Neumann)
    int mg_bc[2*BL_SPACEDIM];
    for ( int i = 0; i < BL_SPACEDIM; ++i ) {
        if ( mg_geom[0].isPeriodic(i) ) {
            mg_bc[i*2 + 0] = MGT_BC_INT;
            mg_bc[i*2 + 1] = MGT_BC_INT;
        }
        else {
//            mg_bc[i*2 + 0] = phys_bc->lo(i)==Outflow? MGT_BC_DIR : MGT_BC_NEU;
//            mg_bc[i*2 + 1] = phys_bc->hi(i)==Outflow? MGT_BC_DIR : MGT_BC_NEU;
        }
    }


}

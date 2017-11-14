
#include <Maestro.H>

using namespace amrex;

// Given a multifab of data (phi), average down to a base state quantity, phibar.
// If we are in plane-parallel, the averaging is at constant height.  
// If we are spherical, then the averaging is done at constant radius.  

void Maestro::Average (Vector<MultiFab>& phi,
                       Vector<Real> phibar)
{

//    Vector<Real> phisum     ((max_level+1)*nr_fine,0.0);
//    Vector<Real> phisum_proc((max_level+1)*nr_fine,0.0);












}

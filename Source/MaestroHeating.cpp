
#include <Maestro.H>

using namespace amrex;

// compute heating term, rho_Hext
void
Maestro::ComputeHeating (Vector<MultiFab>& rho_Hext) 
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::ComputeHeating()",ComputeHeating);
   
    // FIXME
    for (int lev=0; lev<=finest_level; ++lev) {
        rho_Hext[lev].setVal(0.);
    }        

}

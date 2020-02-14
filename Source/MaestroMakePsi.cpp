#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void 
Maestro::MakePsiPlanar(RealVector& psi) 
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakePsiPlanar()", MakePsiPlanar);

    const int max_lev = max_radial_level+1;
    get_numdisjointchunks(numdisjointchunks.dataPtr());
    get_r_start_coord(r_start_coord.dataPtr());
    get_r_end_coord(r_end_coord.dataPtr());
    get_finest_radial_level(&finest_radial_level);

    RealVector grav_edge((finest_radial_level+1)*(nr_fine+1));
    RealVector p0old((finest_radial_level+1)*(nr_fine+1));

    
}

void 
Maestro::MakePsiSphr(RealVector& psi) 
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakePsiSphr()", MakePsiSphr);

    const int max_lev = max_radial_level+1;
    get_numdisjointchunks(numdisjointchunks.dataPtr());
    get_r_start_coord(r_start_coord.dataPtr());
    get_r_end_coord(r_end_coord.dataPtr());
    get_finest_radial_level(&finest_radial_level);

    RealVector grav_edge((finest_radial_level+1)*(nr_fine+1));
    RealVector p0old((finest_radial_level+1)*(nr_fine+1));

}

void 
Maestro::MakePsiIrreg(RealVector& psi) 
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakePsiIrreg()", MakePsiIrreg);

    const int max_lev = max_radial_level+1;
    get_numdisjointchunks(numdisjointchunks.dataPtr());
    get_r_start_coord(r_start_coord.dataPtr());
    get_r_end_coord(r_end_coord.dataPtr());
    get_finest_radial_level(&finest_radial_level);

    RealVector grav_edge((finest_radial_level+1)*(nr_fine+1));
    RealVector p0old((finest_radial_level+1)*(nr_fine+1));

}
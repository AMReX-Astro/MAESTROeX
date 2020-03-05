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

    std::fill(psi.begin(), psi.end(), 0.0);

    for (auto n = 0; n <= finest_radial_level; ++n) {
        for (auto i = 1; i <= numdisjointchunks[n]; ++i){
            for (auto r = r_start_coord[n+max_lev*i]; 
                 r<= r_end_coord[n+max_lev*i]; ++r) {
                if (r < base_cutoff_density_coord[n]) {
                    psi[n+max_lev*r] = etarho_cc[n+max_lev*r] * fabs(grav_const);
                }
            }
        }
    }

    RestrictBase(psi, true);
    FillGhostBase(psi, true);
}

void 
Maestro::MakePsiSphr(RealVector& psi, const RealVector& w0_in, 
                     const RealVector& gamma1bar, 
                     const RealVector& p0_avg,
                     const RealVector& Sbar_in) 
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakePsiSphr()", MakePsiSphr);

    const int max_lev = max_radial_level+1;

    std::fill(psi.begin(), psi.end(), 0.0);

    Real dr0 = dr[0];

    const Real * AMREX_RESTRICT r_cc_loc_p = r_cc_loc.dataPtr();
    const Real * AMREX_RESTRICT r_edge_loc_p = r_edge_loc.dataPtr();
    const Real * AMREX_RESTRICT w0_p = w0_in.dataPtr();
    const Real * AMREX_RESTRICT gamma1bar_p = gamma1bar.dataPtr();
    const Real * AMREX_RESTRICT p0_avg_p = p0_avg.dataPtr();
    const Real * AMREX_RESTRICT Sbar_p = Sbar_in.dataPtr();
    Real * AMREX_RESTRICT psi_p = psi.dataPtr();

    const auto npts = base_cutoff_density_coord[0];
    AMREX_PARALLEL_FOR_1D(npts, r, {
        Real div_w0_sph = 1.0 / (r_cc_loc_p[max_lev*r]*r_cc_loc_p[max_lev*r]) * 
            (r_edge_loc_p[max_lev*(r+1)]*r_edge_loc_p[max_lev*(r+1)] *
             w0_p[max_lev*(r+1)] - 
             r_edge_loc_p[max_lev*r]*r_edge_loc_p[max_lev*r] * 
             w0_p[max_lev*r]) / dr0;

        psi_p[max_lev*r] = gamma1bar_p[max_lev*r] * p0_avg_p[max_lev*r] * 
            (Sbar_p[max_lev*r] - div_w0_sph);
    });
}

void 
Maestro::MakePsiIrreg(RealVector& psi, const RealVector& grav_cell) 
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakePsiIrreg()", MakePsiIrreg);

    const int max_lev = max_radial_level+1;

    std::fill(psi.begin(), psi.end(), 0.0);

    const Real * AMREX_RESTRICT etarho_cc_p = etarho_cc.dataPtr();
    const Real * AMREX_RESTRICT grav_cell_p = grav_cell.dataPtr();
    Real * AMREX_RESTRICT psi_p = psi.dataPtr();

    const auto npts = base_cutoff_density_coord[0];
    AMREX_PARALLEL_FOR_1D(npts, r, {
        psi_p[max_lev*r] = etarho_cc_p[max_lev*r] * grav_cell_p[max_lev*r];
    });

    for (auto r = base_cutoff_density_coord[0]+1; r < nr_fine; ++r) {
        psi[max_lev*r] = psi[max_lev*(r-1)];
    }

    RestrictBase(psi, true);
    FillGhostBase(psi, true);
}
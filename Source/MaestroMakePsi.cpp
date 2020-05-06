#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void 
Maestro::MakePsiPlanar()
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakePsiPlanar()", MakePsiPlanar);

    const int max_lev = base_geom.max_radial_level + 1;

    psi.setVal(0.0);

    auto psi_arr = psi.array();
    auto etarho_cc_arr = etarho_cc.array();

    for (auto n = 0; n <= base_geom.finest_radial_level; ++n) {
        for (auto i = 1; i <= base_geom.numdisjointchunks(n); ++i){
            for (auto r = base_geom.r_start_coord(n,i); 
                 r<= base_geom.r_end_coord(n,i); ++r) {
                if (r < base_geom.base_cutoff_density_coord(n)) {
                    psi_arr(n,r) = etarho_cc_arr(n,r) * fabs(grav_const);
                }
            }
        }
    }
    RestrictBase(psi, true);
    FillGhostBase(psi, true);
}

void 
Maestro::MakePsiSphr(const BaseState<Real>& gamma1bar, 
                     const BaseState<Real>& p0_avg,
                     const RealVector& Sbar_in) 
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakePsiSphr()", MakePsiSphr);

    const int max_lev = base_geom.max_radial_level + 1;

    psi.setVal(0.0);

    Real dr0 = base_geom.dr(0);

    const auto& r_cc_loc = base_geom.r_cc_loc;
    const auto& r_edge_loc = base_geom.r_edge_loc;
    const Real * AMREX_RESTRICT w0_p = w0.dataPtr();
    const Real * AMREX_RESTRICT Sbar_p = Sbar_in.dataPtr();
    auto psi_arr = psi.array();
    const auto gamma1bar_arr = gamma1bar.const_array();
    const auto p0_avg_arr = p0_avg.const_array();

    const auto npts = base_geom.base_cutoff_density_coord(0);
    AMREX_PARALLEL_FOR_1D(npts, r, {
        Real div_w0_sph = 1.0 / (r_cc_loc(0,r)*r_cc_loc(0,r)) * 
            (r_edge_loc(0,r+1)*r_edge_loc(0,r+1) *
             w0_p[max_lev*(r+1)] - 
             r_edge_loc(0,r)*r_edge_loc(0,r) * 
             w0_p[max_lev*r]) / dr0;
        psi_arr(0,r) = gamma1bar_arr(0,r) * p0_avg_arr(0,r) * 
            (Sbar_p[max_lev*r] - div_w0_sph);
    });
    Gpu::synchronize();
}

void 
Maestro::MakePsiIrreg(const BaseState<Real>& grav_cell) 
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakePsiIrreg()", MakePsiIrreg);

    const int max_lev = base_geom.max_radial_level+1;

    psi.setVal(0.0);

    const auto etarho_cc_arr = etarho_cc.const_array();
    auto psi_arr = psi.array();
    const auto grav_cell_arr = grav_cell.const_array();

    const auto npts = base_geom.base_cutoff_density_coord(0);
    AMREX_PARALLEL_FOR_1D(npts, r, {
        psi_arr(0,r) = etarho_cc_arr(0,r) * grav_cell_arr(0,r);
    });
    Gpu::synchronize();

    for (auto r = base_geom.base_cutoff_density_coord(0)+1; r < base_geom.nr_fine; ++r) {
        psi_arr(0,r) = psi_arr(0,r-1);
    }

    RestrictBase(psi, true);
    FillGhostBase(psi, true);
}

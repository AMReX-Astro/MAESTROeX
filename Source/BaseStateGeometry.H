#ifndef BaseStateGeometry_H_
#define BaseStateGeometry_H_

#include <AMReX_AmrCore.H>
#include <AMReX_MultiFab.H>
#include <BaseState.H>

class BaseStateGeometry {
   public:
    BaseStateGeometry() noexcept {};

    void Init(const int max_radial_level_in, const int nr_fine_in,
              const amrex::Real dr_fine_in, const int nr_irreg_in,
              const amrex::Vector<amrex::Geometry>& geom, const int max_level,
              amrex::GpuArray<amrex::Real, 3>& center);

    void ComputeCutoffCoords(const BaseStateArray<amrex::Real>& rho0);

    void InitMultiLevel(const int finest_radial_level_in,
                        const BaseStateArray<int>& tag_array);

    int max_radial_level;
    int finest_radial_level;
    int nr_fine;
    amrex::Real dr_fine;
    int nr_irreg;

    BaseStateArray<int> nr;
    BaseStateArray<amrex::Real> dr;
    BaseStateArray<amrex::Real> r_cc_loc;
    BaseStateArray<amrex::Real> r_edge_loc;

    BaseStateArray<int> base_cutoff_density_coord;
    BaseStateArray<int> anelastic_cutoff_density_coord;
    BaseStateArray<int> burning_cutoff_density_lo_coord;
    BaseStateArray<int> burning_cutoff_density_hi_coord;

    ///these provide information about the multilevel base state configuration
    BaseStateArray<int> numdisjointchunks;
    BaseStateArray<int> r_start_coord;
    BaseStateArray<int> r_end_coord;

   private:
    // these contain the actual data and should not be
    // accessed directly other than by BaseStateGeometry class routines

    BaseState<int> nr_d;
    BaseState<amrex::Real> dr_d;
    BaseState<amrex::Real> r_cc_loc_d;
    BaseState<amrex::Real> r_edge_loc_d;

    BaseState<int> base_cutoff_density_coord_d;
    BaseState<int> anelastic_cutoff_density_coord_d;
    BaseState<int> burning_cutoff_density_lo_coord_d;
    BaseState<int> burning_cutoff_density_hi_coord_d;

    /// these provide information about the multilevel base state configuration
    BaseState<int> numdisjointchunks_d;
    BaseState<int> r_start_coord_d;
    BaseState<int> r_end_coord_d;
};

#endif

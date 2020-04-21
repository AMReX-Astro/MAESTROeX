#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void 
Maestro::AdvectBaseDens(RealVector& rho0_predicted_edge)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::AdvectBaseDens()", AdvectBaseDens); 

    if (spherical == 0) {
        AdvectBaseDensPlanar(rho0_predicted_edge);
        RestrictBase(rho0_new, true);
        FillGhostBase(rho0_new, true);
    } else {
        AdvectBaseDensSphr(rho0_predicted_edge);
    }
}

void 
Maestro::AdvectBaseDensPlanar(RealVector& rho0_predicted_edge)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::AdvectBaseDensPlanar()", AdvectBaseDensPlanar); 

    RealVector force_vec((max_radial_level+1)*nr_fine);

    // zero the new density so we don't leave a non-zero density in fine radial
    // regions that no longer have a corresponding full state
    std::fill(rho0_new.begin(), rho0_new.end(), 0.0);

    const int max_lev = max_radial_level+1;

    // Predict rho_0 to vertical edges

    for (int n = 0; n <= max_radial_level; ++n) {

        Real * AMREX_RESTRICT force = force_vec.dataPtr();
        Real * AMREX_RESTRICT rho0_old_p = rho0_old.dataPtr(); 
        Real * AMREX_RESTRICT w0_p = w0.dataPtr();  

        const Real dr_lev = dr.array()(n);

        for (int i = 1; i <= numdisjointchunks.array()(n); ++i) {

            int lo = r_start_coord.array()(n,i);
            int hi = r_end_coord.array()(n,i);

            AMREX_PARALLEL_FOR_1D(hi-lo+1, j, {
                int r = j + lo;
                int p = n+max_lev*r;

                force[p] = -rho0_old_p[p] * (w0_p[n+max_lev*(r+1)] - w0_p[p]) / dr_lev;
            });
        }
    }

    MakeEdgeState1d(rho0_old, rho0_predicted_edge, force_vec);

    for (int n = 0; n <= max_radial_level; ++n) {

        Real * AMREX_RESTRICT edge = rho0_predicted_edge.dataPtr();
        Real * AMREX_RESTRICT rho0_old_p = rho0_old.dataPtr(); 
        Real * AMREX_RESTRICT rho0_new_p = rho0_new.dataPtr(); 
        Real * AMREX_RESTRICT w0_p = w0.dataPtr();  

        const Real dr_lev = dr.array()(n);
        const Real dt_loc = dt;
        
        for (int i = 1; i <= numdisjointchunks.array()(n); ++i) {

            int lo = r_start_coord.array()(n,i);
            int hi = r_end_coord.array()(n,i);

            AMREX_PARALLEL_FOR_1D(hi-lo+1, j, {
                int r = j + lo;
                int p = n+max_lev*r;

                rho0_new_p[p] = rho0_old_p[p] 
                    - dt_loc/dr_lev * (edge[n+max_lev*(r+1)]*w0_p[n+max_lev*(r+1)] - edge[p]*w0_p[p]);
            });
        }
    }
}

void 
Maestro::AdvectBaseDensSphr(RealVector& rho0_predicted_edge)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::AdvectBaseDensSphr()", AdvectBaseDensSphr);

    const Real dr0 = dr.array()(0);
    const Real dtdr = dt / dr0;
    RealVector force_vec(nr_fine);
    const int max_lev = max_radial_level+1;

    // Predict rho_0 to vertical edges
    Real * AMREX_RESTRICT force = force_vec.dataPtr();
    Real * AMREX_RESTRICT rho0_old_p = rho0_old.dataPtr(); 
    Real * AMREX_RESTRICT rho0_new_p = rho0_new.dataPtr();
    Real * AMREX_RESTRICT w0_p = w0.dataPtr(); 
    const BaseStateArray<Real> r_cc_loc_p = r_cc_loc_b.array();
    const BaseStateArray<Real> r_edge_loc_p = r_edge_loc_b.array();
    Real * AMREX_RESTRICT edge = rho0_predicted_edge.dataPtr();

    AMREX_PARALLEL_FOR_1D(nr_fine, r, {
        int p = max_lev*r;
        force[p] = -rho0_old_p[p] * (w0_p[max_lev*(r+1)] - w0_p[p]) / dr0 - 
            rho0_old_p[p]*(w0_p[p] + w0_p[max_lev*(r+1)])/r_cc_loc_p(0,r);
    });

    MakeEdgeState1d(rho0_old, rho0_predicted_edge, force_vec);

    AMREX_PARALLEL_FOR_1D(nr_fine, r, {
        int p = max_lev*r;
        rho0_new_p[p] = rho0_old_p[p] - dtdr/(r_cc_loc_p(0,r)*r_cc_loc_p(0,r)) * 
            (r_edge_loc_p(0,r+1)*r_edge_loc_p(0,r+1) * edge[max_lev*(r+1)] * w0_p[max_lev*(r+1)] - 
            r_edge_loc_p(0,r)*r_edge_loc_p(0,r) * edge[p] * w0_p[p]);
    });
}

void 
Maestro::AdvectBaseEnthalpy(RealVector& rhoh0_predicted_edge)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::AdvectBaseEnthalpy()", AdvectBaseEnthalpy); 

    if (spherical == 0) {
        AdvectBaseEnthalpyPlanar(rhoh0_predicted_edge);
        RestrictBase(rhoh0_new, true);
        FillGhostBase(rhoh0_new, true);
    } else {
        AdvectBaseEnthalpySphr(rhoh0_predicted_edge);
    }
}

void 
Maestro::AdvectBaseEnthalpyPlanar(RealVector& rhoh0_predicted_edge)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::AdvectBaseEnthalpyPlanar()", AdvectBaseEnthalpyPlanar); 

    RealVector force_vec((max_radial_level+1)*nr_fine);

    // zero the new enthalpy so we don't leave a non-zero enthalpy in fine radial
    // regions that no longer have a corresponding full state
    std::fill(rhoh0_new.begin(), rhoh0_new.end(), 0);

    const int max_lev = max_radial_level+1;

    // Update (rho h)_0

    for (int n = 0; n <= max_radial_level; ++n) {

        Real * AMREX_RESTRICT force = force_vec.dataPtr();
        Real * AMREX_RESTRICT rhoh0_old_p = rhoh0_old.dataPtr(); 
        Real * AMREX_RESTRICT w0_p = w0.dataPtr();  
        Real * AMREX_RESTRICT psi_p = psi.dataPtr();  

        const Real dr_lev = dr.array()(n);

        for (int i = 1; i <= numdisjointchunks.array()(n); ++i) {

            int lo = r_start_coord.array()(n,i);
            int hi = r_end_coord.array()(n,i);

            // here we predict (rho h)_0 on the edges
            AMREX_PARALLEL_FOR_1D(hi-lo+1, j, {
                int r = j + lo;
                int p = n+max_lev*r;
                force[p] = -rhoh0_old_p[p] * (w0_p[n+max_lev*(r+1)] - w0_p[p]) / dr_lev 
                    + psi_p[p];
            });
        }
    }

    MakeEdgeState1d(rhoh0_old, rhoh0_predicted_edge, force_vec);

    for (int n = 0; n <= max_radial_level; ++n) {

        Real * AMREX_RESTRICT edge = rhoh0_predicted_edge.dataPtr();
        Real * AMREX_RESTRICT rhoh0_old_p = rhoh0_old.dataPtr(); 
        Real * AMREX_RESTRICT rhoh0_new_p = rhoh0_new.dataPtr(); 
        Real * AMREX_RESTRICT w0_p = w0.dataPtr();  
        Real * AMREX_RESTRICT psi_p = psi.dataPtr(); 

        const Real dr_lev = dr.array()(n);
        const Real dt_loc = dt;

        for (int i = 1; i <= numdisjointchunks.array()(n); ++i) {

            int lo = r_start_coord.array()(n,i);
            int hi = r_end_coord.array()(n,i);

            // update (rho h)_0
            AMREX_PARALLEL_FOR_1D(hi-lo+1, j, {
                int r = j + lo;
                int p = n+max_lev*r;
                rhoh0_new_p[p] = rhoh0_old_p[p] 
                    - dt_loc/dr_lev * (edge[n+max_lev*(r+1)]*w0_p[n+max_lev*(r+1)] - edge[p]*w0_p[p]) + dt_loc * psi_p[p];
            });
        }
    }
}

void 
Maestro::AdvectBaseEnthalpySphr(RealVector& rhoh0_predicted_edge)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::AdvectBaseEnthalpySphr()", AdvectBaseEnthalpySphr);

    const Real dr0 = dr.array()(0);
    const Real dtdr = dt / dr0;
    const Real dt_loc = dt;

    const int max_lev = max_radial_level + 1;
    RealVector force_vec(nr_fine);

    // predict (rho h)_0 on the edges
    Real * AMREX_RESTRICT force = force_vec.dataPtr();
    Real * AMREX_RESTRICT rhoh0_old_p = rhoh0_old.dataPtr(); 
    Real * AMREX_RESTRICT rhoh0_new_p = rhoh0_new.dataPtr();
    Real * AMREX_RESTRICT w0_p = w0.dataPtr(); 
    const BaseStateArray<Real> r_cc_loc_p = r_cc_loc_b.array();
    const BaseStateArray<Real> r_edge_loc_p = r_edge_loc_b.array();
    Real * AMREX_RESTRICT edge = rhoh0_predicted_edge.dataPtr();
    Real * AMREX_RESTRICT psi_p = psi.dataPtr(); 

    AMREX_PARALLEL_FOR_1D(nr_fine, r, {
        int p = max_lev*r;
        force[p] = -rhoh0_old_p[p] * (w0_p[max_lev*(r+1)] - w0_p[p]) / dr0 - 
            rhoh0_old_p[p]*(w0_p[p] + w0_p[max_lev*(r+1)])/r_cc_loc_p(0,r) + psi_p[p];
    });

    MakeEdgeState1d(rhoh0_old, rhoh0_predicted_edge, force_vec);

    AMREX_PARALLEL_FOR_1D(nr_fine, r, {
        int p = max_lev*r;

        rhoh0_new_p[p] = rhoh0_old_p[p] - dtdr/(r_cc_loc_p(0,r)*r_cc_loc_p(0,r)) * 
            (r_edge_loc_p(0,r+1)*r_edge_loc_p(0,r+1) * edge[max_lev*(r+1)] * w0_p[max_lev*(r+1)] - 
            r_edge_loc_p(0,r)*r_edge_loc_p(0,r) * edge[p] * w0_p[p]) + 
            dt_loc * psi_p[p];
    });
}


void 
Maestro::AdvectBaseDens(BaseState<Real>& rho0_predicted_edge)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::AdvectBaseDens()", AdvectBaseDens); 
    
    rho0_predicted_edge.setVal(0.0);

    if (spherical == 0) {
        AdvectBaseDensPlanar(rho0_predicted_edge);
        RestrictBase(rho0_new, true);
        FillGhostBase(rho0_new, true);
    } else {
        AdvectBaseDensSphr(rho0_predicted_edge);
    }
}

void 
Maestro::AdvectBaseDensPlanar(BaseState<Real>& rho0_predicted_edge)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::AdvectBaseDensPlanar()", AdvectBaseDensPlanar); 

    RealVector force_vec((max_radial_level+1)*nr_fine);

    // zero the new density so we don't leave a non-zero density in fine radial
    // regions that no longer have a corresponding full state
    std::fill(rho0_new.begin(), rho0_new.end(), 0.0);

    const int max_lev = max_radial_level+1;

    // Predict rho_0 to vertical edges

    for (int n = 0; n <= max_radial_level; ++n) {

        Real * AMREX_RESTRICT force = force_vec.dataPtr();
        Real * AMREX_RESTRICT rho0_old_p = rho0_old.dataPtr(); 
        Real * AMREX_RESTRICT w0_p = w0.dataPtr();  

        const Real dr_lev = dr.array()(n);

        for (int i = 1; i <= numdisjointchunks.array()(n); ++i) {

            int lo = r_start_coord.array()(n,i);
            int hi = r_end_coord.array()(n,i);

            AMREX_PARALLEL_FOR_1D(hi-lo+1, j, {
                int r = j + lo;
                int p = n+max_lev*r;

                force[p] = -rho0_old_p[p] * (w0_p[n+max_lev*(r+1)] - w0_p[p]) / dr_lev;
            });
        }
    }

    MakeEdgeState1d(rho0_old, rho0_predicted_edge, force_vec);

    const BaseStateArray<Real> rho0_predicted_edge_arr = rho0_predicted_edge.array();
    
    for (int n = 0; n <= max_radial_level; ++n) {

        Real * AMREX_RESTRICT rho0_old_p = rho0_old.dataPtr(); 
        Real * AMREX_RESTRICT rho0_new_p = rho0_new.dataPtr(); 
        Real * AMREX_RESTRICT w0_p = w0.dataPtr();  

        const Real dr_lev = dr.array()(n);
        const Real dt_loc = dt;
        
        for (int i = 1; i <= numdisjointchunks.array()(n); ++i) {

            int lo = r_start_coord.array()(n,i);
            int hi = r_end_coord.array()(n,i);

            AMREX_PARALLEL_FOR_1D(hi-lo+1, j, {
                int r = j + lo;
                int p = n+max_lev*r;

                rho0_new_p[p] = rho0_old_p[p] 
                    - dt_loc/dr_lev * (rho0_predicted_edge_arr(n,r+1)*w0_p[n+max_lev*(r+1)] - rho0_predicted_edge_arr(n,r)*w0_p[p]);
            });
        }
    }
}

void 
Maestro::AdvectBaseDensSphr(BaseState<Real>& rho0_predicted_edge)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::AdvectBaseDensSphr()", AdvectBaseDensSphr);

    const Real dr0 = dr.array()(0);
    const Real dtdr = dt / dr0;
    RealVector force_vec(nr_fine);
    const int max_lev = max_radial_level+1;

    // Predict rho_0 to vertical edges
    Real * AMREX_RESTRICT force = force_vec.dataPtr();
    Real * AMREX_RESTRICT rho0_old_p = rho0_old.dataPtr(); 
    Real * AMREX_RESTRICT rho0_new_p = rho0_new.dataPtr();
    Real * AMREX_RESTRICT w0_p = w0.dataPtr(); 
    const auto r_cc_loc_p = r_cc_loc_b.array();
    const auto r_edge_loc_p = r_edge_loc_b.array();

    AMREX_PARALLEL_FOR_1D(nr_fine, r, {
        int p = max_lev*r;
        force[p] = -rho0_old_p[p] * (w0_p[max_lev*(r+1)] - w0_p[p]) / dr0 - 
            rho0_old_p[p]*(w0_p[p] + w0_p[max_lev*(r+1)])/r_cc_loc_p(0,r);
    });

    MakeEdgeState1d(rho0_old, rho0_predicted_edge, force_vec);

    const BaseStateArray<Real> rho0_predicted_edge_arr = rho0_predicted_edge.array();
    
    AMREX_PARALLEL_FOR_1D(nr_fine, r, {
        int p = max_lev*r;
        rho0_new_p[p] = rho0_old_p[p] - dtdr/(r_cc_loc_p(0,r)*r_cc_loc_p(0,r)) * 
            (r_edge_loc_p(0,r+1)*r_edge_loc_p(0,r+1) * rho0_predicted_edge_arr(0,r+1) * w0_p[max_lev*(r+1)] - 
            r_edge_loc_p(0,r)*r_edge_loc_p(0,r) * rho0_predicted_edge_arr(0,r) * w0_p[p]);
    });
}

void 
Maestro::AdvectBaseEnthalpy(BaseState<Real>& rhoh0_predicted_edge)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::AdvectBaseEnthalpy()", AdvectBaseEnthalpy); 
    
    rhoh0_predicted_edge.setVal(0.0);

    if (spherical == 0) {
        AdvectBaseEnthalpyPlanar(rhoh0_predicted_edge);
        RestrictBase(rhoh0_new, true);
        FillGhostBase(rhoh0_new, true);
    } else {
        AdvectBaseEnthalpySphr(rhoh0_predicted_edge);
    }
}

void 
Maestro::AdvectBaseEnthalpyPlanar(BaseState<Real>& rhoh0_predicted_edge)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::AdvectBaseEnthalpyPlanar()", AdvectBaseEnthalpyPlanar); 

    RealVector force_vec((max_radial_level+1)*nr_fine);

    // zero the new enthalpy so we don't leave a non-zero enthalpy in fine radial
    // regions that no longer have a corresponding full state
    std::fill(rhoh0_new.begin(), rhoh0_new.end(), 0);

    const int max_lev = max_radial_level+1;

    // Update (rho h)_0

    for (int n = 0; n <= max_radial_level; ++n) {

        Real * AMREX_RESTRICT force = force_vec.dataPtr();
        Real * AMREX_RESTRICT rhoh0_old_p = rhoh0_old.dataPtr(); 
        Real * AMREX_RESTRICT w0_p = w0.dataPtr();  
        Real * AMREX_RESTRICT psi_p = psi.dataPtr();  

        const Real dr_lev = dr.array()(n);

        for (int i = 1; i <= numdisjointchunks.array()(n); ++i) {

            int lo = r_start_coord.array()(n,i);
            int hi = r_end_coord.array()(n,i);

            // here we predict (rho h)_0 on the edges
            AMREX_PARALLEL_FOR_1D(hi-lo+1, j, {
                int r = j + lo;
                int p = n+max_lev*r;
                force[p] = -rhoh0_old_p[p] * (w0_p[n+max_lev*(r+1)] - w0_p[p]) / dr_lev 
                    + psi_p[p];
            });
        }
    }

    MakeEdgeState1d(rhoh0_old, rhoh0_predicted_edge, force_vec);

    for (int n = 0; n <= max_radial_level; ++n) {

        Real * AMREX_RESTRICT rhoh0_old_p = rhoh0_old.dataPtr(); 
        Real * AMREX_RESTRICT rhoh0_new_p = rhoh0_new.dataPtr(); 
        Real * AMREX_RESTRICT w0_p = w0.dataPtr();  
        Real * AMREX_RESTRICT psi_p = psi.dataPtr();
        
        const BaseStateArray<Real> rhoh0_predicted_edge_arr = rhoh0_predicted_edge.array();

        const Real dr_lev = dr.array()(n);
        const Real dt_loc = dt;

        for (int i = 1; i <= numdisjointchunks.array()(n); ++i) {

            int lo = r_start_coord.array()(n,i);
            int hi = r_end_coord.array()(n,i);

            // update (rho h)_0
            AMREX_PARALLEL_FOR_1D(hi-lo+1, j, {
                int r = j + lo;
                int p = n+max_lev*r;
                rhoh0_new_p[p] = rhoh0_old_p[p] 
                    - dt_loc/dr_lev * (rhoh0_predicted_edge_arr(n,r+1)*w0_p[n+max_lev*(r+1)] - rhoh0_predicted_edge_arr(n,r)*w0_p[p]) + dt_loc * psi_p[p];
            });
        }
    }
}

void 
Maestro::AdvectBaseEnthalpySphr(BaseState<Real>& rhoh0_predicted_edge)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::AdvectBaseEnthalpySphr()", AdvectBaseEnthalpySphr);

    const Real dr0 = dr.array()(0);
    const Real dtdr = dt / dr0;
    const Real dt_loc = dt;

    const int max_lev = max_radial_level + 1;
    RealVector force_vec(nr_fine);

    // predict (rho h)_0 on the edges
    Real * AMREX_RESTRICT force = force_vec.dataPtr();
    Real * AMREX_RESTRICT rhoh0_old_p = rhoh0_old.dataPtr(); 
    Real * AMREX_RESTRICT rhoh0_new_p = rhoh0_new.dataPtr();
    Real * AMREX_RESTRICT w0_p = w0.dataPtr(); 
    const auto r_cc_loc_p = r_cc_loc_b.array();
    const auto r_edge_loc_p = r_edge_loc_b.array();
    Real * AMREX_RESTRICT psi_p = psi.dataPtr(); 

    AMREX_PARALLEL_FOR_1D(nr_fine, r, {
        int p = max_lev*r;
        force[p] = -rhoh0_old_p[p] * (w0_p[max_lev*(r+1)] - w0_p[p]) / dr0 - 
            rhoh0_old_p[p]*(w0_p[p] + w0_p[max_lev*(r+1)])/r_cc_loc_p(0,r) + psi_p[p];
    });

    MakeEdgeState1d(rhoh0_old, rhoh0_predicted_edge, force_vec);

    BaseStateArray<Real> rhoh0_predicted_edge_arr = rhoh0_predicted_edge.array();
    
    AMREX_PARALLEL_FOR_1D(nr_fine, r, {
        int p = max_lev*r;

        rhoh0_new_p[p] = rhoh0_old_p[p] - dtdr/(r_cc_loc_p(0,r)*r_cc_loc_p(0,r)) * 
            (r_edge_loc_p(0,r+1)*r_edge_loc_p(0,r+1) * rhoh0_predicted_edge_arr(0,r+1) * w0_p[max_lev*(r+1)] - 
            r_edge_loc_p(0,r)*r_edge_loc_p(0,r) * rhoh0_predicted_edge_arr(0,r) * w0_p[p]) + 
            dt_loc * psi_p[p];
    });
}

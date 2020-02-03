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
        // restrict base 
        // fill_ghost_base
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
    // RealVector edge_vec((max_radial_level+1)*(nr_fine+1));

    // zero the new density so we don't leave a non-zero density in fine radial
    // regions that no longer have a corresponding full state
    std::fill(rho0_new.begin(), rho0_new.end(), 0);

    const int max_lev = max_radial_level+1;

    // Predict rho_0 to vertical edges

    for (int n = 0; n <= max_radial_level; ++n) {

        Real * AMREX_RESTRICT force = force_vec.dataPtr();
        Real * AMREX_RESTRICT rho0_old_p = rho0_old.dataPtr(); 
        Real * AMREX_RESTRICT w0_p = w0.dataPtr();  

        const Real dr_lev = dr_fine * pow(2.0, n);

        for (int i = 0; i < numdisjointchunks[n]; ++i) {

            int lo = r_start_coord[n + max_lev*i];
            int hi = r_end_coord[n + max_lev*i];

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

        const Real dr_lev = dr_fine * pow(2.0, n);
        const Real dt_loc = dt;

        for (int i = 0; i < numdisjointchunks[n]; ++i) {

            int lo = r_start_coord[n + max_lev*i];
            int hi = r_end_coord[n + max_lev*i];

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

    const Real dr = dr_fine * pow(2.0, max_radial_level);
    const Real dtdr = dt / dr;
    RealVector force_vec(nr_fine);
    const int max_lev = max_radial_level+1;

    // Predict rho_0 to vertical edges
    Real * AMREX_RESTRICT force = force_vec.dataPtr();
    Real * AMREX_RESTRICT rho0_old_p = rho0_old.dataPtr(); 
    Real * AMREX_RESTRICT rho0_new_p = rho0_new.dataPtr();
    Real * AMREX_RESTRICT w0_p = w0.dataPtr(); 
    Real * AMREX_RESTRICT r_cc_loc_p = r_cc_loc.dataPtr();
    Real * AMREX_RESTRICT r_edge_loc_p = r_edge_loc.dataPtr();
    Real * AMREX_RESTRICT edge = rho0_predicted_edge.dataPtr();

    AMREX_PARALLEL_FOR_1D(nr_fine, r, {
        int p = max_lev*r;
        force[p] = -rho0_old_p[p] * (w0_p[max_lev*(r+1)] - w0_p[p]) / dr - 
            rho0_old_p[p]*(w0_p[p] + w0_p[max_lev*(r+1)])/r_cc_loc_p[p];
    });

    MakeEdgeState1d(rho0_old, rho0_predicted_edge, force_vec);

    AMREX_PARALLEL_FOR_1D(nr_fine, r, {
        int p = max_lev*r;
        rho0_new_p[p] = rho0_old_p[p] - dtdr/(r_cc_loc_p[p]*r_cc_loc_p[p]) * 
            (r_edge_loc_p[max_lev*(r+1)]*r_edge_loc_p[max_lev*(r+1)] * edge[max_lev*(r+1)] * w0_p[max_lev*(r+1)] - 
            r_edge_loc_p[p]*r_edge_loc_p[p] * edge[p] * w0_p[p]);
    });
}

void 
Maestro::AdvectBaseEnthalpy(RealVector& rhoh0_predicted_edge)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::AdvectBaseEnthalpy()", AdvectBaseEnthalpy); 

    if (spherical == 0) {
        AdvectBaseEnthalpyPlanar(rhoh0_predicted_edge);
        // restrict base 
        // fill_ghost_base
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
    // RealVector edge_vec((max_radial_level+1)*(nr_fine+1));

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

        const Real dr_lev = dr_fine * pow(2.0, n);

        for (int i = 0; i < numdisjointchunks[n]; ++i) {

            int lo = r_start_coord[n + max_lev*i];
            int hi = r_end_coord[n + max_lev*i];

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

        const Real dr_lev = dr_fine * pow(2.0, n);
        const Real dt_loc = dt;

        for (int i = 0; i < numdisjointchunks[n]; ++i) {

            int lo = r_start_coord[n + max_lev*i];
            int hi = r_end_coord[n + max_lev*i];

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

    const Real dr = dr_fine * pow(2.0, max_radial_level);
    const Real dtdr = dt / dr;
    const Real dt_loc = dt;

    const int max_lev = max_radial_level + 1;
    RealVector force_vec(nr_fine);

    // predict (rho h)_0 on the edges
    Real * AMREX_RESTRICT force = force_vec.dataPtr();
    Real * AMREX_RESTRICT rhoh0_old_p = rhoh0_old.dataPtr(); 
    Real * AMREX_RESTRICT rhoh0_new_p = rhoh0_new.dataPtr();
    Real * AMREX_RESTRICT w0_p = w0.dataPtr(); 
    Real * AMREX_RESTRICT r_cc_loc_p = r_cc_loc.dataPtr();
    Real * AMREX_RESTRICT r_edge_loc_p = r_edge_loc.dataPtr();
    Real * AMREX_RESTRICT edge = rhoh0_predicted_edge.dataPtr();
    Real * AMREX_RESTRICT psi_p = psi.dataPtr(); 

    AMREX_PARALLEL_FOR_1D(nr_fine, r, {
        int p = max_lev*r;

        force[p] = -rhoh0_old_p[p] * (w0_p[max_lev*(r+1)] - w0_p[p]) / dr - 
            rhoh0_old_p[p]*(w0_p[p] + w0_p[max_lev*(r+1)])/r_cc_loc_p[p] + psi_p[p];
    });

    MakeEdgeState1d(rhoh0_old, rhoh0_predicted_edge, force_vec);

    AMREX_PARALLEL_FOR_1D(nr_fine, r, {
        int p = max_lev*r;
        rhoh0_new_p[p] = rhoh0_old_p[p] - dtdr/(r_cc_loc_p[p]*r_cc_loc_p[p]) * 
            (r_edge_loc_p[max_lev*(r+1)]*r_edge_loc_p[max_lev*(r+1)] * edge[max_lev*(r+1)] * w0_p[max_lev*(r+1)] - 
            r_edge_loc_p[p]*r_edge_loc_p[p] * edge[p] * w0_p[p]) + 
            dt_loc * psi_p[p];
    });
}
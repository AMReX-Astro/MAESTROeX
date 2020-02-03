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
    rho0_new.setVal(0.);

    // Predict rho_0 to vertical edges

    for (int n = 0; n <= max_radial_level; ++n) {

        Real * AMREX_RESTRICT force = force_vec.dataPtr(n);
        Real * AMREX_RESTRICT rho0_old_p = rho0_old.dataPtr(n); 
        Real * AMREX_RESTRICT w0_p = w0.dataPtr(n);  

        const Real dr_lev = dr[n];

        for (int i = 0; i < numdisjointchunks[n]; ++i) {

            int lo = r_start_coord[n + i*(max_radial_level+1)];
            int hi = r_end_coord[n + i*(max_radial_level+1)];

            AMREX_PARALLEL_1D(hi-lo+1, j, {
                int r = j + lo;

                force[r] = -rho0_old_p[r] * (w0_p[r+1] - w0_p[r]) / dr_lev;
            });
        }
    }

    MakeEdgeState1d(rho0_old, rho0_predicted_edge, w0, force_vec);

    for (int n = 0; n <= max_radial_level; ++n) {

        Real * AMREX_RESTRICT edge = rho0_predicted_edge.dataPtr(n);
        Real * AMREX_RESTRICT rho0_old_p = rho0_old.dataPtr(n); 
        Real * AMREX_RESTRICT rho0_new_p = rho0_new.dataPtr(n); 
        Real * AMREX_RESTRICT w0_p = w0.dataPtr(n);  

        const Real dr_lev = dr[n];
        const Real dt_loc = dt;

        for (int i = 0; i < numdisjointchunks[n]; ++i) {

            int lo = r_start_coord[n + i*(max_radial_level+1)];
            int hi = r_end_coord[n + i*(max_radial_level+1)];

            AMREX_PARALLEL_1D(hi-lo+1, j, {
                int r = j + lo;

                rho0_new_p[r] = rho0_old_p[r] 
                    - dt_loc/dr_lev * (edge[r+1]*w0_p[r+1] - edge[r]*w0_p[r]);
            });
        }
    }
}

void 
Maestro::AdvectBaseDensSphr(RealVector& rho0_predicted_edge)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::AdvectBaseDensSphr()", AdvectBaseDensSphr); 

    Real dtdr = dt / dr[0];
    Real dr_loc = dr[0];

    // Predict rho_0 to vertical edges
    Real * AMREX_RESTRICT force = force_vec.dataPtr();
    Real * AMREX_RESTRICT rho0_old_p = rho0_old.dataPtr(); 
    Real * AMREX_RESTRICT rho0_new_p = rho0_new.dataPtr();
    Real * AMREX_RESTRICT w0_p = w0.dataPtr(); 
    Real * AMREX_RESTRICT r_cc_loc_p = r_cc_loc.dataPtr();
    Real * AMREX_RESTRICT edge = rho0_predicted_edge.dataPtr();

    AMREX_PARALLEL_FOR_1D(nr_fine, r, {
        force[r] = -rho0_old_p[r] * (w0_p[r+1] - w0_p[r]) / dr_loc - 
            rho0_old_p[r]*(w0_p[r] + w0_p[r+1])/r_cc_loc_p[r];
    });

    MakeEdgeState1d(rho0_old, rho0_predicted_edge, w0, force_vec);

    AMREX_PARALLEL_FOR_1D(nr_fine, r, {
        rho0_new_p[r] = rho0_old_p[r] - dtdr/(r_cc_loc_p[r]*r_cc_loc_p[r]) * 
            (r_edge_loc_p[r+1]*r_edge_loc_p[r+1] * edge[r+1] * w0_p[r+1] - 
            r_edge_loc_p[r]*r_edge_loc_p[r] * edge[r] * w0_p[r]);
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
    rhoh0_new.setVal(0.);

    // Update (rho h)_0

    for (int n = 0; n <= max_radial_level; ++n) {

        Real * AMREX_RESTRICT force = force_vec.dataPtr(n);
        Real * AMREX_RESTRICT rhoh0_old_p = rhoh0_old.dataPtr(n); 
        Real * AMREX_RESTRICT w0_p = w0.dataPtr(n);  
        Real * AMREX_RESTRICT psi_p = psi.dataPtr(n);  

        const Real dr_lev = dr[n];

        for (int i = 0; i < numdisjointchunks[n]; ++i) {

            int lo = r_start_coord[n + i*(max_radial_level+1)];
            int hi = r_end_coord[n + i*(max_radial_level+1)];

            // here we predict (rho h)_0 on the edges
            AMREX_PARALLEL_1D(hi-lo+1, j, {
                int r = j + lo;
                force[r] = -rhoh0_old_p[r] * (w0_p[r+1] - w0_p[r]) / dr_lev 
                    + psi_p[r];
            });
        }
    }

    MakeEdgeState1d(rhoh0_old, rhoh0_predicted_edge, w0, force_vec);

    for (int n = 0; n <= max_radial_level; ++n) {

        Real * AMREX_RESTRICT edge = rhoh0_predicted_edge.dataPtr(n);
        Real * AMREX_RESTRICT rhoh0_old_p = rhoh0_old.dataPtr(n); 
        Real * AMREX_RESTRICT rhoh0_new_p = rhoh0_new.dataPtr(n); 
        Real * AMREX_RESTRICT w0_p = w0.dataPtr(n);  
        Real * AMREX_RESTRICT psi_p = psi.dataPtr(n); 

        const Real dr_lev = dr[n];
        const Real dt_loc = dt;

        for (int i = 0; i < numdisjointchunks[n]; ++i) {

            int lo = r_start_coord[n + i*(max_radial_level+1)];
            int hi = r_end_coord[n + i*(max_radial_level+1)];

            // update (rho h)_0
            AMREX_PARALLEL_1D(hi-lo+1, j, {
                int r = j + lo;
                rhoh0_new_p[r] = rhoh0_old_p[r] 
                    - dt_loc/dr_lev * (edge[r+1]*w0_p[r+1] - edge[r]*w0_p[r]) + dt_loc * psi_p[r];
            });
        }
    }
}

void 
Maestro::AdvectBaseEnthalpySphr(RealVector& rhoh0_predicted_edge)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::AdvectBaseEnthalpySphr()", AdvectBaseEnthalpySphr);

    Real dtdr = dt / dr[0];
    Real dr_loc = dr[0];
    Real dt_loc = dt;

    // predict (rho h)_0 on the edges
    Real * AMREX_RESTRICT force = force_vec.dataPtr();
    Real * AMREX_RESTRICT rhho0_old_p = rhoh0_old.dataPtr(); 
    Real * AMREX_RESTRICT rhoh0_new_p = rhoh0_new.dataPtr();
    Real * AMREX_RESTRICT w0_p = w0.dataPtr(); 
    Real * AMREX_RESTRICT r_cc_loc_p = r_cc_loc.dataPtr();
    Real * AMREX_RESTRICT edge = rho0_predicted_edge.dataPtr();
    Real * AMREX_RESTRICT psi_p = psi.dataPtr(); 

    AMREX_PARALLEL_FOR_1D(nr_fine, r, {
        force[r] = -rhoh0_old_p[r] * (w0_p[r+1] - w0_p[r]) / dr_loc - 
            rhoh0_old_p[r]*(w0_p[r] + w0_p[r+1])/r_cc_loc_p[r] + psi_p[r];
    });

    MakeEdgeState1d(rhoh0_old, rhoh0_predicted_edge, w0, force_vec);

    AMREX_PARALLEL_FOR_1D(nr_fine, r, {
        rhoh0_new_p[r] = rhoh0_old_p[r] - dtdr/(r_cc_loc_p[r]*r_cc_loc_p[r]) * 
            (r_edge_loc_p[r+1]*r_edge_loc_p[r+1] * edge[r+1] * w0_p[r+1] - 
            r_edge_loc_p[r]*r_edge_loc_p[r] * edge[r] * w0_p[r]) + 
            dt_loc * psi_p[r];
    });
}
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void 
Maestro::RestrictBase(RealVector& s0_vec, bool is_cell_centered)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::RestrictBase()", RestrictBase); 

    get_numdisjointchunks(numdisjointchunks.dataPtr());
    get_r_start_coord(r_start_coord.dataPtr());
    get_r_end_coord(r_end_coord.dataPtr());

    const int max_lev = max_radial_level + 1;

    for (int n = max_radial_level; n >= 1; --n) {
        Real * AMREX_RESTRICT s0 = s0_vec.dataPtr();
        // Real * AMREX_RESTRICT s0_m = s0_vec[n-1].dataPtr();
        for (int i = 0; i < numdisjointchunks[n]; ++i) {

            if (is_cell_centered) {
                // for level n, make the coarser cells underneath simply the average of the fine
                const int lo = r_start_coord[n+i*max_lev];
                const int hi = r_end_coord[n+i*max_lev]-1;

                AMREX_PARALLEL_FOR_1D((hi-lo+1)/2+1, j, {
                    // int r = j*2 + lo;
                    // s0_m[r/2] = 0.5 * (s0[r] + s0[r+1]);
                    int r = n + max_lev*(j*2 + lo);
                    int p = n + max_lev*(j*2 + lo + 1);
                    int m = n-1 + max_lev*(j*2 + lo)/2;

                    s0[m] = 0.5 * (s0[r] + s0[p]);
                });
            } else {
                // for level n, make the coarse edge underneath equal to the fine edge value
                const int lo = r_start_coord[n+i*max_lev];
                const int hi = r_end_coord[n+i*max_lev]+1;

                AMREX_PARALLEL_FOR_1D((hi-lo+1)/2+1, j, {
                    // int r = j*2 + lo;
                    // s0_m[r/2] = s0[r];
                    int r = n + max_lev*(j*2 + lo);
                    int m = n-1 + max_lev*(j*2 + lo)/2;

                    s0[m] = s0[r];
                });
            }
        }
    }
}

void 
Maestro::FillGhostBase(RealVector& s0_vec, bool is_cell_centered)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::FillGhostBase()", FillGhostBase); 

    get_numdisjointchunks(numdisjointchunks.dataPtr());
    get_r_start_coord(r_start_coord.dataPtr());
    get_r_end_coord(r_end_coord.dataPtr());

    const int max_lev = max_radial_level + 1;

    for (int n = max_radial_level; n >= 1; --n) {
        Real * AMREX_RESTRICT s0 = s0_vec.dataPtr();
        const int nr = nr_fine / pow(2,max_radial_level-n);

        // Real * AMREX_RESTRICT s0_m = s0_vec[n-1].dataPtr();
        for (int i = 0; i < numdisjointchunks[n]; ++i) {

            const int lo = r_start_coord[n+i*max_lev];
            const int hi = r_end_coord[n+i*max_lev];

            if (is_cell_centered) {

                if (lo != 0) {
                    int r_crse = lo/2-1;
                    Real del = 0.5*(s0[n-1+max_lev*(r_crse+1)]-s0[n-1+max_lev*(r_crse-1)]);
                    Real dpls = 2.0*(s0[n-1+max_lev*(r_crse+1)]-s0[n-1+max_lev*r_crse]);
                    Real dmin = 2.0*(s0[n-1+max_lev*r_crse]-s0[n-1+max_lev*(r_crse-1)]);
                    Real slim = min(fabs(dpls),fabs(dmin));
                    slim = dpls*dmin > 0.0 ? slim : 0.0;
                    Real slope = copysign(1.0,del)*min(slim,fabs(del));
                    s0[n+max_lev*(lo-1)] = s0[n-1+max_lev*r_crse] + 0.25*slope;
                    s0[n+max_lev*(lo-2)] = s0[n-1+max_lev*r_crse] - 0.25*slope;

                    r_crse = lo/2-2;
                    del = 0.5*(s0[n-1+max_lev*(r_crse+1)]-s0[n-1+max_lev*(r_crse-1)]);
                    dpls = 2.0*(s0[n-1+max_lev*(r_crse+1)]-s0[n-1+max_lev*r_crse]);
                    dmin = 2.0*(s0[n-1+max_lev*r_crse]-s0[n-1+max_lev*(r_crse-1)]);
                    slim = min(fabs(dpls),fabs(dmin));
                    slim = dpls*dmin > 0.0 ? slim : 0.0;
                    slope = copysign(1.0,del)*min(slim,fabs(del));
                    s0[n+max_lev*(lo-3)] = s0[n-1+max_lev*r_crse] + 0.25*slope;
                    s0[n+max_lev*(lo-4)] = s0[n-1+max_lev*r_crse] - 0.25*slope;
                }

                if (hi != nr-1) {
                    int r_crse = (hi+1)/2;
                    Real del = 0.5*(s0[n-1+max_lev*(r_crse+1)]-s0[n-1+max_lev*(r_crse-1)]);
                    Real dpls = 2.0*(s0[n-1+max_lev*(r_crse+1)]-s0[n-1+max_lev*r_crse]);
                    Real dmin = 2.0*(s0[n-1+max_lev*r_crse]-s0[n-1+max_lev*(r_crse-1)]);
                    Real slim = min(fabs(dpls),fabs(dmin));
                    slim = dpls*dmin > 0.0 ? slim : 0.0;
                    Real slope = copysign(1.0,del)*min(slim,fabs(del));
                    s0[n+max_lev*(hi+1)] = s0[n-1+max_lev*r_crse] - 0.25*slope;
                    s0[n+max_lev*(hi+2)] = s0[n-1+max_lev*r_crse] + 0.25*slope;

                    r_crse = (hi+1)/2+1;
                    del = 0.5*(s0[n-1+max_lev*(r_crse+1)]-s0[n-1+max_lev*(r_crse-1)]);
                    dpls = 2.0*(s0[n-1+max_lev*(r_crse+1)]-s0[n-1+max_lev*r_crse]);
                    dmin = 2.0*(s0[n-1+max_lev*r_crse]-s0[n-1+max_lev*(r_crse-1)]);
                    slim = min(fabs(dpls),fabs(dmin));
                    slim = dpls*dmin > 0.0 ? slim : 0.0;
                    slope = copysign(1.0,del)*min(slim,fabs(del));
                    s0[n+max_lev*(hi+3)] = s0[n-1+max_lev*r_crse] - 0.25*slope;
                    s0[n+max_lev*(hi+4)] = s0[n-1+max_lev*r_crse] + 0.25*slope;
                }
                    
            } else {

                if (lo != 0) {
                    // quadratic interpolation from the three closest points
                    s0[n+max_lev*(lo-1)] = -s0[n+max_lev*(lo+1)]/3.0
                        + s0[n+max_lev*lo] + s0[n-1+max_lev*(lo/2-1)]/3.0;
                    // copy the next ghost cell value directly in from the coarser level
                    s0[n+max_lev*(lo-2)] = s0[n-1+max_lev*((lo-2)/2)];
                }

                if (hi+1 != nr) {
                    // quadratic interpolation from the three closest points
                    s0[n+max_lev*(hi+2)] = -s0[n+max_lev*hi]/3.0
                        + s0[n+max_lev*(hi+1)] + s0[n-1+max_lev*((hi+3)/2)]/3.0;
                    // copy the next ghost cell value directly in from the coarser level
                    s0[n+max_lev*(hi+3)] = s0[n-1+max_lev*((hi+3)/2)];
                }
            }
        }
    }
}
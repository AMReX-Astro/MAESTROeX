#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void Maestro::MakeEdgeState1d(RealVector& s, RealVector& sedge, 
                              RealVector& w0, RealVector& force_1d)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeEdgeState1d()",MakeEdgeState1d);

    get_numdisjointchunks(numdisjointchunks.dataPtr());

    

    if (spherical) {
        MakeEdgeState1dSphr(s, sedge, w0, force_1d);
    } else {
        // MakeEdgeState1dPlanar(s_p, sedge_p, w0_p, force_p);
    }
}

void Maestro::MakeEdgeState1dSphr(RealVector& s_vec, RealVector& sedge_vec, 
                              RealVector& w0, RealVector& force_1d)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeEdgeState1dSphr()",MakeEdgeState1dSphr);

    Real rel_eps = 0.0;
    get_rel_eps(&rel_eps);

    const Real dth = 0.5 * dt;
    const Real C = 1.25;
    const int cen = 0;
    const int lim = 1;
    const int flag = 2;
    const int fromm = 3;
    const Real dr = dr_fine;

    const int nr = nr_fine;

    RealVector sedgel_vec(nr_fine+1);
    RealVector sedger_vec(nr_fine+1);

    // copy valid data into array with ghost cells
    RealVector s_ghost_vec(nr_fine+2*3);
    for (int i = 0; i < s_vec.size(); ++i) {
        s_ghost_vec[i+3] = s_vec[i];
    }

    for (int i = 0; i < 3; i++) {
        // symmetry boundary condition at center 
        s_ghost_vec[2-i] = s_vec[i];
        // first-order extrapolation at top of star
        s_ghost_vec[3+nr_fine+i] = s_vec[nr_fine-1];
    }

    Real * AMREX_RESTRICT sedgel = sedgel_vec.dataPtr();
    Real * AMREX_RESTRICT sedger = sedger_vec.dataPtr();
    Real * AMREX_RESTRICT s_ghost = s_ghost_vec.dataPtr();

    Real * AMREX_RESTRICT s = s_vec.dataPtr();
    Real * AMREX_RESTRICT sedge = sedge_vec.dataPtr();
    Real * AMREX_RESTRICT w0_p = w0.dataPtr();
    Real * AMREX_RESTRICT force_p = force_1d.dataPtr();

    if (ppm_type == 0) {

        // this will hold values at r-1, r and r+1
        GpuArray<GpuArray<Real,4>,3> dsscr;

        AMREX_PARALLEL_FOR_1D(nr_fine, r, {

            Real slope = 0.0;

            if (slope_order == 0) {
                slope = 0.0;
            } else if (slope_order == 2) {
                // index of ghost array offset by 3
                int g = r + 3;

                Real del = 0.5 * (s_ghost[g+1] - s_ghost[g-1]);
                Real dpls = 2.0*(s_ghost[g+1] - s_ghost[g]);
                Real dmin = 2.0*(s_ghost[g] - s_ghost[g-1]);
                Real slim = min(fabs(dpls), fabs(dmin));
                slim = dpls*dmin > 0.0 ? slim : 0.0;
                Real sflag = copysign(1.0,del);
                slope = sflag*min(slim,fabs(del));
            } else if (slope_order == 4) {
                for (int i = 0; i < 3; ++i) {
                    int g = r + 3 - 1 + i;
                    // do standard limiting to compute temporary slopes
                    dsscr[i][cen] = 0.5*(s_ghost[g+1]-s_ghost[g-1]);
                    Real dpls = 2.0*(s_ghost[g+1]-s_ghost[g]);
                    Real dmin = 2.0*(s_ghost[g]-s_ghost[g-1]);
                    dsscr[i][lim]= min(fabs(dmin),fabs(dpls));
                    dsscr[i][lim] = dpls*dmin > 0.0 ? dsscr[i][lim] : 0.0;
                    dsscr[i][flag] = copysign(1.0,dsscr[i][cen]);
                    dsscr[i][fromm]= dsscr[i][flag]*min(dsscr[i][lim],fabs(dsscr[i][cen]));
                }

                // fourth-order limited slopes
                Real ds = 4.0/3.0*dsscr[1][cen] - (dsscr[2][fromm]+dsscr[0][fromm]) / 6.0;
                slope = dsscr[1][flag]*min(fabs(ds),dsscr[1][lim]);
            }

            // compute sedgel and sedger
            Real u = 0.5*(w0_p[r]+w0_p[r+1]);
            Real ubardth = dth*u/dr;  // NOTE: ubardth=0 for use_exact_base_state case
            sedgel[r+1 ]= s[r] + (0.5-ubardth)*slope + dth*force_p[r];
            sedger[r  ] = s[r] - (0.5+ubardth)*slope + dth*force_p[r];
        });

    } else if (ppm_type == 1) {

        AMREX_PARALLEL_FOR_1D(nr_fine, r, {

            // interpolate s to radial edges

            // sm

            // right side
            int g = r + 3; 
            // compute van Leer slopes
            Real del  = 0.5 * (s_ghost[g+1] - s_ghost[g-1]);
            Real dmin = 2.0  * (s_ghost[g] - s_ghost[g-1]);
            Real dpls = 2.0  * (s_ghost[g+1] - s_ghost[g]);
            Real dsscrr = dmin*dpls  >  0.0 ? copysign(1.0,del)*min(fabs(del),fabs(dmin),fabs(dpls)) : 0.0;

            // left side 
            g = r + 3 - 1;
            // compute van Leer slopes
            del  = 0.5 * (s_ghost[g+1] - s_ghost[g-1]);
            dmin = 2.0  * (s_ghost[g] - s_ghost[g-1]);
            dpls = 2.0  * (s_ghost[g+1] - s_ghost[g]);
            Real dsscrl = dmin*dpls  >  0.0 ? copysign(1.0,del)*min(fabs(del),fabs(dmin),fabs(dpls)) : 0.0;

            // sm
            g = r + 3;
            // 4th order interpolation of s to radial faces
            Real sm = 0.5*(s_ghost[g]+s_ghost[g-1])-(dsscrr-dsscrl) / 6.0;
            // make sure sedgel lies in between adjacent cell-centered values
            sm = max(sm,min(s_ghost[g],s_ghost[g-1]));
            sm = min(sm,max(s_ghost[g],s_ghost[g-1]));

            // sp

            // left side
            dsscrl = dsscrr;

            // right side
            g = r + 3 + 1; 
            // compute van Leer slopes
            del  = 0.5 * (s_ghost[g+1] - s_ghost[g-1]);
            dmin = 2.0  * (s_ghost[g] - s_ghost[g-1]);
            dpls = 2.0  * (s_ghost[g+1] - s_ghost[g]);
            dsscrr = dmin*dpls  >  0.0 ? copysign(1.0,del)*min(fabs(del),fabs(dmin),fabs(dpls)) : 0.0;

            // sp
            g = r + 3 + 1;
            // 4th order interpolation of s to radial faces
            Real sp = 0.5*(s_ghost[g]+s_ghost[g-1])-(dsscrr-dsscrl) / 6.0;
            // make sure sedgel lies in between adjacent cell-centered values
            sp = max(sm,min(s_ghost[g],s_ghost[g-1]));
            sp = min(sm,max(s_ghost[g],s_ghost[g-1]));

            // modify using quadratic limiters
            if ((sp-s[r])*(s[r]-sm) <= 0.0) {
                sp = s[r];
                sm = s[r];
            } else if (fabs(sp-s[r]) >= 2.0*fabs(sm-s[r])) {
                sp = 3.0*s[r] - 2.0*sm;
            } else if (fabs(sm-s[r]) >= 2.0*fabs(sp-s[r])) {
                sm = 3.0*s[r] - 2.0*sp;
            }

            // compute Ip and Im
            Real sigmap = fabs(w0_p[r+1])*dt/dr;  // NOTE: sigmap=0 for use_exact_base_state case
            Real sigmam = fabs(w0_p[r])*dt/dr;  // NOTE: sigmam=0 for use_exact_base_state case
            Real s6 = 6.0*s[r] - 3.0*(sm+sp);
            
            Real Ip = w0_p[r+1] > rel_eps ? sp - 0.5*sigmap*(sp-sm-(1.0-2.0/3.0*sigmap)*s6) : s[r];

            Real Im = w0_p[r] < -rel_eps ? sm + 0.5*sigmam*(sp-sm+(1.0-2.0/3.0*sigmam)*s6) : s[r];

            // ! compute sedgel and sedger
            sedgel[r+1] = Ip + dth*force_p[r];
            sedger[r] = Im + dth*force_p[r];
        });
    } else if (ppm_type == 2) {

        // this will hold values at r-1, r, r+1 and r+2
        GpuArray<Real,4> sedget;

        AMREX_PARALLEL_FOR_1D(nr_fine, r, {
            // interpolate s to radial edges, store these temporary values into sedgel

            for (int i = 0; i < 4; i++) {
                int g = r + 3 - 1 + i;
                int rr = r - 1 + i;

                sedget[r] = (7.0/12.0)*(s_ghost[g-1]+s_ghost[g]) 
                        - (1.0/12.0)*(s_ghost[g-2]+s_ghost[g+1]);

                // limit sedge
                if ((sedget[rr]-s_ghost[g-1])*(s_ghost[g]-sedget[rr]) < 0.0) {
                    Real D2  = 3.0*(s_ghost[g-1]-2.0*sedgel[r]+s_ghost[g]);
                    Real D2L = s_ghost[g-2]-2.0*s_ghost[g-1]+s_ghost[g];
                    Real D2R = s_ghost[g-1]-2.0*s_ghost[g]+s_ghost[g+1];
                    Real sgn = copysign(1.0,D2);
                    Real D2LIM = sgn*max(min(C*sgn*D2L,min(C*sgn*D2R,sgn*D2)),0.0);
                    sedget[rr] = 0.5*(s_ghost[g-1]+s_ghost[g]) - D2LIM/6.0;
                }
            }

            int g = r + 3;

            // use Colella 2008 limiters
            // This is a new version of the algorithm
            // to eliminate sensitivity to roundoff.
            Real alphap = sedget[2]-s_ghost[g];
            Real alpham = sedget[1]-s_ghost[g];
            bool bigp = fabs(alphap) > 2.0*fabs(alpham);
            bool bigm = fabs(alpham) > 2.0*fabs(alphap);
            bool extremum = false;

            if (alpham*alphap >= 0.0) {
                extremum = true;
            } else if (bigp || bigm) {
                // Possible extremum. We look at cell centered values and face
                // centered values for a change in copysign in the differences adjacent to
                // the cell. We use the pair of differences whose minimum magnitude is
                // the largest, and thus least susceptible to sensitivity to roundoff.
                Real dafacem = sedget[1] - sedget[0];
                Real dafacep = sedget[3] - sedget[2];
                Real dabarm = s_ghost[g] - s_ghost[g-1];
                Real dabarp = s_ghost[g+1] - s_ghost[g];
                Real dafacemin = min(fabs(dafacem),fabs(dafacep));
                Real dabarmin = min(fabs(dabarm),fabs(dabarp));
                Real dachkm = 0.0;
                Real dachkp = 0.0;
                if (dafacemin >= dabarmin) {
                    dachkm = dafacem;
                    dachkp = dafacep;
                } else {
                    dachkm = dabarm;
                    dachkp = dabarp;
                }
                extremum = (dachkm*dachkp <= 0.0);
            }

            if (extremum) {
                Real D2  = 6.0*(alpham + alphap);
                Real D2L = s_ghost[g-2]-2.0*s_ghost[g-1]+s_ghost[g];
                Real D2R = s_ghost[g]-2.0*s_ghost[g+1]+s_ghost[g+2];
                Real D2C = s_ghost[g-1]-2.0*s_ghost[g]+s_ghost[g+1];
                Real sgn = copysign(1.0,D2);
                Real D2LIM = max(min(sgn*D2,min(C*sgn*D2L,min(C*sgn*D2R,C*sgn*D2C))),0.0);
                Real D2ABS = max(fabs(D2),1.e-10);
                alpham = alpham*D2LIM/D2ABS;
                alphap = alphap*D2LIM/D2ABS;
            } else {
                if (bigp) {
                    Real sgn = copysign(1.0,alpham);
                    Real amax = -alphap*alphap / (4.0*(alpham + alphap));
                    Real delam = s_ghost[g-1] - s_ghost[g];
                    if (sgn*amax >= sgn*delam) {
                        if (sgn*(delam - alpham) >= 1.e-10) {
                            alphap = (-2.0*delam - 2.0*sgn*sqrt(delam*delam - delam*alpham));
                        } else {
                            alphap = -2.0*alpham;
                        }
                    }
                }
                if (bigm) {
                    Real sgn = copysign(1.0,alphap);
                    Real amax = -alpham*alpham / (4.0*(alpham + alphap));
                    Real delap = s_ghost[g+1] - s_ghost[g];
                    if (sgn*amax >= sgn*delap) {
                        if (sgn*(delap - alphap)>=1.e-10) {
                            alpham = (-2.0*delap - 2.0*sgn*sqrt(delap*delap - delap*alphap));
                        } else {
                            alpham = -2.0*alphap;
                        }
                    }
                }
            }

            Real sm = s_ghost[g] + alpham;
            Real sp = s_ghost[g] + alphap;

            // compute Ip and Im
            Real sigmap = fabs(w0_p[r+1])*dt/dr;  // NOTE: sigmap=0 for use_exact_base_state case
            Real sigmam = fabs(w0_p[r])*dt/dr;  // NOTE: sigmam=0 for use_exact_base_state case
            Real s6 = 6.0*s[r] - 3.0*(sm+sp);
            
            Real Ip = w0_p[r+1] > rel_eps ? sp - 0.5*sigmap*(sp-sm-(1.0-2.0/3.0*sigmap)*s6) : s[r];

            Real Im = w0_p[r] < -rel_eps ? sm + 0.5*sigmam*(sp-sm+(1.0-2.0/3.0*sigmam)*s6) : s[r];

            // ! compute sedgel and sedger
            sedgel[r+1] = Ip + dth*force_p[r];
            sedger[r] = Im + dth*force_p[r];
        });
    }

    AMREX_PARALLEL_FOR_1D(nr_fine+1, r, {
        // Fix center and edge of star by reflecting the extrapolated state.
        // An alternate way would be to compute these values using the entire algorithm,
        // but that would require more ghost cells at several stages.
        // By symmetry arguments, this would make no difference at the center of the star
        // and the accuracy at the edge of the star is not important here
        if (r == 0) {
            sedgel[r] = sedge[r];
        } else if (r == nr_fine) {
            sedger[r] = sedgel[r];
        }

        // solve Riemann problem to get final edge state
        sedge[r]= w0_p[r] > 0.0 ? sedgel[r] : sedger[r];
        sedge[r] = fabs(w0_p[r])<rel_eps ? 0.5*(sedger[r]+sedgel[r]) : sedge[r];
    });

}

void Maestro::MakeEdgeState1dPlanar(Real* s, Real* sedge, 
                                    Real* w0, Real* force)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeEdgeState1dPlanar()",MakeEdgeState1dPlanar);
    
}
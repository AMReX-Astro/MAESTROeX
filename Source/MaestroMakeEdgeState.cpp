#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void Maestro::MakeEdgeState1d(RealVector& s, RealVector& sedge, 
                              RealVector& w0, RealVector& force_1d)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeEdgeState1d()",MakeEdgeState1d);   

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

            // // compute sedgel and sedger
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

            // // compute sedgel and sedger
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

void Maestro::MakeEdgeState1dPlanar(RealVector& s_vec, RealVector& sedge_vec, 
                              RealVector& w0, RealVector& force_1d)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeEdgeState1dPlanar()",MakeEdgeState1dPlanar);

    get_numdisjointchunks(numdisjointchunks.dataPtr());
    get_r_start_coord(r_start_coord.dataPtr());
    get_r_end_coord(r_end_coord.dataPtr());

    Real rel_eps = 0.0;
    get_rel_eps(&rel_eps);

    const Real dth = 0.5 * dt;
    const Real C = 1.25;
    const int cen = 0;
    const int lim = 1;
    const int flag = 2;
    const int fromm = 3;
    

    Vector<RealVector> sedgel_vec(max_radial_level+1);
    Vector<RealVector> sedger_vec(max_radial_level+1);

    Real * AMREX_RESTRICT s = s_vec.dataPtr();
    Real * AMREX_RESTRICT sedge = sedge_vec.dataPtr();
    Real * AMREX_RESTRICT w0_p = w0.dataPtr();
    Real * AMREX_RESTRICT force_p = force_1d.dataPtr();


    for (int n = 0; n <= max_radial_level; ++n) {

        sedgel_vec[n].resize(nr_fine+1);
        sedger_vec[n].resize(nr_fine+1);

        Real * AMREX_RESTRICT sedgel = sedgel_vec[n].dataPtr();
        Real * AMREX_RESTRICT sedger = sedger_vec[n].dataPtr();   

        const int nr_lev = nr[n];
        const Real dr_lev = dr[n];
        const Real dtdr = dt / dr[n];

        for (int i = 0; i < numdisjointchunks[n]; ++i) {

            const int lo = r_start_coord[n+i*max_radial_level];
            const int hi = r_end_coord[n+i*max_radial_level]; 

            // error checking to make sure that there is a 2 cell buffer at the top and bottom
            // of the domain for finer levels in planar geometry.  This can be removed if
            // blocking_factor is implemented at set > 1.
            if (ppm_type == 1 || ppm_type == 2) {
                
                if (r_start_coord[n+i*max_radial_level] == 2) {
                    Abort("make_edge_state assumes blocking_factor > 1 at lo boundary");
                } else if (r_end_coord[n+i*max_radial_level] == nr_lev-3) {
                    Abort("make_edge_state assumes blocking_factor > 1 at hi boundary");
                }
            }

            if (ppm_type == 0) {

                GpuArray<Real, 3> dxscr;

                // compute slopes

                AMREX_PARALLEL_FOR_1D(hi-lo+1, i, {

                    Real slope = 0.0;
                    int r = i + lo;

                    if (slope_order == 0) {
                        slope = 0.0;
                    } else if (slope_order == 2) {
                        if (r == 0) {
                            // one-sided difference
                            slope = s(n,r+1)-s(n,r)
                        } else if (r == nr_lev-1) {
                            // one-sided difference
                            slope = s(n,r)-s(n,r-1)
                        } else {
                            // do standard limiting on interior cells
                            Real del = 0.5*(s(n,r+1) - s(n,r-1))
                            Real dpls = 2.0*(s(n,r+1) - s(n,r  ))
                            Real dmin = 2.0*(s(n,r  ) - s(n,r-1))
                            Real slim = min(fabs(dpls), fabs(dmin))
                            slim = dpls*dmin>0.0 ? slim : 0.0;
                            Real sflag = copysign(1.0,del)
                            slope = sflag*min(slim,fabs(del))
                        }

                    } else if (slope_order == 4) {

                        // we need to calculate dxscr(fromm) for r-1 and r+1
                        Real dxscrm = 0.0;
                        Real dxscrp = 0.0;
                        // r-1
                        int rm = r-1;

                        if (rm == 0) {
                            // one-sided difference
                            dxscrm = s(n,rm+1)-s(n,rm);
                        } else if (rm == nr_lev-1) {
                            // one-sided difference
                            dxscrm = s(n,rm)-s(n,rm-1)
                        } else if (rm > 0 && rm < nr_lev-1) {
                            // do standard limiting to compute temporary slopes
                            dxscr(cen) = 0.5*(s(n,rm+1)-s(n,rm-1))
                            Real dpls = 2.0*(s(n,rm+1)-s(n,rm  ))
                            Real dmin = 2.0*(s(n,rm  )-s(n,rm-1))
                            dxscr(lim)= min(fabs(dmin),fabs(dpls));
                            dxscr(lim) = dpls*dmin>0.0 ? dxscr(lim) : 0.0;
                            dxscr(flag) = copysign(1.0,dxscr(cen));
                            dxscrm = dxscr(flag) * min(dxscr(lim),fabs(dxscr(cen)));
                        }

                        int rp = r+1;

                        if (rp == 0) {
                            // one-sided difference
                            dxscrp = s(n,rp+1)-s(n,rp)
                        } else if (rp == nr_lev-1) {
                            // one-sided difference
                            dxscrp = s(n,rp)-s(n,rp-1)
                        } else if (rp > 0 && rp < nr_lev-1) {
                            // do standard limiting to compute temporary slopes
                            dxscr(cen) = 0.5*(s(n,rp+1)-s(n,rp-1))
                            Real dpls = 2.0*(s(n,rp+1)-s(n,rp  ))
                            Real dmin = 2.0*(s(n,rp  )-s(n,rp-1))
                            dxscr(lim)= min(fabs(dmin),fabs(dpls));
                            dxscr(lim) = dpls*dmin>0.0 ? dxscr(lim) : 0.0;
                            dxscr(flag) = copysign(1.0,dxscr(cen));
                            dxscrp = dxscr(flag) * min(dxscr(lim),fabs(dxscr(cen)));
                        }
                        // end do

                        // now find dxscr for r
                        if (r > 0 && r < nr_lev-1) {
                            // do standard limiting to compute temporary slopes
                            dxscr(cen) = 0.5*(s(n,rp+1)-s(n,rp-1))
                            Real dpls = 2.0*(s(n,rp+1)-s(n,rp  ))
                            Real dmin = 2.0*(s(n,rp  )-s(n,rp-1))
                            dxscr(lim)= min(fabs(dmin),fabs(dpls));
                            dxscr(lim) = dpls*dmin>0.0 ? dxscr(lim) : 0.0;
                            dxscr(flag) = copysign(1.0,dxscr(cen));
                        } 

                        if (r == 0) {
                            // one-sided difference
                            slope = s(n,r+1)-s(n,r)
                        } else if (r == nr_lev-1) {
                            // one-sided difference
                            slope = s(n,r)-s(n,r-1)
                        } else {
                            // fourth-order limited slopes on interior
                            Real ds = 4.0/3.0*dxscr(cen) - (dxscrp + dxscrm)/6.0;
                            slope = dxscr(flag)*min(fabs(ds),dxscr(lim));
                        }

                    } // which slope order

                    // compute sedgel and sedger
                    Real u = 0.5*(w0_p(n,r)+w0_p(n,r+1))
                    Real ubardth = dth*u/dr_lev;
                    sedgel[r+1] = s(n,r) + (0.5-ubardth)*slope + dth * force(n,r)
                    sedger[r] = s(n,r) - (0.5+ubardth)*slope + dth * force(n,r)
                });

            } else if (ppm_type == 1) {

        // interpolate s to radial edges, store these temporary values into sedgel

                AMREX_PARALLEL_FOR_1D(hi-lo+1, i, {

                    int r = i + lo;

                    // calculate sm

                    // compute van Leer slopes
                    // r - 1
                    Real dsvlm = 0.0
                    int rm = r - 1;
                    if (rm == 0) {
                        // one-sided difference
                        dsvlm = s(n,rm+1)-s(n,rm)
                    } else if (rm == nr_lev-1) {
                        // one-sided difference
                        dsvlm = s(n,rm)-s(n,rm-1)
                    } else if (r > 0 && r < nr_lev-1) {
                        Real del  = 0.5 * (s(n,rm+1) - s(n,rm-1))
                        Real dmin = 2.0  * (s(n,rm  ) - s(n,rm-1))
                        Real dpls = 2.0  * (s(n,rm+1) - s(n,rm  ))
                        dsvlm = dmin*dpls > 0.0 ? copysign(1.0,del)*min(fabs(del),min(fabs(dmin),fabs(dpls))) : 0.0;
                    }

                    // r
                    Real dsvl = 0.0
                    if (r == 0) {
                        // one-sided difference
                        dsvl = s(n,r+1)-s(n,r)
                    } else if (r == nr_lev-1) {
                        // one-sided difference
                        dsvl = s(n,r)-s(n,r-1)
                    } else if (r > 0 && r < nr_lev-1) {
                        Real del  = 0.5 * (s(n,r+1) - s(n,r-1))
                        Real dmin = 2.0  * (s(n,r  ) - s(n,r-1))
                        Real dpls = 2.0  * (s(n,r+1) - s(n,r  ))
                        dsvl = dmin*dpls > 0.0 ? copysign(1.0,del)*min(fabs(del),min(fabs(dmin),fabs(dpls))) : 0.0;
                    }

                    Real sm = 0.0;
                    if (r == 0) {
                        // 2nd order interpolation to boundary face
                        sm = s(n,r) - 0.5*dsvl
                    } else if (r == nr_lev) {
                        // 2nd order interpolation to boundary face
                        sm = s(n,r-1) + 0.5*dsvl
                    } else {
                        // 4th order interpolation of s to radial faces
                        sm = 0.5*(s(n,r)+s(n,r-1)) - (dsvl-dsvlm)/6.0
                        // make sure sedgel lies in between adjacent cell-centered values
                        sm = max(sm,min(s(n,r),s(n,r-1)))
                        sm = min(sm,max(s(n,r),s(n,r-1)))
                    }

                    // calculate sp
                    // compute van Leer slopes
                    // r + 1 - 1
                    dsvlm = dsvl;

                    // r + 1
                    int rp = r + 1;
                    if (rp == 0) {
                        // one-sided difference
                        dsvl = s(n,rp+1)-s(n,rp)
                    } else if (rp == nr_lev-1) {
                        // one-sided difference
                        dsvl = s(n,rp)-s(n,rp-1)
                    } else if (rp > 0 && rp < nr_lev-1) {
                        Real del  = 0.5 * (s(n,rp+1) - s(n,rp-1))
                        Real dmin = 2.0  * (s(n,rp  ) - s(n,rp-1))
                        Real dpls = 2.0  * (s(n,rp+1) - s(n,rp  ))
                        dsvl = dmin*dpls > 0.0 ? copysign(1.0,del)*min(fabs(del),min(fabs(dmin),fabs(dpls))) : 0.0;
                    }

                    Real sp = 0.0;
                    if (rp == 0) {
                        // 2nd order interpolation to boundary face
                        sp = s(n,rp) - 0.5*dsvl(n,rp)
                    } else if (rp == nr_lev) {
                        // 2nd order interpolation to boundary face
                        sp = s(n,rp-1) + 0.5*dsvl
                    } else {
                        // 4th order interpolation of s to radial faces
                        sp = 0.5*(s(n,rp)+s(n,rp-1)) - (dsvl-dsvlm)/6.0
                        // make sure sedgel lies in between adjacent cell-centered values
                        sp = max(sp,min(s(n,rp),s(n,rp-1)))
                        sp = min(sp,max(s(n,r),s(n,rp-1)))
                    }

                    // modify using quadratic limiters
                    if ((sp-s(n,r))*(s(n,r)-sm) <= 0.0) {
                        sp = s(n,r)
                        sm = s(n,r)
                    } else if (fabs(sp-s(n,r)) >= 2.0*fabs(sm-s(n,r))) {
                        sp = 3.0*s(n,r) - 2.0*sm
                    } else if (fabs(sm-s(n,r)) >= 2.0*fabs(sp-s(n,r))) {
                        sm = 3.0*s(n,r) - 2.0*sp
                    }
            //   end do

                    // compute Ip and Im
                    Real sigmap = fabs(w0_p(n,r+1))*dtdr
                    Real sigmam = fabs(w0_p(n,r  ))*dtdr
                    Real s6 = 6.0*s(n,r) - 3.0*(sm+sp)
                    Real Ip = 0.0;
                    Real Im = 0.0;
                    if (w0_p(n,r+1) > rel_eps) {
                        Ip = sp - (sigmap/2.0)*(sp-sm-(1.0-2.0/3.0*sigmap)*s6)
                    } else {
                        Ip = s(n,r)
                    }
                    if (w0_p(n,r) < -rel_eps) {
                        Im = sm + (sigmam/2.0)*(sp-sm+(1.0-2.0/3.0*sigmam)*s6)
                    } else {
                        Im = s(n,r)
                    }

                    // compute sedgel and sedger
                    sedgel[r+1] = Ip + dth * force_p(n,r)
                    sedger[r] = Im + dth * force_p(n,r)
              
                });

            } else if (ppm_type == 2) {

        // interpolate s to radial edges

        // store centered differences in dsvl
        //       do r=lo-3,hi+3

                AMREX_PARALLEL_FOR_1D(hi-lo+2, i, {
                    int r = i + lo;

                    // left side 
                    Real dsvl = 0.0;
                    if (r == 0) {
                        // one-sided difference
                        dsvl = s(n,r)-s(n,r-1)
                    } else if (r == nr_lev-1) {
                        // one-sided difference
                        dsvl = s(n,r-1)-s(n,r-2)
                    } else if (r > 0 && r < nr_lev-1) {
                        // centered difference
                        dsvl = 0.5 * (s(n,r) - s(n,r-2))
                    }

                    // right side 
                    Real dsvr = 0.0;
                    if (r == 0) {
                        // one-sided difference
                        dsvr = s(n,r+1)-s(n,r)
                    } else if (r == nr_lev-1) {
                        // one-sided difference
                        dsvr = s(n,r)-s(n,r-1)
                    } else if (r > 0 && r < nr_lev-1) {
                        // centered difference
                        dsvr = 0.5 * (s(n,r+1) - s(n,r-1))
                    }

            //   end do

            //   do r=lo-2,hi+3
                    if (r == 0) {
                        // 2nd order interpolation to boundary face
                        sedgel(n,r) = s(n,r) - 0.5*dsvr
                    } else if (r == nr_lev) {
                        // 2nd order interpolation to boundary face
                        sedgel(n,r) = s(n,r-1) + 0.5*dsvr
                    } else if (r > 0 && r < nr_lev) {
                        // 4th order interpolation of s to radial faces
                        sedgel(n,r) = 0.5*(s(n,r)+s(n,r-1)) - (dsvr-dsvl)/6.0
                        if (r >= 2 && r <= nr_lev-2) {
                            // limit sedge
                            if ((sedgel(n,r)-s(n,r-1))*(s(n,r)-sedgel(n,r)) < 0.0) {
                                Real D2  = 3.0*(s(n,r-1)-2.0*sedgel(n,r)+s(n,r))
                                Real D2L = s(n,r-2)-2.0*s(n,r-1)+s(n,r)
                                Real D2R = s(n,r-1)-2.0*s(n,r)+s(n,r+1)
                                Real sgn = copysign(1.0,D2);
                                Real D2LIM = sgn*max(min(C*sgn*D2L,min(C*sgn*D2R,sgn*D2)),0.0);
                                sedgel(n,r) = 0.5*(s(n,r-1)+s(n,r)) - D2LIM/6.0
                            }
                        }
                    }
        

            //   do r=lo,hi

                // use Colella 2008 limiters
                // This is a new version of the algorithm
                // to eliminate sensitivity to roundoff.

                    Real sm = 0.0;
                    Real sp = 0.0;

                    if (r >= 2 && r <= nr_lev-3) {

                        Real alphap = sedgel(n,r+1)-s(n,r)
                        Real alpham = sedgel(n,r  )-s(n,r)
                        bool bigp = fabs(alphap)>2.0*fabs(alpham)
                        bool bigm = fabs(alpham)>2.0*fabs(alphap)
                        bool extremum = false

                        if (alpham*alphap >= 0.0) {
                            extremum = true
                        } else if (bigp || bigm) {
                            // Possible extremum. We look at cell centered values and face
                            // centered values for a change in copysign in the differences adjacent to
                            // the cell. We use the pair of differences whose minimum magnitude is
                            // the largest, and thus least susceptible to sensitivity to roundoff.
                            Real dafacem = sedgel(n,r) - sedgel(n,r-1)
                            Real dafacep = sedgel(n,r+2) - sedgel(n,r+1)
                            Real dabarm = s(n,r) - s(n,r-1)
                            Real dabarp = s(n,r+1) - s(n,r)
                            Real dafacemin = min(fabs(dafacem),fabs(dafacep))
                            Real dabarmin= min(fabs(dabarm),fabs(dabarp))
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
                            Real D2L = s(n,r-2)-2.0*s(n,r-1)+s(n,r)
                            Real D2R = s(n,r)-2.0*s(n,r+1)+s(n,r+2)
                            Real D2C = s(n,r-1)-2.0*s(n,r)+s(n,r+1)
                            Real sgn = copysign(1.0,D2);
                            Real D2LIM = max(min(sgn*D2,min(C*sgn*D2L,min(C*sgn*D2R,C*sgn*D2C))),0.0);
                            Real D2ABS = max(fabs(D2),1.e-10);
                            alpham = alpham*D2LIM/D2ABS;
                            alphap = alphap*D2LIM/D2ABS;
                        } else {
                            if (bigp) {
                                Real sgn = copysign(1.0,alpham);
                                Real amax = -alphap*alphap / (4.0*(alpham + alphap));
                                Real delam = s(n,r-1) - s(n,r)
                                if (sgn*amax >= sgn*delam) {
                                    if (sgn*(delam - alpham)>=1.e-10) {
                                        alphap = (-2.0*delam - 2.0*sgn*sqrt(delam*delam - delam*alpham));
                                    } else {
                                        alphap = -2.0*alpham;
                                    }
                                }
                            }
                            if (bigm) {
                                Real sgn = copysign(1.0,alphap);
                                Real amax = -alpham*alpham / (4.0*(alpham + alphap));
                                Real delap = s(n,r+1) - s(n,r)
                                if (sgn*amax >= sgn*delap) {
                                    if (sgn*(delap - alphap)>=1.e-10) {
                                        alpham = (-2.0*delap - 2.0*sgn*sqrt(delap*delap - delap*alphap));
                                    } else {
                                        alpham = -2.0*alphap;
                                    }
                                }
                            }
                        }

                        sm = s(n,r) + alpham
                        sp = s(n,r) + alphap

                    } else {

                        sp = sedgel(n,r+1)
                        sm = sedgel(n,r  )

                    } // test (r >= 2 && r <= nr_lev-3)

                    // compute Ip and Im
                    Real sigmap = fabs(w0_p(n,r+1))*dtdr
                    Real sigmam = fabs(w0_p(n,r  ))*dtdr
                    Real s6 = 6.0*s(n,r) - 3.0*(sm+sp)
                    Real Ip = 0.0;
                    Real Im = 0.0;
                    if (w0_p(n,r+1) > rel_eps) {
                        Ip = sp - (sigmap/2.0)*(sp-sm-(1.0-2.0/3.0*sigmap)*s6)
                    } else {
                        Ip = s(n,r)
                    }
                    if (w0_p(n,r) < -rel_eps) {
                        Im = sm + (sigmam/2.0)*(sp-sm+(1.0-2.0/3.0*sigmam)*s6)
                    } else {
                        Im = s(n,r)
                    }

                    // compute sedgel and sedger
                    sedgel(n,r+1) = Ip + dth * force(n,r)
                    sedger(n,r  ) = Im + dth * force(n,r)

                }); // loop over r
            }
        }
    }

    for (int n = 0; n <= max_radial_level; ++n) {

        Real * AMREX_RESTRICT sedgel = sedgel_vec[n].dataPtr();
        Real * AMREX_RESTRICT sedger = sedger_vec[n].dataPtr(); 

        for (int i = 0; i < numdisjointchunks[n]; ++i) {

     // sync up edge states at coarse-fine interface

            // if we are not at the finest level, copy in the sedger and sedgel states
            // from the next finer level at the c-f interface
            if (n != finest_radial_level) {
                sedger(n,r_start_coord(n+1,i)/2) = sedger(n+1,r_start_coord(n+1,i))
                sedgel(n,(r_end_coord(n+1,i)+1)/2) = sedgel(n+1,r_end_coord(n+1,i)+1)
            }

            // if we are not at the coarsest level, copy in the sedgel and sedger states
            // from the next coarser level at the c-f interface
            if (n != 0) {
                sedgel(n,lo) = sedgel(n-1,lo/2)
                sedger(n,hi+1) = sedger(n-1,(hi+1)/2)
            }


     // solve Riemann problem to get final edge state

    //        do r=lo,hi+1

            AMREX_PARALLEL_FOR_1D(hi-lo+2, i, {
                int r = i + lo;

                if (r == 0) {
                    // pick interior state at lo domain boundary
                    sedge[r] = sedger[r]
                } else if (r == nr_lev) {
                    // pick interior state at hi domain boundary
                    sedge[r] = sedgel[r]
                } else {
                    // upwind
                    sedge[r] = w0_p(n,r) > 0.0 ? sedgel[r] : sedger[r]
                    sedge[r] = fabs(w0_p(n,r)) < rel_eps ? 0.5*(sedger[r] + sedgel[r]) : sedge[r]
                }
            });
        //    end do

        }  // loop over disjointchunks
    } // loop over levels
    
}
#include <Maestro.H>
#include <MaestroHydro_F.H>

using namespace amrex;

#if (AMREX_SPACEDIM == 2)

void 
Maestro::VelPredInterface(const MFIter& mfi, 
                          Array4<const Real> const utilde,
                          Array4<const Real> const ufull,
                          Array4<const Real> const utrans,
                          Array4<const Real> const vtrans,
                          Array4<const Real> const Imu,
                          Array4<const Real> const Ipu,
                          Array4<const Real> const Imv,
                          Array4<const Real> const Ipv,
                          Array4<Real> const ulx,
                          Array4<Real> const urx,
                          Array4<Real> const uimhx,
                          Array4<Real> const uly,
                          Array4<Real> const ury,
                          Array4<Real> const uimhy,
                          const Box& domainBox,
                          const Real* dx)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::VelPredInterface()",VelPredInterface);

    // NOTE: for ppm_type_local == 0, slopex == Ipu, slopey == Imv

    ////////////////////////////////////
    // Create u_{\i-\half\e_x}^x, etc.
    ////////////////////////////////////

    const Box& tileBox = mfi.tilebox();
    const Box& obx = amrex::grow(tileBox, 1);
    const Box& mxbx = amrex::growLo(obx, 0, -1);
    const Box& mybx = amrex::growLo(obx, 1, -1);

    Real rel_eps;
    get_rel_eps(&rel_eps);

    Real dt2 = 0.5 * dt;

    Real hx = dx[0];
    Real hy = dx[1];

    // x-direction
    int ilo = domainBox.loVect()[0];
    int ihi = domainBox.hiVect()[0];

    int bclo = phys_bc[0];
    int bchi = phys_bc[AMREX_SPACEDIM];

    int ppm_type_local = ppm_type;

    AMREX_PARALLEL_FOR_3D(mxbx, i, j, k, 
    {
        if (ppm_type_local == 0) {
            Real maxu = max(0.0,ufull(i-1,j,k,0));
            Real minu = min(0.0,ufull(i  ,j,k,0));
            // extrapolate both components of velocity to left face
            ulx(i,j,k,0) = utilde(i-1,j,k,0) + (0.5 - (dt2/hx)*maxu)*Ipu(i-1,j,k,0);
            ulx(i,j,k,1) = utilde(i-1,j,k,1) + (0.5 - (dt2/hx)*maxu)*Ipu(i-1,j,k,1);
            // extrapolate both components of velocity to right face
            urx(i,j,k,0) = utilde(i  ,j,k,0) - (0.5 + (dt2/hx)*minu)*Ipu(i  ,j,k,0);
            urx(i,j,k,1) = utilde(i  ,j,k,1) - (0.5 + (dt2/hx)*minu)*Ipu(i  ,j,k,1);
        } else if (ppm_type_local == 1 || ppm_type_local == 2) {
            // extrapolate both components of velocity to left face
            ulx(i,j,k,0) = Ipu(i-1,j,k,0);
            ulx(i,j,k,1) = Ipv(i-1,j,k,0);
            // extrapolate both components of velocity to right face
            urx(i,j,k,0) = Imu(i,j,k,0);
            urx(i,j,k,1) = Imv(i,j,k,0);
        }

        // impose lo side bc's
        if (i == ilo) {
            switch(bclo) {
                case Inflow:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        ulx(i,j,k,n) = utilde(i-1,j,k,n);
                        urx(i,j,k,n) = utilde(i-1,j,k,n);
                    }
                case SlipWall:
                case Symmetry: 
                    ulx(i,j,k,0) = 0.0;
                    urx(i,j,k,0) = 0.0;
                    ulx(i,j,k,1) = urx(i,j,k,1);
                case NoSlipWall:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        ulx(i,j,k,n) = 0.0;
                        urx(i,j,k,n) = 0.0;
                    }
                case Outflow:
                    urx(i,j,k,0) = min(urx(i,j,k,0),0.0);
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        urx(i,j,k,n) = ulx(i,j,k,n);
                    }
                case Interior:
                    break;
                default:
                    Print() << "velpred_2d: invalid boundary type phys_bc(1,1)" << std::endl;
            }

        // impose hi side bc's
        } else if (i == ihi+1) {
            switch (bchi) {
                case Inflow:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        ulx(i,j,k,n) = utilde(i,j,k,n);
                        urx(i,j,k,n) = utilde(i,j,k,n);
                    }
                case SlipWall:
                case Symmetry:
                    ulx(i,j,k,0) = 0.0;
                    urx(i,j,k,0) = 0.0;
                    urx(i,j,k,1) = ulx(i,j,k,1);
                case NoSlipWall:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        ulx(i,j,k,n) = 0.0;
                        urx(i,j,k,n) = 0.0;
                    }
                case Outflow:
                    ulx(i,j,k,0) = max(ulx(i,j,k,0),0.0);
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        urx(i,j,k,n) = ulx(i,j,k,n);
                    }
                case Interior:
                    break;
                default:
                    Print() << "velpred_2d: invalid boundary type phys_bc(1,2)" << std::endl;
            }
        }

        // No need to compute uimh(:,:,0) since it's equal to utrans-w0
        // upwind using full velocity to get transverse component of uimhx
        // Note: utrans already contains w0
        uimhx(i,j,k,1) = utrans(i,j,k) > 0.0 ? ulx(i,j,k,1) : urx(i,j,k,1);
        uimhx(i,j,k,1) = fabs(utrans(i,j,k)) < rel_eps ? 
            0.5*(ulx(i,j,k,1)+urx(i,j,k,1)) : uimhx(i,j,k,1);
    });

    // y-direction
    int jlo = domainBox.loVect()[1];
    int jhi = domainBox.hiVect()[1];
    bclo = phys_bc[1];
    bchi = phys_bc[AMREX_SPACEDIM+1];

    AMREX_PARALLEL_FOR_3D(mybx, i, j, k, 
    {
        if (ppm_type_local == 0) {
            Real maxu = max(0.0,ufull(i,j-1,k,1));
            Real minu = min(0.0,ufull(i,j,k,1));
            // extrapolate both components of velocity to left face
            uly(i,j,k,0) = utilde(i,j-1,k,0) + (0.5-(dt2/hy)*maxu)*Imv(i,j-1,k,0);
            uly(i,j,k,1) = utilde(i,j-1,k,1) + (0.5-(dt2/hy)*maxu)*Imv(i,j-1,k,1);
            // extrapolate both components of velocity to right face
            ury(i,j,k,0) = utilde(i,j,k,0) - (0.5+(dt2/hy)*minu)*Imv(i,j,k,0);
            ury(i,j,k,1) = utilde(i,j,k,1) - (0.5+(dt2/hy)*minu)*Imv(i,j,k,1);
        } else if (ppm_type_local == 1 || ppm_type_local == 2) {
            // extrapolate both components of velocity to left face
            uly(i,j,k,0) = Ipu(i,j-1,k,0);
            uly(i,j,k,1) = Ipv(i,j-1,k,1);
            // extrapolate both components of velocity to right face
            ury(i,j,k,0) = Imu(i,j,k,0);
            ury(i,j,k,1) = Imv(i,j,k,1);
        }

        // impose lo side bc's
        if (j == jlo) {
            switch (bclo) {
                case Inflow:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        uly(i,j,k,n) = utilde(i,j-1,k,n);
                        ury(i,j,k,n) = utilde(i,j-1,k,n);
                    }
                case SlipWall:
                case Symmetry:
                    uly(i,j,k,0) = ury(i,j,k,0);
                    uly(i,j,k,1) = 0.0;
                    ury(i,j,k,1) = 0.0;
                case NoSlipWall:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        uly(i,j,k,n) = 0.0;
                        ury(i,j,k,n) = 0.0;
                    }
                case Outflow:
                    ury(i,j,k,1) = min(ury(i,j,k,1),0.0);
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        uly(i,j,k,n) = ury(i,j,k,n);
                    }
                case Interior:
                    break;
                default:
                    Print() << "velpred_2d: invalid boundary type phys_bc(2,1)" << std::endl;
            }

        // impose hi side bc's
        } else if (j == jhi+1) {
            switch (bchi) {
                case Inflow:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        uly(i,j,k,n) = utilde(i,j,k,n);
                        ury(i,j,k,n) = utilde(i,j,k,n);
                    }
                case SlipWall:
                case Symmetry:
                    ury(i,j,k,0) = uly(i,j,k,0);
                    uly(i,j,k,1) = 0.0;
                    ury(i,j,k,1) = 0.0;
                case NoSlipWall:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        uly(i,j,k,n) = 0.0;
                        ury(i,j,k,n) = 0.0;
                    }
                case Outflow:
                    uly(i,j,k,1) = max(uly(i,j,k,1),0.0);
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        ury(i,j,k,n) = uly(i,j,k,n);
                    }
                case Interior:
                    break;
                default:
                    Print() << "velpred_2d: invalid boundary type phys_bc(2,2)" << std::endl;
            }
        }
        // No need to compute uimh(:,:,1) since it's equal to utrans-w0
        // upwind using full velocity to get transverse component of uimhy
        // Note: utrans already contains w0
        uimhy(i,j,k,0) = vtrans(i,j,k) > 0.0 ? uly(i,j,k,0) : ury(i,j,k,0);
        uimhy(i,j,k,0) = fabs(vtrans(i,j,k)) < rel_eps ? 
            0.5*(uly(i,j,k,0)+ury(i,j,k,0)) : uimhy(i,j,k,0);
    });
}

void 
Maestro::VelPredVelocities(const MFIter& mfi, 
                            Array4<const Real> const utilde,
                            Array4<const Real> const utrans,
                            Array4<const Real> const vtrans,
                            Array4<Real> const umac,
                            Array4<Real> const vmac,
                            Array4<const Real> const Imfx,
                            Array4<const Real> const Ipfx,
                            Array4<const Real> const Imfy,
                            Array4<const Real> const Ipfy,
                            Array4<const Real> const ulx,
                            Array4<const Real> const urx,
                            Array4<const Real> const uimhx,
                            Array4<const Real> const uly,
                            Array4<const Real> const ury,
                            Array4<const Real> const uimhy,
                            Array4<const Real> const force,
                            Array4<const Real> const w0_cart,
                            const Box& domainBox,
                            const Real* dx,
                            const Vector<BCRec>& bcs)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::VelPredVelocities()",VelPredVelocities);

    //******************************************************************
    // Create umac and vmac
    //******************************************************************

    const Box& xbx = mfi.nodaltilebox(0);
    const Box& ybx = mfi.nodaltilebox(1);

    Real rel_eps;
    get_rel_eps(&rel_eps);

    Real dt2 = 0.5 * dt;
    Real dt4 = 0.25 * dt;

    Real hx = dx[0];
    Real hy = dx[1];

    int ppm_trace_forces_local = ppm_trace_forces;

    // x-direction
    int ilo = domainBox.loVect()[0];
    int ihi = domainBox.hiVect()[0];

    int bclo = phys_bc[0];
    int bchi = phys_bc[AMREX_SPACEDIM];

    AMREX_PARALLEL_FOR_3D(xbx, i, j, k, 
    {
        // use the traced force if ppm_trace_forces = 1
        Real fl = ppm_trace_forces_local == 0 ? force(i-1,j,k,0) : Ipfx(i-1,j,k,0);
        Real fr = ppm_trace_forces_local == 0 ? force(i,j,k,0) : Imfx(i,j,k,0);

        // extrapolate to edges
        Real umacl = ulx(i,j,k,0) 
            - (dt4/hy)*(vtrans(i-1,j+1,k)+vtrans(i-1,j,k)) 
            * (uimhy(i-1,j+1,k,0)-uimhy(i-1,j,k,0)) + dt2*fl;
        Real umacr = urx(i,j,k,0) 
            - (dt4/hy)*(vtrans(i  ,j+1,k)+vtrans(i  ,j,k)) 
            * (uimhy(i  ,j+1,k,0)-uimhy(i  ,j,k,0)) + dt2*fr;

        // solve Riemann problem using full velocity
        umac(i,j,k) = 0.5*(umacl+umacr) > 0.0 ? umacl : umacr; 
        umac(i,j,k) = (umacl <= 0.0 && umacr >= 0.0) || (fabs(umacl+umacr) < rel_eps) ? 
            0.0 : umac(i,j,k);

        // impose lo side bc's
        if (i == ilo) {
            switch (bclo) {
                case Inflow:
                    umac(i,j,k) = utilde(i-1,j,k,0);
                case SlipWall:
                case NoSlipWall:
                case Symmetry:
                    umac(i,j,k) = 0.0;
                case Outflow:
                    umac(i,j,k) = min(umacr,0.0);
                case Interior:
                    break;
                default:
                    Print() << "velpred_2d: invalid boundary type phys_bc(1,1)" << std::endl;
            }

        // impose hi side bc's
        } else if (i == ihi+1) {
            switch (bchi) {
                case Inflow:
                    umac(i,j,k) = utilde(i,j,k,0);
                case SlipWall:
                case NoSlipWall:
                case Symmetry:
                    umac(i,j,k) = 0.0;
                case Outflow:
                    umac(i,j,k) = max(umacl,0.0);
                case Interior:
                default:
                    Print() << "velpred_2d: invalid boundary type phys_bc(1,2)" << std::endl;
            }
        }
    });

    // y-direction
    int jlo = domainBox.loVect()[1];
    int jhi = domainBox.hiVect()[1];
    bclo = phys_bc[1];
    bchi = phys_bc[AMREX_SPACEDIM+1];

    AMREX_PARALLEL_FOR_3D(ybx, i, j, k, 
    {
        // use the traced force if ppm_trace_forces = 1
        Real fl = ppm_trace_forces_local == 0 ? force(i,j-1,k,1) : Ipfy(i,j-1,k,1);
        Real fr = ppm_trace_forces_local == 0 ? force(i,j,k,1) : Imfy(i,j,k,1);

        // extrapolate to edges
        Real vmacl = uly(i,j,k,1) 
            - (dt4/hx)*(utrans(i+1,j-1,k)+utrans(i,j-1,k)) 
            * (uimhx(i+1,j-1,k,1)-uimhx(i,j-1,k,1)) + dt2*fl;
        Real vmacr = ury(i,j,k,1) 
            - (dt4/hx)*(utrans(i+1,j,k)+utrans(i,j,k)) 
            * (uimhx(i+1,j,k,1)-uimhx(i,j,k,1)) + dt2*fr;

        // solve Riemann problem using full velocity
        bool test = (vmacl+w0_cart(i,j,k,AMREX_SPACEDIM-1) <= 0.0 && 
                     vmacr+w0_cart(i,j,k,AMREX_SPACEDIM-1) >= 0.0) || 
                    (fabs(vmacl+vmacr+2*w0_cart(i,j,k,AMREX_SPACEDIM-1)) < rel_eps);
        vmac(i,j,k) = 0.5*(vmacl+vmacr)+w0_cart(i,j,k,AMREX_SPACEDIM-1) > 0.0 ? 
                      vmacl : vmacr;
        vmac(i,j,k) = test ? 0.0 : vmac(i,j,k);

        // impose lo side bc's
        if (j == jlo) {
            switch (bclo) {
                case Inflow:
                    vmac(i,j,k) = utilde(i,j-1,k,1);
                case SlipWall:
                case NoSlipWall:
                case Symmetry:
                    vmac(i,j,k) = 0.0;
                case Outflow:
                    vmac(i,j,k) = min(vmacr,0.0);
                case Interior:
                    break;
                default:
                    Print() << "velpred_2d: invalid boundary type phys_bc(2,1)" << std::endl;
            }

        // impose hi side bc's
        } else if (j == jhi+1) {
            switch (bchi) {
                case Inflow:
                    vmac(i,j,k) = utilde(i,j,k,1);
                case SlipWall:
                case NoSlipWall:
                case Symmetry:
                    vmac(i,j,k) = 0.0;
                case Outflow:
                    vmac(i,j,k) = max(vmacl,0.0);
                case Interior:
                    break;
                default:
                    Print() << "velpred_2d: invalid boundary type phys_bc(2,2)" << std::endl;
            }
        }
    });
}

#else 

void 
Maestro::VelPredInterface(const MFIter& mfi, 
                          Array4<const Real> const utilde,
                          Array4<const Real> const ufull,
                          Array4<const Real> const utrans,
                          Array4<const Real> const vtrans,
                          Array4<const Real> const wtrans,
                          Array4<const Real> const Imu,
                          Array4<const Real> const Ipu,
                          Array4<const Real> const Imv,
                          Array4<const Real> const Ipv,
                          Array4<const Real> const Imw,
                          Array4<const Real> const Ipw,
                          Array4<Real> const ulx,
                          Array4<Real> const urx,
                          Array4<Real> const uimhx,
                          Array4<Real> const uly,
                          Array4<Real> const ury,
                          Array4<Real> const uimhy,
                          Array4<Real> const ulz,
                          Array4<Real> const urz,
                          Array4<Real> const uimhz,
                          const Box& domainBox,
                          const Real* dx)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::VelPredInterface()",VelPredInterface);

    ////////////////////////////////////
    // Create u_{\i-\half\e_x}^x, etc.
    ////////////////////////////////////

    // NOTE: for ppm_type_local == 0, slopex == Ipu, slopey == Imv

    // normal predictor states
     // Allocated from lo:hi+1 in the normal direction
     // lo-1:hi+1 in the transverse directions

    const Box& tileBox = mfi.tilebox();
    const Box& obx = amrex::grow(tileBox, 1);
    const Box& mxbx = amrex::growLo(obx, 0, -1);
    const Box& mybx = amrex::growLo(obx, 1, -1);
    const Box& mzbx = amrex::growLo(obx,2, -1);

    Real rel_eps;
    get_rel_eps(&rel_eps);

    Real dt2 = 0.5 * dt;

    Real hx = dx[0];
    Real hy = dx[1];
    Real hz = dx[2];

    // x-direction
    int ilo = domainBox.loVect()[0];
    int ihi = domainBox.hiVect()[0];

    int bclo = phys_bc[0];
    int bchi = phys_bc[AMREX_SPACEDIM];

    int ppm_type_local = ppm_type;

    AMREX_PARALLEL_FOR_3D(mxbx, i, j, k, 
    {
        if (ppm_type_local == 0) {
            Real maxu = (0.5 - dt2*max(0.0,ufull(i-1,j,k,0))/hx);
            Real minu = (0.5 + dt2*min(0.0,ufull(i  ,j,k,0))/hx);

            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                // extrapolate all components of velocity to left face
                ulx(i,j,k,n) = utilde(i-1,j,k,n) + maxu * Ipu(i-1,j,k,n);
                // extrapolate all components of velocity to right face
                urx(i,j,k,n) = utilde(i,j,k,n) - minu * Ipu(i,j,k,n);
            }
        } else if (ppm_type_local == 1 || ppm_type_local == 2) {
            // extrapolate all components of velocity to left face
            ulx(i,j,k,0) = Ipu(i-1,j,k,0);
            ulx(i,j,k,1) = Ipv(i-1,j,k,0);
            ulx(i,j,k,2) = Ipw(i-1,j,k,0);

            // extrapolate all components of velocity to right face
            urx(i,j,k,0) = Imu(i,j,k,0);
            urx(i,j,k,1) = Imv(i,j,k,0);
            urx(i,j,k,2) = Imw(i,j,k,0);
        }

        // impose lo side bc's
        if (i == ilo) {
            switch (bclo) {
                case Inflow:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        ulx(i,j,k,n) = utilde(i-1,j,k,n);
                        urx(i,j,k,n) = utilde(i-1,j,k,n);
                    }
                case SlipWall:
                case Symmetry:
                    ulx(i,j,k,0) = 0.0;
                    urx(i,j,k,0) = 0.0;
                    ulx(i,j,k,1) = urx(i,j,k,1);
                    ulx(i,j,k,2) = urx(i,j,k,2);
                case NoSlipWall:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        ulx(i,j,k,n) = 0.0;
                        urx(i,j,k,n) = 0.0;
                    }
                case Outflow:
                    urx(i,j,k,0) = min(urx(i,j,k,0),0.0);
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        ulx(i,j,k,n) = ulx(i,j,k,n);
                    }
                case Interior:
                    break;
                default:
                    Print() << "velpred_3d: invalid boundary type phys_bc(1,1)" << std::endl;
            }

        // impose hi side bc's
        } else if (i == ihi+1) {
            switch (bchi) {
                case Inflow:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        ulx(i,j,k,n) = utilde(i,j,k,n);
                        urx(i,j,k,n) = utilde(i,j,k,n);
                    }
                case SlipWall:
                case Symmetry:
                    ulx(i,j,k,0) = 0.0;
                    urx(i,j,k,0) = 0.0;
                    urx(i,j,k,1) = ulx(i,j,k,1);
                    urx(i,j,k,2) = ulx(i,j,k,2);
                case NoSlipWall:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        ulx(i,j,k,n) = 0.0;
                        urx(i,j,k,n) = 0.0;
                    }
                case Outflow:
                    ulx(i,j,k,0) = max(ulx(i,j,k,0),0.0);
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        urx(i,j,k,n) = ulx(i,j,k,n);
                    }
                case Interior:
                    break;
                default:
                    Print() << "velpred_3d: invalid boundary type phys_bc(1,2)" << std::endl;
            }
        }

        // No need to compute uimhx(:,:,:,0) since it's equal to utrans-w0
        // upwind using full velocity to get transverse components of uimhx
        // Note: utrans already contains w0
        uimhx(i,j,k,1) = utrans(i,j,k) > 0.0 ? ulx(i,j,k,1) : urx(i,j,k,1);
        uimhx(i,j,k,1) = fabs(utrans(i,j,k)) < rel_eps ? 
            0.5*(ulx(i,j,k,1)+urx(i,j,k,1)) : uimhx(i,j,k,1);

        uimhx(i,j,k,2) = utrans(i,j,k) > 0.0 ? ulx(i,j,k,2) : urx(i,j,k,2);
        uimhx(i,j,k,2) = fabs(utrans(i,j,k)) < rel_eps ? 
            0.5*(ulx(i,j,k,2)+urx(i,j,k,2)) : uimhx(i,j,k,2);
    });

    // y-direction
    int jlo = domainBox.loVect()[1];
    int jhi = domainBox.hiVect()[1];
    bclo = phys_bc[1];
    bchi = phys_bc[AMREX_SPACEDIM+1];

    AMREX_PARALLEL_FOR_3D(mybx, i, j, k, 
    {
        if (ppm_type_local == 0) {
            Real maxu = (0.5 - dt2*max(0.0,ufull(i,j-1,k,1))/hy);
            Real minu = (0.5 + dt2*min(0.0,ufull(i,j  ,k,1))/hy);

            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                // extrapolate all components of velocity to left face
                uly(i,j,k,n) = utilde(i,j-1,k,n) + maxu * Imv(i,j-1,k,n);
                // extrapolate all components of velocity to right face
                ury(i,j,k,n) = utilde(i,j,k,n) - minu * Imv(i,j,k,n);
            }

        } else if (ppm_type_local == 1 || ppm_type_local == 2) {
            // extrapolate all components of velocity to left face
            uly(i,j,k,0) = Ipu(i,j-1,k,1);
            uly(i,j,k,1) = Ipv(i,j-1,k,1);
            uly(i,j,k,2) = Ipw(i,j-1,k,1);

            // extrapolate all components of velocity to right face
            ury(i,j,k,0) = Imu(i,j,k,1);
            ury(i,j,k,1) = Imv(i,j,k,1);
            ury(i,j,k,2) = Imw(i,j,k,1);
        }

        // impose lo side bc's
        if (j == jlo) {
            switch (bclo) {
                case Inflow:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        uly(i,j,k,n) = utilde(i,j-1,k,n);
                        ury(i,j,k,n) = utilde(i,j-1,k,n);
                    }
                case SlipWall:
                case Symmetry:
                    uly(i,j,k,0) = ury(i,j,k,0);
                    uly(i,j,k,1) = 0.0;
                    ury(i,j,k,1) = 0.0;
                    uly(i,j,k,2) = ury(i,j,k,2);
                case NoSlipWall:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        uly(i,j,k,n) = 0.0;
                        ury(i,j,k,n) = 0.0;
                    }
                case Outflow:
                    ury(i,j,k,1) = min(ury(i,j,k,1),0.0);
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        uly(i,j,k,n) = ury(i,j,k,n);
                    }
                case Interior:
                    break;
                default:
                    Print() << "velpred_3d: invalid boundary type phys_bc(2,1)" << std::endl;
            }

        // impose hi side bc's
        } else if (j == jhi+1) {
            switch (bchi) {
                case Inflow:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        uly(i,j,k,n) = utilde(i,j,k,n);
                        ury(i,j,k,n) = utilde(i,j,k,n);
                    }
                case SlipWall:
                case Symmetry:
                    ury(i,j,k,0) = uly(i,j,k,0);
                    uly(i,j,k,1) = 0.0;
                    ury(i,j,k,1) = 0.0;
                    ury(i,j,k,2) = uly(i,j,k,2);
                case NoSlipWall:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        uly(i,j,k,n) = 0.0;
                        ury(i,j,k,n) = 0.0;
                    }
                case Outflow:
                    uly(i,j,k,1) = max(uly(i,j,k,1),0.0);
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        ury(i,j,k,n) = uly(i,j,k,n);
                    }
                case Interior:
                        break;
                default:
                    Print() << "velpred_3d: invalid boundary type phys_bc(2,2)" << std::endl;
            }
        }

        // No need to compute uimhy(:,:,:,1) since it's equal to vtrans-w0
        // upwind using full velocity to get transverse components of uimhy
        // Note: vtrans already contains w0
        uimhy(i,j,k,0) = vtrans(i,j,k) > 0.0 ? uly(i,j,k,0) : ury(i,j,k,0);
        uimhy(i,j,k,0) = fabs(vtrans(i,j,k)) < rel_eps ? 
            0.5*(uly(i,j,k,0)+ury(i,j,k,0)) : uimhy(i,j,k,0);

        uimhy(i,j,k,2) = vtrans(i,j,k) > 0.0 ? uly(i,j,k,2) : ury(i,j,k,2);
        uimhy(i,j,k,2) = fabs(vtrans(i,j,k)) < rel_eps ? 
            0.5*(uly(i,j,k,2)+ury(i,j,k,2)) : uimhy(i,j,k,2);
    });

    // z-direction
    int klo = domainBox.loVect()[2];
    int khi = domainBox.hiVect()[2];
    bclo = phys_bc[2];
    bchi = phys_bc[AMREX_SPACEDIM+2];

    AMREX_PARALLEL_FOR_3D(mzbx, i, j, k, 
    {
        if (ppm_type_local == 0) {
            Real maxu = (0.5 - dt2*max(0.0,ufull(i,j,k-1,2))/hz);
            Real minu = (0.5 + dt2*min(0.0,ufull(i,j,k  ,2))/hz);

            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                // extrapolate all components of velocity to left face
                ulz(i,j,k,n) = utilde(i,j,k-1,n) + maxu * Imw(i,j,k-1,n);
                // extrapolate all components of velocity to right face
                urz(i,j,k,n) = utilde(i,j,k,n) - minu * Imw(i,j,k,n);
            }

        } else if (ppm_type_local == 1 || ppm_type_local == 2) {
            // extrapolate all components of velocity to left face
            ulz(i,j,k,0) = Ipu(i,j,k-1,2);
            ulz(i,j,k,1) = Ipv(i,j,k-1,2);
            ulz(i,j,k,2) = Ipw(i,j,k-1,2);

            // extrapolate all components of velocity to right face
            urz(i,j,k,0) = Imu(i,j,k,2);
            urz(i,j,k,1) = Imv(i,j,k,2);
            urz(i,j,k,2) = Imw(i,j,k,2);
        }

        // impose lo side bc's
        if (k == klo) {
            switch (bclo) {
                case Inflow:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        ulz(i,j,k,n) = utilde(i,j,k-1,n);
                        urz(i,j,k,n) = utilde(i,j,k-1,n);
                    }
                case SlipWall:
                case Symmetry:
                    ulz(i,j,k,0) = urz(i,j,k,0);
                    ulz(i,j,k,1) = urz(i,j,k,1);
                    ulz(i,j,k,2) = 0.0;
                    urz(i,j,k,2) = 0.0;
                case NoSlipWall:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        ulz(i,j,k,n) = 0.0;
                        urz(i,j,k,n) = 0.0;
                    }
                case Outflow:
                    urz(i,j,k,2) = min(urz(i,j,k,2),0.0);
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        ulz(i,j,k,n) = urz(i,j,k,n);
                    }
                case Interior:
                        break;
                default:
                    Print() << "velpred_3d: invalid boundary type phys_bc(3,1)" << std::endl;
            }

        // impose hi side bc's
        } else if (k == khi+1) {
            switch (bchi) {
                case Inflow:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        ulz(i,j,k,n) = utilde(i,j,k,n);
                        urz(i,j,k,n) = utilde(i,j,k,n);
                    }
                case SlipWall:
                case Symmetry:
                    urz(i,j,k,0) = ulz(i,j,k,0);
                    urz(i,j,k,1) = ulz(i,j,k,1);
                    ulz(i,j,k,2) = 0.0;
                    urz(i,j,k,2) = 0.0;
                case NoSlipWall:
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        ulz(i,j,k,n) = 0.0;
                        urz(i,j,k,n) = 0.0;
                    }
                case Outflow:
                    ulz(i,j,k,2) = max(ulz(i,j,k,2),0.0);
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        urz(i,j,k,n) = ulz(i,j,k,n);
                    }
                case Interior:
                        break;
                default:
                    Print() << "velpred_3d: invalid boundary type phys_bc(3,2)" << std::endl;
            }
        }

        // No need to compute uimhz(:,:,:,2) since it's equal to wtrans-w0
        // upwind using full velocity to get transverse components of uimhz
        // Note: wtrans already contains w0
        uimhz(i,j,k,0) = wtrans(i,j,k) > 0.0 ? ulz(i,j,k,0) : urz(i,j,k,0);
        uimhz(i,j,k,0) = fabs(wtrans(i,j,k)) < rel_eps ? 
            0.5*(ulz(i,j,k,0)+urz(i,j,k,0)) : uimhz(i,j,k,0);

        uimhz(i,j,k,1) = wtrans(i,j,k) > 0.0 ? ulz(i,j,k,1) : urz(i,j,k,1);
        uimhz(i,j,k,1) = fabs(wtrans(i,j,k)) < rel_eps ? 
            0.5*(ulz(i,j,k,1)+urz(i,j,k,1)) : uimhz(i,j,k,1);
    });

}

#endif
#include <Maestro.H>
#include <MaestroHydro_F.H>

using namespace amrex;

#if (AMREX_SPACEDIM == 2)

void 
Maestro::VelPredInterface(const MFIter& mfi, 
                        Array4<Real> const utilde,
                        Array4<Real> const ufull,
                        Array4<Real> const utrans,
                        Array4<Real> const Imu,
                        Array4<Real> const Ipu,
                        Array4<Real> const Imv,
                        Array4<Real> const Ipv,
                        Array4<Real> const ul,
                        Array4<Real> const ur,
                        Array4<Real> const uimh,
                        const Box& domainBox,
                        const Real* dx)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::VelPredInterface()",VelPredInterface);

    const Box& tileBox = mfi.tilebox();
    const Box& obx = amrex::grow(tileBox, 1);
    const Box& mxbx = amrex::growLo(obx, 0, -1);
    const Box& mybx = amrex::growLo(obx, 1, -1);

    Real rel_eps;
    get_rel_eps(&rel_eps);

    IntVector bclo(ncomp);
    IntVector bchi(ncomp);

    int * AMREX_RESTRICT bclo_p = bclo.dataPtr();
    int * AMREX_RESTRICT bchi_p = bchi.dataPtr();

    Real dt2 = 0.5 * dt;
    Real dt4 = 0.25 * dt;

    Real hx = dx[0];
    Real hy = dx[1];

    // x-direction
    int ilo = domainBox.loVect()[0];
    int ihi = domainBox.hiVect()[0];

    int bclo = physbc.lo()[0];
    int bchi = physbc.hi()[0];

    int ppm_type_local = ppm_type;

    // NOTE: for ppm_type == 0, slopex == Ipu, slopey == Imv

    ////////////////////////////////////
    // Create u_{\i-\half\e_x}^x, etc.
    ////////////////////////////////////

    AMREX_PARALLEL_FOR_3D(mxbx, i, j, k, 
    {
        if (ppm_type_local == 0) {
            Real maxu = max(0.0,ufull(i-1,j,k,0));
            Real minu = min(0.0,ufull(i  ,j,k,0));
            // extrapolate both components of velocity to left face
            ul(i,j,k,0) = utilde(i-1,j,k,0) + (0.5 - (dt2/hx)*maxu)*Ipu(i-1,j,k,0);
            ul(i,j,k,1) = utilde(i-1,j,k,1) + (0.5 - (dt2/hx)*maxu)*Ipu(i-1,j,k,1);
            // extrapolate both components of velocity to right face
            ur(i,j,k,0) = utilde(i  ,j,k,0) - (0.5 + (dt2/hx)*minu)*Ipu(i  ,j,k,0);
            ur(i,j,k,1) = utilde(i  ,j,k,1) - (0.5 + (dt2/hx)*minu)*Ipu(i  ,j,k,1);
        } else if (ppm_type_local == 1 || ppm_type_local == 2) {
            // extrapolate both components of velocity to left face
            ul(i,j,k,0) = Ipu(i-1,j,k,0);
            ul(i,j,k,1) = Ipv(i-1,j,k,0);
            // extrapolate both components of velocity to right face
            ur(i,j,k,0) = Imu(i,j,k,0);
            ur(i,j,k,1) = Imv(i,j,k,0);
        }

        // impose lo side bc's
        if (i == ilo) {
            switch(bclo) {
                case Inflow:
                    ul(i,j,k,0:1) = utilde(i-1,j,k,0:1);
                    ur(i,j,k,0:1) = utilde(i-1,j,k,0:1);
                case SlipWall:
                case Symmetry: 
                    ul(i,j,k,0) = 0.0;
                    ur(i,j,k,0) = 0.0;
                    ul(i,j,k,1) = ur(i,j,k,1);
                case NoSlipWall:
                    ul(i,j,k,0:1) = 0.0;
                    ur(i,j,k,0:1) = 0.0;
                case Outflow:
                    ur(i,j,k,0) = min(ur(i,j,k,0),0.0);
                    ur(i,j,k,0:1) = ul(i,j,k,0:1);
                case Interior:
                default:
                    Print() << "velpred_2d: invalid boundary type phys_bc(1,1)" << std::endl;
            }

        // impose hi side bc's
        } else if (i == ihi+1) {
            switch (bchi) {
                case Inflow:
                    ul(i,j,k,0:1) = utilde(i,j,k,0:1);
                    ur(i,j,k,0:1) = utilde(i,j,k,0:1);
                case SlipWall:
                case Symmetry:
                    ul(i,j,k,0) = 0.0;
                    ur(i,j,k,0) = 0.0;
                    ur(i,j,k,1) = ul(i,j,k,1);
                case NoSlipWall:
                    ul(i,j,k,0:1) = 0.0;
                    ur(i,j,k,0:1) = 0.0;
                case Outflow:
                    ul(i,j,k,0) = max(ul(i,j,k,0),0.0);
                    ur(i,j,k,0:1) = ul(i,j,k,0:1);
                case Interior:
                default:
                    Print() << "velpred_2d: invalid boundary type phys_bc(1,2)" << std::endl;
            }
        }

        // No need to compute uimh(:,:,0) since it's equal to utrans-w0
        // upwind using full velocity to get transverse component of uimhx
        // Note: utrans already contains w0
        uimh(i,j,k,1) = utrans(i,j,k) > 0.0 ? ul(i,j,k,1) : ur(i,j,k,1);
        uimh(i,j,k,1) = fabs(utrans(i,j,k)) < rel_eps ? 0.5*(ul(i,j,k,1)+ur(i,j,k,1)) : uimh(i,j,k,1);

    });

    // y-direction
    int jlo = domainBox.loVect()[1];
    int jhi = domainBox.hiVect()[1];
    bclo = physbc.lo()[1];
    bchi = physbc.hi()[1];

    AMREX_PARALLEL_FOR_3D(mybx, i, j, k, 
    {
        if (ppm_type_local == 0) {
            Real maxu = max(0.0,ufull(i,j-1,k,1));
            Real minu = min(0.0,ufull(i,j,k,1));
            // extrapolate both components of velocity to left face
            ul(i,j,k,0) = utilde(i,j-1,k,0) + (0.5-(dt2/hy)*maxu)*Imv(i,j-1,k,0);
            ul(i,j,k,1) = utilde(i,j-1,k,1) + (0.5-(dt2/hy)*maxu)*Imv(i,j-1,k,1);
            // extrapolate both components of velocity to right face
            ur(i,j,k,0) = utilde(i,j,k,0) - (0.5+(dt2/hy)*minu)*Imv(i,j,k,0);
            ur(i,j,k,1) = utilde(i,j,k,1) - (0.5+(dt2/hy)*minu)*Imv(i,j,k,1);
        } else if (ppm_type_local == 1 || ppm_type_local == 2) {
            // extrapolate both components of velocity to left face
            ul(i,j,k,0) = Ipu(i,j-1,k,0);
            ul(i,j,k,1) = Ipv(i,j-1,k,1);
            // extrapolate both components of velocity to right face
            ur(i,j,k,0) = Imu(i,j,k,0);
            ur(i,j,k,1) = Imv(i,j,k,1);
        }

        // impose lo side bc's
        if (j == jlo) {
            switch (bclo) {
                case (Inflow)
                    ul(i,j,k,1:2) = utilde(i,j-1,k,1:2)
                    ur(i,j,k,1:2) = utilde(i,j-1,k,1:2)
                case (SlipWall, Symmetry)
                    ul(i,j,k,1) = ur(i,j,k,1)
                    ul(i,j,k,2) = 0.0
                    ur(i,j,k,2) = 0.0
                case (NoSlipWall)
                    ul(i,j,k,1:2) = 0.0
                    ur(i,j,k,1:2) = 0.0
                case (Outflow)
                    ur(i,j,k,2) = min(ur(i,j,k,2),0.0)
                    ul(i,j,k,1:2) = ur(i,j,k,1:2)
                case (Interior)
                default:
                    Print() << "velpred_2d: invalid boundary type phys_bc(2,1)" << std::endl;
            }

        // impose hi side bc's
        } else if (j == jhi+1) {
            switch (bchi) {
                case (Inflow)
                    ul(i,j,k,1:2) = utilde(i,j,k,1:2)
                    ur(i,j,k,1:2) = utilde(i,j,k,1:2)
                case (SlipWall, Symmetry)
                    ur(i,j,k,1) = ul(i,j,k,1)
                    ul(i,j,k,2) = 0.0
                    ur(i,j,k,2) = 0.0
                case (NoSlipWall)
                    ul(i,j,k,1:2) = 0.0
                    ur(i,j,k,1:2) = 0.0
                case (Outflow)
                    ul(i,j,k,2)   = max(ul(i,j,k,2),0.0)
                    ur(i,j,k,1:2) = ul(i,j,k,1:2)
                case (Interior)
                case  default
    #ifndef AMREX_USE_CUDA
                    call amrex_error("velpred_2d: invalid boundary type phys_bc(2,2)")
    #endif
            }
        }
        // No need to compute uimh(:,:,1) since it's equal to utrans-w0
        // upwind using full velocity to get transverse component of uimhy
        // Note: utrans already contains w0
        uimh(i,j,k,0) = merge(ul(i,j,k,0),ur(i,j,k,0),utrans(i,j,k)>0.0)
        uavg = 0.5*(ul(i,j,k,0)+ur(i,j,k,0))
        uimh(i,j,k,0) = merge(uavg,uimh(i,j,k,0),fabs(utrans(i,j,k))<rel_eps)

    });
}

#else 

void 
Maestro::VelPredInterface(const MFIter& mfi, 
                Array4<Real> const s,
                Array4<Real> const slx,
                const Box& domainBox,
                const Vector<BCRec>& bcs,
                int ncomp,
                int bc_start_comp)
{

}

#endif
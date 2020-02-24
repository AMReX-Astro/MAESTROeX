
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void
Maestro::FillUmacGhost(const Box& domainBox, const Box& bx, 
                       const Array4<Real> umac,
                       const Array4<Real> vmac
#if (AMREX_SPACEDIM == 3)
                       , const Array4<Real> wmac
#endif
                       )
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::FillUmacGhost()", FillUmacGhost);

    const auto domlo = domainBox.loVect3d();
    const auto domhi = domainBox.hiVect3d();

    const auto xbx = mfi.grownnodaltilebox(0, 1);
    const auto ybx = mfi.grownnodaltilebox(1, 1);
    const auto zbx = mfi.grownnodaltilebox(2, 1);

    AMREX_PARALLEL_FOR_3D(xbx, i, j, k, 
    {
        // lo x-faces
        if (i == domlo[0] - 1) {
            switch(phys_bc[0]) {
                case Inflow:
                    umac(i,j,k) = umac(i+1,j,k);
                    vmac(i,j,k) = 0.0;
#if (AMREX_SPACEDIM == 3)
                    wmac(i,j,k) = 0.0;
#endif
                    break;
                case SlipWall:
                case NoSlipWall:
                    umac(i,j,k) = 0.0;
                    vmac(i,j,k) = 0.0;
#if (AMREX_SPACEDIM == 3)
                    wmac(i,j,k) = 0.0;
#endif
                    break;
                case Outflow:
                    umac(i,j,k) = umac(i+1,j,k);
                    vmac(i,j,k) = vmac(i+1,j,k);
#if (AMREX_SPACEDIM == 3)
                    wmac(i,j,k) = wmac(i+1,j,k);
#endif
                    break;
                case Symmetry:
                    umac(i,j,k) = -umac(i+2,j,k);
                    vmac(i,j,k) = vmac(i+1,j,k);
#if (AMREX_SPACEDIM == 3)
                    wmac(i,j,k) = wmac(i+1,j,k);
#endif
                    break;
                case Interior:
                // do nothing
                    break;
            }
        }

        // hi x-faces
        if (i == domhi[0]+2) {
            switch(phys_bc[AMREX_SPACEDIM]) {
                case Inflow:
                    umac(i,j,k) = umac(i-1,j,k);
                    vmac(i-1,j,k) = 0.0;
#if (AMREX_SPACEDIM == 3)
                    wmac(i-1,j,k) = 0.0;
#endif
                    break;
                case SlipWall:
                case NoSlipWall:
                    umac(i,j,k) = 0.0;
                    vmac(i-1,j,k) = 0.0;
#if (AMREX_SPACEDIM == 3)
                    wmac(i-1,j,k) = 0.0;
#endif
                    break;
                case Outflow:
                    umac(i,j,k) = umac(i-1,j,k);
                    vmac(i-1,j,k) = vmac(i-2,j,k);
#if (AMREX_SPACEDIM == 3)
                    wmac(i-1,j,:) = wmac(i-2,j,k);
#endif
                    break;
                case Symmetry:
                    umac(i,j,k) = -umac(i-2,j,k);
                    vmac(i-1,j,k) = vmac(i-2,j,k);
#if (AMREX_SPACEDIM == 3)
                    wmac(i-1,j,k) = wmac(i-2,j,k);
#endif
                    break;
                case Interior:
                    // do nothing
                    break;        
            }
        }
    });

    AMREX_PARALLEL_FOR_3D(ybx, i, j, k, 
    {
        // lo y-faces
        if (j == domlo[1]-1) {
            switch (phys_bc[1]) {
                case Inflow:
                    umac(i,j,k) = 0.0;
                    vmac(i,j,k) = vmac(i,j+1,k)
#if (AMREX_SPACEDIM == 3)
                    wmac(i,j,k) = 0.0;
#endif
                    break;
                case SlipWall:
                case NoSlipWall:
                    umac(i,j,k) = 0.0;
                    vmac(i,j,k) = 0.0;
#if (AMREX_SPACEDIM == 3)
                    wmac(i,j,k) = 0.0;
#endif
                    break;
                case Outflow:
                    umac(i,j,k) = umac(i,j+1,k);
                    vmac(i,j,k) = vmac(i,j+1,k);
#if (AMREX_SPACEDIM == 3)
                    wmac(i,j,k) = wmac(i,j+1,k);
#endif
                    break;
                case Symmetry:
                    umac(i,j,k) = umac(i,j+1,k);
                    vmac(i,j,k) = -vmac(i,j+2,k);
#if (AMREX_SPACEDIM == 3)
                    wmac(i,j,k) = wmac(i,j+1,k);
#endif
                    break;
                case Interior:
                    // do nothing
                    break;            
            }
        }

        // hi y-faces
        if (bx.hi(1) == domhi[1]+2) {
            switch (phys_bc[AMREX_SPACEDIM+1]) {
                case Inflow:
                    umac(i,j-1,k) = 0.0;
                    vmac(i,j,k) = vmac(i,j-1,k);
#if (AMREX_SPACEDIM == 3)
                    wmac(i,j-1,k) = 0.0;
#endif
                    break;
                case SlipWall:
                case NoSlipWall:
                    umac(i,j-1,k) = 0.0;
                    vmac(i,j,k) = 0.0;
#if (AMREX_SPACEDIM == 3)
                    wmac(i,j-1,k) = 0.0;
#endif
                    break;
                case Outflow:
                    umac(i,j-1,k) = umac(i,j-2,k);
                    vmac(i,j,k) = vmac(i,j-1,k);
#if (AMREX_SPACEDIM == 3)
                    wmac(i,j-1,k) = wmac(i,j-2,k);
#endif
                    break;
                case Symmetry:
                    umac(i,j-1,k) = umac(i,j-2,k);
                    vmac(i,j,k) = -vmac(i,j-2,k);
#if (AMREX_SPACEDIM == 3)
                    wmac(i,j-1,k) = wmac(i,j-2,k);
#endif
                    break;
                case Interior:
                    // do nothing
                    break;            
            }
        }
    });

#if (AMREX_SPACEDIM == 3)

    AMREX_PARALLEL_FOR_3D(zbx, i, j, k, 
    {
        // lo z-faces
        if (bx.lo(2) == domlo[2]-1) {
            switch (phys_bc[2]) {
                case Inflow:
                    umac(i,j,k) = 0.0;
                    vmac(i,j,k) = 0.0;
                    wmac(i,j,k) = wmac(i,j,k+1);
                    break;
                case SlipWall:
                case NoSlipWall:
                    umac(i,j,k) = 0.0;
                    vmac(i,j,k) = 0.0;
                    wmac(i,j,k) = 0.0;
                    break;
                case Outflow:
                    umac(i,j,k) = umac(i,j,k+1);
                    vmac(i,j,k) = vmac(i,j,k+1);
                    wmac(i,j,k) = wmac(i,j,k+1);
                    break;
                case Symmetry:
                    umac(i,j,k) = umac(i,j,k+1);
                    vmac(i,j,k) = vmac(i,j,k+1);
                    wmac(i,j,k) = -wmac(i,j,k+2);
                    break;
                case Interior:
                    // do nothing
                    break;        
            }
        }

        // hi z-faces
        if (bx.hi(2) == domhi[2]+2) {
            switch (phys_bc[2+AMREX_SPACEDIM]) {
                case Inflow:
                    umac(i,j,k-1) = 0.0;
                    vmac(i,j,k-1) = 0.0;
                    wmac(i,j,k) = wmac(i,j,k-1);
                    break;
                case SlipWall:
                case NoSlipWall:
                    umac(i,j,k-1) = 0.0;
                    vmac(i,j,k-1) = 0.0;
                    wmac(i,j,k) = 0.0;
                    break;
                case Outflow:
                    umac(i,j,k-1) = umac(i,j,k-2);
                    vmac(i,j,k-1) = vmac(i,j,k-2);
                    wmac(i,j,k) = wmac(i,j,k-1);
                    break;
                case Symmetry:
                    umac(i,j,k-1) = umac(i,j,k-2);
                    vmac(i,j,k-1) = vmac(i,j,k-2);
                    wmac(i,j,k) = -wmac(i,j,k-2);
                    break;
                case Interior:
                    // do nothing
                    break;        
            }
        }
    });
#endif
}

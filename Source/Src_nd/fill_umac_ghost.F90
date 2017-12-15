
#include "AMReX_BC_TYPES.H"

module fill_umac_ghost_module
  
  implicit none

  private

contains

  subroutine fill_umac_ghost(domlo, domhi, lo, hi, &
                             umac, umac_lo, umac_hi, &
                             vmac, vmac_lo, vmac_hi, &
#if (AMREX_SPACEDIM == 3)
                             wmac, wmac_lo, wmac_hi, &
#endif
                             lo_bc, hi_bc) bind(C, name="fill_umac_ghost")

    
    integer         , intent(in   ) :: domlo(3), domhi(3), lo(3), hi(3)
    integer         , intent(in   ) :: umac_lo(3), umac_hi(3)
    integer         , intent(in   ) :: vmac_lo(3), vmac_hi(3)
    double precision, intent(inout) :: umac(umac_lo(1):umac_hi(1), &
                                            umac_lo(2):umac_hi(2), &
                                            umac_lo(3):umac_hi(3))
    double precision, intent(inout) :: vmac(vmac_lo(1):vmac_hi(1), &
                                            vmac_lo(2):vmac_hi(2), &
                                            vmac_lo(3):vmac_hi(3))
#if (AMREX_SPACEDIM == 3)
    integer         , intent(in   ) :: wmac_lo(3), wmac_hi(3)
    double precision, intent(inout) :: wmac(wmac_lo(1):wmac_hi(1), &
                                            wmac_lo(2):wmac_hi(2), &
                                            wmac_lo(3):wmac_hi(3))
#endif
    integer         , intent(in   ) :: lo_bc(AMREX_SPACEDIM), hi_bc(AMREX_SPACEDIM)

    ! lo x-faces
    if (lo(1) .eq. domlo(1)) then
       select case (lo_bc(1))
       case (Inflow)
          umac(lo(1)-1,:,:) = umac(lo(1),:,:)
          vmac(lo(1)-1,:,:) = 0.d0
#if (AMREX_SPACEDIM == 3)
          wmac(lo(1)-1,:,:) = 0.d0
#endif
       case (SlipWall, NoSlipWall)
          umac(lo(1)-1,:,:) = 0.d0
          vmac(lo(1)-1,:,:) = 0.d0
#if (AMREX_SPACEDIM == 3)
          wmac(lo(1)-1,:,:) = 0.d0
#endif
       case (Outflow)
          umac(lo(1)-1,:,:) = umac(lo(1),:,:)
          vmac(lo(1)-1,:,:) = vmac(lo(1),:,:)
#if (AMREX_SPACEDIM == 3)
          wmac(lo(1)-1,:,:) = wmac(lo(1),:,:)
#endif
       case (Symmetry)
          umac(lo(1)-1,:,:) = -umac(lo(1)+1,:,:)
          vmac(lo(1)-1,:,:) = vmac(lo(1),:,:)
#if (AMREX_SPACEDIM == 3)
          wmac(lo(1)-1,:,:) = wmac(lo(1),:,:)
#endif
       case (Interior)
          ! do nothing
       case default
       end select
    end if

    ! hi x-faces
    if (hi(1) .eq. domhi(1)) then
       select case (hi_bc(1))
       case (Inflow)
          umac(hi(1)+2,:,:) = umac(hi(1)+1,:,:)
          vmac(hi(1)+1,:,:) = 0.d0
#if (AMREX_SPACEDIM == 3)
          wmac(hi(1)+1,:,:) = 0.d0
#endif
       case (SlipWall, NoSlipWall)
          umac(hi(1)+2,:,:) = 0.d0
          vmac(hi(1)+1,:,:) = 0.d0
#if (AMREX_SPACEDIM == 3)
          wmac(hi(1)+1,:,:) = 0.d0
#endif
       case (Outflow)
          umac(hi(1)+2,:,:) = umac(hi(1)+1,:,:)
          vmac(hi(1)+1,:,:) = vmac(hi(1),:,:)
#if (AMREX_SPACEDIM == 3)
          wmac(hi(1)+1,:,:) = wmac(hi(1),:,:)
#endif
       case (Symmetry)
          umac(hi(1)+2,:,:) = -umac(hi(1),:,:)
          vmac(hi(1)+1,:,:) = vmac(hi(1),:,:)
#if (AMREX_SPACEDIM == 3)
          wmac(hi(1)+1,:,:) = wmac(hi(1),:,:)
#endif
       case (Interior)
          ! do nothing
       case default
       end select
    end if

    ! lo y-faces
    if (lo(2) .eq. domlo(2)) then
       select case (lo_bc(2))
       case (Inflow)
          umac(:,lo(2)-1,:) = 0.d0
          vmac(:,lo(2)-1,:) = vmac(:,lo(2),:)
#if (AMREX_SPACEDIM == 3)
          wmac(:,lo(2)-1,:) = 0.d0
#endif
       case (SlipWall, NoSlipWall)
          umac(:,lo(2)-1,:) = 0.d0
          vmac(:,lo(2)-1,:) = 0.d0
#if (AMREX_SPACEDIM == 3)
          wmac(:,lo(2)-1,:) = 0.d0
#endif
       case (Outflow)
          umac(:,lo(2)-1,:) = umac(:,lo(2),:)
          vmac(:,lo(2)-1,:) = vmac(:,lo(2),:)
#if (AMREX_SPACEDIM == 3)
          wmac(:,lo(2)-1,:) = wmac(:,lo(2),:)
#endif
       case (Symmetry)
          umac(:,lo(2)-1,:) = umac(:,lo(2),:)
          vmac(:,lo(2)-1,:) = -vmac(:,lo(2)+1,:)
#if (AMREX_SPACEDIM == 3)
          wmac(:,lo(2)-1,:) = wmac(:,lo(2),:)
#endif
       case (Interior)
          ! do nothing
       case default
       end select
    end if

    ! hi y-faces
    if (hi(2) .eq. domhi(2)) then
       select case (hi_bc(2))
       case (Inflow)
          umac(:,hi(2)+1,:) = 0.d0
          vmac(:,hi(2)+2,:) = vmac(:,hi(2)+1,:)
#if (AMREX_SPACEDIM == 3)
          wmac(:,hi(2)+1,:) = 0.d0
#endif
       case (SlipWall, NoSlipWall)
          umac(:,hi(2)+1,:) = 0.d0
          vmac(:,hi(2)+2,:) = 0.d0
#if (AMREX_SPACEDIM == 3)
          wmac(:,hi(2)+1,:) = 0.d0
#endif
       case (Outflow)
          umac(:,hi(2)+1,:) = umac(:,hi(2),:)
          vmac(:,hi(2)+2,:) = vmac(:,hi(2)+1,:)
#if (AMREX_SPACEDIM == 3)
          wmac(:,hi(2)+1,:) = wmac(:,hi(2),:)
#endif
       case (Symmetry)
          umac(:,hi(2)+1,:) = umac(:,hi(2),:)
          vmac(:,hi(2)+2,:) = -vmac(:,hi(2),:)
#if (AMREX_SPACEDIM == 3)
          wmac(:,hi(2)+1,:) = wmac(:,hi(2),:)
#endif
       case (Interior)
          ! do nothing
       case default
       end select
    end if

#if (AMREX_SPACEDIM == 3)

    ! lo z-faces
    if (lo(3) .eq. domlo(3)) then
       select case (lo_bc(3))
       case (Inflow)
          umac(:,:,lo(3)-1) = 0.d0
          vmac(:,:,lo(3)-1) = 0.d0
          wmac(:,:,lo(3)-1) = wmac(:,:,lo(3))
       case (SlipWall, NoSlipWall)
          umac(:,:,lo(3)-1) = 0.d0
          vmac(:,:,lo(3)-1) = 0.d0
          wmac(:,:,lo(3)-1) = 0.d0
       case (Outflow)
          umac(:,:,lo(3)-1) = umac(:,:,lo(3))
          vmac(:,:,lo(3)-1) = vmac(:,:,lo(3))
          wmac(:,:,lo(3)-1) = wmac(:,:,lo(3))
       case (Symmetry)
          umac(:,:,lo(3)-1) = umac(:,:,lo(3))
          vmac(:,:,lo(3)-1) = vmac(:,:,lo(3))
          wmac(:,:,lo(3)-1) = -wmac(:,:,lo(3)+1)
       case (Interior)
          ! do nothing
       case default
       end select
    end if

    ! hi z-faces
    if (hi(3) .eq. domhi(3)) then
       select case (hi_bc(3))
       case (Inflow)
          umac(:,:,hi(3)+1) = 0.d0
          vmac(:,:,hi(3)+1) = 0.d0
          wmac(:,:,hi(3)+2) = wmac(:,:,hi(3)+1)
       case (SlipWall, NoSlipWall)
          umac(:,:,hi(3)+1) = 0.d0
          vmac(:,:,hi(3)+1) = 0.d0
          wmac(:,:,hi(3)+2) = 0.d0
       case (Outflow)
          umac(:,:,hi(3)+1) = umac(:,:,hi(3))
          vmac(:,:,hi(3)+1) = vmac(:,:,hi(3))
          wmac(:,:,hi(3)+2) = wmac(:,:,hi(3)+1)
       case (Symmetry)
          umac(:,:,hi(3)+1) = umac(:,:,hi(3))
          vmac(:,:,hi(3)+1) = vmac(:,:,hi(3))
          wmac(:,:,hi(3)+2) = -wmac(:,:,hi(3))
       case (Interior)
          ! do nothing
       case default
       end select
    end if

#endif
    
  end subroutine fill_umac_ghost

end module fill_umac_ghost_module

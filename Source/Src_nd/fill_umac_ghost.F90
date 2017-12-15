
module fill_umac_ghost_module
  
  implicit none

  private

contains

  subroutine fill_umac_ghost(lo, hi, &
                             umac, umac_lo, umac_hi, &
                             vmac, vmac_lo, vmac_hi, &
#if (AMREX_SPACEDIM == 3)
                             wmac, wmac_lo, wmac_hi, &
#endif
                             lo_bc, hi_bc) bind(C, name="fill_umac_ghost")

    
    integer         , intent(in   ) :: lo(3), hi(3)
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

    
  end subroutine fill_umac_ghost

end module fill_umac_ghost_module

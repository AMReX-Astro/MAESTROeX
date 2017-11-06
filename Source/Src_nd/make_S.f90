subroutine make_S_cc(lo, hi,S_cc, s_lo, s_hi, &
                     ng, ncomp) bind (C,name="make_S_cc")

  use eos_type_module
  use eos_module
  use network, only: nspec

  integer         , intent (in   ) :: lo(3), hi(3), s_lo(3), s_hi(3), ng, ncomp
  double precision, intent (inout) :: S_cc(s_lo(1):s_hi(1), &
                                           s_lo(2):s_hi(2), &
                                           s_lo(3):s_hi(3),1)

  integer i,j,k
  integer :: pt_index(3)
  type(eos_t) :: eos_state

  ! loop over the data
  do k = lo(3),hi(3)
  do j = lo(2),hi(2)
  do i = lo(1),hi(1)

     eos_state%rho   = 1.d0
     eos_state%T     = 300.d0
     eos_state%xn(1:nspec-1) = 0.01d0
     eos_state%xn(nspec) = 1- sum(eos_state%xn(1:nspec-1))

     ! dens, temp, and xmass are inputs
     call eos(eos_input_rt, eos_state)

     S_cc(i,j,k,1) = 0.d0

  enddo
  enddo
  enddo

end subroutine make_S_cc

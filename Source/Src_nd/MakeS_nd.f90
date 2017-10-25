subroutine make_S_cc_work(lo, hi, S_cc, &
                          s_lo, s_hi) bind (C,name="make_S_cc_work")

  integer         , intent (in   ) :: lo(3), hi(3)
  integer         , intent (in   ) :: s_lo(3), s_hi(3)
  double precision, intent (inout) :: S_cc(s_lo(1):s_hi(1), &
                                           s_lo(2):s_hi(2), &
                                           s_lo(3):s_hi(3),1)

  integer i,j,k

  ! loop over the data
  do k = lo(3),hi(3)
  do j = lo(2),hi(2)
  do i = lo(1),hi(1)
     S_cc(i,j,k,1) = 0.d0
  enddo
  enddo
  enddo

end subroutine make_S_cc_work

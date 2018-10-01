module tagging_module

  use parallel, only: parallel_IOProcessor
  use meth_params_module, only: temp_comp, nscal

  implicit none

  private

contains

  ! ::: -----------------------------------------------------------
  ! ::: This routine will tag high error cells based on the state
  ! :::
  ! ::: INPUTS/OUTPUTS:
  ! :::
  ! ::: tag        <=  integer tag array
  ! ::: tag_lo,hi   => index extent of tag array
  ! ::: state       => state array
  ! ::: state_lo,hi => index extent of state array
  ! ::: set         => integer value to tag cell for refinement
  ! ::: clear       => integer value to untag cell
  ! ::: lo,hi       => work region we are allowed to change
  ! ::: dx          => cell size
  ! ::: time        => problem evolution time
  ! ::: level       => refinement level of this array
  ! ::: -----------------------------------------------------------

  subroutine state_error(tag,tag_lo,tag_hi, &
       state,state_lo,state_hi, &
       set,clear,&
       lo,hi,&
       dx,time,tag_err) bind(C, name="state_error")

    integer          :: lo(3),hi(3)
    integer          :: state_lo(3),state_hi(3)
    integer          :: tag_lo(3),tag_hi(3)
    double precision :: state(state_lo(1):state_hi(1), &
         state_lo(2):state_hi(2), &
         state_lo(3):state_hi(3), 1:nscal)
    integer          :: tag(tag_lo(1):tag_hi(1),tag_lo(2):tag_hi(2),tag_lo(3):tag_hi(3))
    double precision :: dx(3),time
    double precision :: tag_err(2)
    integer          :: set,clear

    ! local
    integer          :: i, j, k
    double precision :: temperr, denserr


    if (parallel_IOProcessor()) then
       print*,"Create a local copy of tagging.f90 in your build directory"
       print*,"Here is a sample that tags the temperature above temperr"
    end if

    ! abort program
    call bl_error()


    ! set temperature and density flags
    temperr = tag_err(1)
    denserr = tag_err(2)

    ! Tag on regions of high temperature
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (state(i,j,k,temp_comp) .ge. temperr) then
                tag(i,j,k) = set
             endif
          enddo
       enddo
    enddo

  end subroutine state_error

end module tagging_module

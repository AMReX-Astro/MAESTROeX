module tagging_module

  use amrex_error_module
  use meth_params_module, only: temp_comp, rho_comp, nscal
  use probin_module, only: tag_density_1, tag_density_2, tag_density_3
  use base_state_geometry_module, only: nr_fine, max_radial_level

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
       dx,time,&
       lev,tag_array) bind(C, name="state_error")

    integer          :: lo(3),hi(3)
    integer          :: state_lo(3),state_hi(3)
    integer          :: tag_lo(3),tag_hi(3)
    double precision :: state(state_lo(1):state_hi(1), &
         state_lo(2):state_hi(2), &
         state_lo(3):state_hi(3), 1:nscal)
    integer          :: tag(tag_lo(1):tag_hi(1),tag_lo(2):tag_hi(2),tag_lo(3):tag_hi(3))
    double precision :: dx(3),time
    integer          :: tag_array(0:max_radial_level,0:nr_fine-1)
    integer          :: set,clear,lev

    ! local
    integer          :: i, j, k, r

    ! Tag on regions of high density
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (((lev .eq. 0) .and. (state(i,j,k,rho_comp) .ge. tag_density_1)) .or. &
	     	 ((lev .eq. 1) .and. (state(i,j,k,rho_comp) .ge. tag_density_2)) .or. &
		 ((lev .eq. 2) .and. (state(i,j,k,rho_comp) .ge. tag_density_3))) then
                tag(i,j,k) = set
	     end if	     
          enddo
       enddo
    enddo

  end subroutine state_error

  subroutine tag_boxes(tag,tag_lo,tag_hi, &
                        set,clear,&
                        lo,hi,&
                        dx,time,&
                        lev,tag_array) bind(C, name="tag_boxes")

    integer          :: lo(3),hi(3)
    integer          :: tag_lo(3),tag_hi(3)
    integer          :: tag(tag_lo(1):tag_hi(1),tag_lo(2):tag_hi(2),tag_lo(3):tag_hi(3))
    double precision :: dx(3),time
    integer          :: tag_array(0:max_radial_level,0:nr_fine-1)
    integer          :: set,clear,lev

    call amrex_error("tag_boxes not needed for spherical")
       
  end subroutine tag_boxes
  
  subroutine retag_array(set,clear,&
                          lo,hi,&
                          lev,tag_array) bind(C, name="retag_array")

    integer          :: lo(3),hi(3)
    integer          :: tag_array(0:max_radial_level,0:nr_fine-1)
    integer          :: set,clear,lev

    call amrex_error("retag_array not needed for spherical")
       
  end subroutine retag_array
  
end module tagging_module

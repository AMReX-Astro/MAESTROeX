
module update_scal_module

  use amrex_constants_module

  implicit none

  private
  
contains
  
#if (AMREX_SPACEDIM == 1)
  subroutine update_rhoX_1d(lo, hi, &
                              sold,   so_lo, so_hi, nc_so, &
                              snew,   sn_lo, sn_hi, nc_sn, &
                              sfluxx, x_lo, x_hi, nc_x, &
                              force,  f_lo, f_hi, nc_f, & 
                              dx, dt, &
                              startcomp, endcomp) bind(C,name="update_rhoX_1d")

    integer         , intent(in   ) :: lo(1), hi(1)
    integer         , intent(in   ) :: so_lo(1), so_hi(1), nc_so
    double precision, intent(in   ) :: sold  (so_lo(1):so_hi(1),nc_so)
    integer         , intent(in   ) :: sn_lo(1), sn_hi(1), nc_sn
    double precision, intent(inout) :: snew  (sn_lo(1):sn_hi(1),nc_sn)
    integer         , intent(in   ) :: x_lo(1), x_hi(1), nc_x
    double precision, intent(in   ) :: sfluxx(x_lo(1):x_hi(1),nc_x)
    integer         , intent(in   ) :: f_lo(1), f_hi(1), nc_f
    double precision, intent(in   ) :: force (f_lo(1):f_hi(1),nc_f)
    double precision, intent(in   ) :: dx(1), dt
    integer         , intent(in   ) :: startcomp, endcomp

    ! Local variables
    integer          :: i, comp, comp2
    double precision :: divterm 
    double precision :: delta, frac, sumX
    logical          :: has_negative_species
    ! from MAESTRO parameters file
    double precision :: base_cutoff_density = 3.d6

    do comp = startcomp, endcomp
       do i=lo(1),hi(1)
          divterm = (sfluxx(i+1,comp) - sfluxx(i,comp))/dx(1)
          snew(i,comp) = sold(i,comp) + dt*(-divterm + force(i,comp))
       end do
    enddo

    ! update density
    snew(:,1) = sold(:,1)
       
    do i = lo(1), hi(1)

       has_negative_species = .false.

       ! define the update to rho as the sum of the updates to (rho X)_i
       do comp = startcomp, endcomp
          snew(i,1) = snew(i,1) + (snew(i,comp)-sold(i,comp))
          if (snew(i,comp) .lt. ZERO) has_negative_species = .true.
       enddo

       ! enforce a density floor
       if (snew(i,1) .lt. HALF*base_cutoff_density) then
          do comp = startcomp, endcomp
             snew(i,comp) = snew(i,comp) * HALF*base_cutoff_density/snew(i,1)
          end do
          snew(i,1) = HALF*base_cutoff_density
       end if

       ! do not allow the species to leave here negative.
       if (has_negative_species) then
          do comp = startcomp, endcomp
             if (snew(i,comp) .lt. ZERO) then
                delta = -snew(i,comp)
                sumX = ZERO 
                do comp2 = startcomp, endcomp
                   if (comp2 .ne. comp .and. snew(i,comp2) .ge. ZERO) then
                      sumX = sumX + snew(i,comp2)
                   end if
                enddo
                do comp2 = startcomp, endcomp
                   if (comp2 .ne. comp .and. snew(i,comp2) .ge. ZERO) then
                      frac = snew(i,comp2) / sumX
                      snew(i,comp2) = snew(i,comp2) - frac * delta
                   end if
                enddo
                snew(i,comp) = ZERO
             end if
          end do
       end if
    end do

  end subroutine update_rhoX_1d
#endif

#if (AMREX_SPACEDIM == 2)
  subroutine update_rhoX_2d(lo, hi, &
                              sold,   so_lo, so_hi, nc_so, &
                              snew,   sn_lo, sn_hi, nc_sn, &
                              sfluxx, x_lo, x_hi, nc_x, &
                              sfluxy, y_lo, y_hi, nc_y, &
                              force,  f_lo, f_hi, nc_f, & 
                              dx, dt, &
                              startcomp, endcomp) bind(C,name="update_rhoX_2d")

    integer         , intent(in   ) :: lo(2), hi(2)
    integer         , intent(in   ) :: so_lo(2), so_hi(2), nc_so
    double precision, intent(in   ) :: sold  (so_lo(1):so_hi(1),so_lo(2):so_hi(2),nc_so)
    integer         , intent(in   ) :: sn_lo(2), sn_hi(2), nc_sn
    double precision, intent(inout) :: snew  (sn_lo(1):sn_hi(1),sn_lo(2):sn_hi(2),nc_sn)
    integer         , intent(in   ) :: x_lo(2), x_hi(2), nc_x
    double precision, intent(in   ) :: sfluxx(x_lo(1):x_hi(1),x_lo(2):x_hi(2),nc_x)
    integer         , intent(in   ) :: y_lo(2), y_hi(2), nc_y
    double precision, intent(in   ) :: sfluxy(y_lo(1):y_hi(1),y_lo(2):y_hi(2),nc_y)
    integer         , intent(in   ) :: f_lo(2), f_hi(2), nc_f
    double precision, intent(in   ) :: force (f_lo(1):f_hi(1),f_lo(2):f_hi(2),nc_f)
    double precision, intent(in   ) :: dx(2), dt
    integer         , intent(in   ) :: startcomp, endcomp

    ! Local variables
    integer          :: i, j, comp, comp2
    double precision :: divterm 
    double precision :: delta, frac, sumX
    logical          :: has_negative_species
    ! from MAESTRO parameters file
    double precision :: base_cutoff_density = 3.d6

    do comp = startcomp, endcomp
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             
             divterm = (sfluxx(i+1,j,comp) - sfluxx(i,j,comp))/dx(1) &
                     + (sfluxy(i,j+1,comp) - sfluxy(i,j,comp))/dx(2)
             snew(i,j,comp) = sold(i,j,comp) + dt*(-divterm + force(i,j,comp))
             
          end do
       end do
    enddo

    ! update density
    snew(:,:,1) = sold(:,:,1)
           
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          has_negative_species = .false.
          
          ! define the update to rho as the sum of the updates to (rho X)_i  
          do comp = startcomp, endcomp
             snew(i,j,1) = snew(i,j,1) + (snew(i,j,comp)-sold(i,j,comp))
             if (snew(i,j,comp) .lt. ZERO) has_negative_species = .true.
          enddo

          ! enforce a density floor
          if (snew(i,j,1) .lt. HALF*base_cutoff_density) then
             do comp = startcomp, endcomp
                snew(i,j,comp) = snew(i,j,comp) * &
                     HALF*base_cutoff_density/snew(i,j,1)
             end do
             snew(i,j,1) = HALF*base_cutoff_density
          end if

          ! do not allow the species to leave here negative.
          if (has_negative_species) then
             do comp = startcomp, endcomp
                if (snew(i,j,comp) .lt. ZERO) then
                   delta = -snew(i,j,comp)
                   sumX = ZERO 
                   do comp2 = startcomp, endcomp
                      if (comp2 .ne. comp .and. snew(i,j,comp2) .ge. ZERO) then
                         sumX = sumX + snew(i,j,comp2)
                      end if
                   enddo
                   do comp2 = startcomp, endcomp
                      if (comp2 .ne. comp .and. snew(i,j,comp2) .ge. ZERO) then
                         frac = snew(i,j,comp2) / sumX
                         snew(i,j,comp2) = snew(i,j,comp2) - frac * delta
                      end if
                   enddo
                   snew(i,j,comp) = ZERO
                end if
             end do
          end if
       enddo
    enddo

  end subroutine update_rhoX_2d
#endif


#if (AMREX_SPACEDIM == 3)
  subroutine update_rhoX_3d(lo, hi, &
                              sold,   so_lo, so_hi, nc_so, &
                              snew,   sn_lo, sn_hi, nc_sn, &
                              sfluxx, x_lo, x_hi, nc_x, &
                              sfluxy, y_lo, y_hi, nc_y, &
                              sfluxz, z_lo, z_hi, nc_z, &
                              force,  f_lo, f_hi, nc_f, & 
                              dx, dt, &
                              startcomp, endcomp) bind(C,name="update_rhoX_3d")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: so_lo(3), so_hi(3), nc_so
    double precision, intent(in   ) :: sold  (so_lo(1):so_hi(1),so_lo(2):so_hi(2),so_lo(3):so_hi(3),nc_so)
    integer         , intent(in   ) :: sn_lo(3), sn_hi(3), nc_sn
    double precision, intent(inout) :: snew  (sn_lo(1):sn_hi(1),sn_lo(2):sn_hi(2),sn_lo(3):sn_hi(3),nc_sn)
    integer         , intent(in   ) :: x_lo(3), x_hi(3), nc_x
    double precision, intent(in   ) :: sfluxx(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3),nc_x)
    integer         , intent(in   ) :: y_lo(3), y_hi(3), nc_y
    double precision, intent(in   ) :: sfluxy(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3),nc_y)
    integer         , intent(in   ) :: z_lo(3), z_hi(3), nc_z
    double precision, intent(in   ) :: sfluxz(z_lo(1):z_hi(1),z_lo(2):z_hi(2),z_lo(3):z_hi(3),nc_z)
    integer         , intent(in   ) :: f_lo(3), f_hi(3), nc_f
    double precision, intent(in   ) :: force (f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),nc_f)
    double precision, intent(in   ) :: dx(3), dt
    integer         , intent(in   ) :: startcomp, endcomp

    ! Local variables
    integer          :: i, j, k, comp, comp2
    double precision :: divterm
    double precision :: delta, frac, sumX
    logical          :: has_negative_species
    ! from MAESTRO parameters file
    double precision :: base_cutoff_density = 3.d6
    
    !$OMP PARALLEL PRIVATE(i,j,k,divterm,comp) 
    do comp = startcomp, endcomp
       !$OMP DO
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                
                divterm = (sfluxx(i+1,j,k,comp) - sfluxx(i,j,k,comp))/dx(1) &
                        + (sfluxy(i,j+1,k,comp) - sfluxy(i,j,k,comp))/dx(2) &
                        + (sfluxz(i,j,k+1,comp) - sfluxz(i,j,k,comp))/dx(3)
                snew(i,j,k,comp) = sold(i,j,k,comp) + dt * (-divterm + force(i,j,k,comp))
                
             enddo
          enddo
       enddo
       !$OMP END DO NOWAIT
    end do
    !$OMP END PARALLEL
    
    ! update density
    snew(:,:,:,1) = sold(:,:,:,1)

    !$OMP PARALLEL DO PRIVATE(i,j,k,has_negative_species,comp,delta,sumX,comp2,frac)       
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             has_negative_species = .false.

             ! define the update to rho as the sum of the updates to (rho X)_i
             do comp = startcomp, endcomp
                snew(i,j,k,1) = snew(i,j,k,1) &
                     + (snew(i,j,k,comp)-sold(i,j,k,comp))
                if (snew(i,j,k,comp) .lt. ZERO) has_negative_species = .true.
             enddo

             ! enforce a density floor
             if (snew(i,j,k,1) .lt. HALF*base_cutoff_density) then
                do comp = startcomp, endcomp
                   snew(i,j,k,comp) = snew(i,j,k,comp) * &
                        HALF*base_cutoff_density/snew(i,j,k,1)
                end do
                snew(i,j,k,1) = HALF*base_cutoff_density
             end if

             ! do not allow the species to leave here negative.
             if (has_negative_species) then
                do comp = startcomp, endcomp
                   if (snew(i,j,k,comp) .lt. ZERO) then
                      delta = -snew(i,j,k,comp)
                      sumX = ZERO 
                      do comp2 = startcomp, endcomp
                         if (comp2 .ne. comp .and. snew(i,j,k,comp2) .ge. ZERO) then
                            sumX = sumX + snew(i,j,k,comp2)
                         end if
                      enddo
                      do comp2 = startcomp, endcomp
                         if (comp2 .ne. comp .and. snew(i,j,k,comp2) .ge. ZERO) then
                            frac = snew(i,j,k,comp2) / sumX
                            snew(i,j,k,comp2) = snew(i,j,k,comp2) - frac * delta
                         end if
                      enddo
                      snew(i,j,k,comp) = ZERO
                   end if
                end do
             end if         
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine update_rhoX_3d
#endif

end module update_scal_module

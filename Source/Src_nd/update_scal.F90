
module update_scal_module

  use amrex_constants_module
  use meth_params_module, only: rho_comp, rhoh_comp,spec_comp, temp_comp, &
                                  do_eos_h_above_cutoff, base_cutoff_density
  use base_state_geometry_module, only: nr_fine, max_radial_level
  use eos_module
  use eos_type_module

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

    do comp = startcomp, endcomp
       do i=lo(1),hi(1)
          divterm = (sfluxx(i+1,comp) - sfluxx(i,comp))/dx(1)
          snew(i,comp) = sold(i,comp) + dt*(-divterm + force(i,comp))
       end do
    enddo

    ! update density
    snew(lo(1):hi(1),rho_comp) = sold(lo(1):hi(1),rho_comp)

    do i = lo(1), hi(1)

       has_negative_species = .false.

       ! define the update to rho as the sum of the updates to (rho X)_i
       do comp = startcomp, endcomp
          snew(i,rho_comp) = snew(i,rho_comp) + (snew(i,comp)-sold(i,comp))
          if (snew(i,comp) .lt. ZERO) has_negative_species = .true.
       enddo

       ! enforce a density floor
       if (snew(i,rho_comp) .lt. HALF*base_cutoff_density) then
          do comp = startcomp, endcomp
             snew(i,comp) = snew(i,comp) * HALF*base_cutoff_density/snew(i,rho_comp)
          end do
          snew(i,rho_comp) = HALF*base_cutoff_density
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
    snew(lo(1):hi(1),lo(2):hi(2),rho_comp) = sold(lo(1):hi(1),lo(2):hi(2),rho_comp)

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          has_negative_species = .false.

          ! define the update to rho as the sum of the updates to (rho X)_i
          do comp = startcomp, endcomp
             snew(i,j,rho_comp) = snew(i,j,rho_comp) + (snew(i,j,comp)-sold(i,j,comp))
             if (snew(i,j,comp) .lt. ZERO) has_negative_species = .true.
          enddo

          ! enforce a density floor
          if (snew(i,j,rho_comp) .lt. HALF*base_cutoff_density) then
             do comp = startcomp, endcomp
                snew(i,j,comp) = snew(i,j,comp) * &
                     HALF*base_cutoff_density/snew(i,j,rho_comp)
             end do
             snew(i,j,rho_comp) = HALF*base_cutoff_density
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
    snew(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),rho_comp) = sold(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),rho_comp)

    !$OMP PARALLEL DO PRIVATE(i,j,k,has_negative_species,comp,delta,sumX,comp2,frac)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             has_negative_species = .false.

             ! define the update to rho as the sum of the updates to (rho X)_i
             do comp = startcomp, endcomp
                snew(i,j,k,rho_comp) = snew(i,j,k,rho_comp) &
                     + (snew(i,j,k,comp)-sold(i,j,k,comp))
                if (snew(i,j,k,comp) .lt. ZERO) has_negative_species = .true.
             enddo

             ! enforce a density floor
             if (snew(i,j,k,rho_comp) .lt. HALF*base_cutoff_density) then
                do comp = startcomp, endcomp
                   snew(i,j,k,comp) = snew(i,j,k,comp) * &
                        HALF*base_cutoff_density/snew(i,j,k,rho_comp)
                end do
                snew(i,j,k,rho_comp) = HALF*base_cutoff_density
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


#if (AMREX_SPACEDIM == 1)
  subroutine update_rhoh_1d(lev, lo, hi, &
                              sold,   so_lo, so_hi, nc_so, &
                              snew,   sn_lo, sn_hi, nc_sn, &
                              sfluxx, x_lo, x_hi, nc_x, &
                              force,  f_lo, f_hi, nc_f, &
                              p0_new, &
                              dx, dt, &
                              nspec) bind(C,name="update_rhoh_1d")

    integer         , intent(in   ) :: lev, lo(1), hi(1)
    integer         , intent(in   ) :: so_lo(1), so_hi(1), nc_so
    double precision, intent(in   ) :: sold  (so_lo(1):so_hi(1),nc_so)
    integer         , intent(in   ) :: sn_lo(1), sn_hi(1), nc_sn
    double precision, intent(inout) :: snew  (sn_lo(1):sn_hi(1),nc_sn)
    integer         , intent(in   ) :: x_lo(1), x_hi(1), nc_x
    double precision, intent(in   ) :: sfluxx(x_lo(1):x_hi(1),nc_x)
    integer         , intent(in   ) :: f_lo(1), f_hi(1), nc_f
    double precision, intent(in   ) :: force (f_lo(1):f_hi(1),nc_f)
    double precision, intent(in   ) :: p0_new(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: dx(1), dt
    integer         , intent(in   ) :: nspec

    ! Local variables
    integer          :: i
    double precision :: divterm

    integer :: pt_index(3)
    type(eos_t)      :: eos_state

    do i=lo(1),hi(1)
       divterm = (sfluxx(i+1,rhoh_comp) - sfluxx(i,rhoh_comp))/dx(1)
       snew(i,rhoh_comp) = sold(i,rhoh_comp) + dt*(-divterm + force(i,rhoh_comp))
    end do

    if (do_eos_h_above_cutoff) then
       do i = lo(1), hi(1)
          if (snew(i,rho_comp) .le. base_cutoff_density) then

             eos_state%rho = snew(i,rho_comp)
             eos_state%T   = sold(i,temp_comp)
             eos_state%p   = p0_new(lev,i)
             eos_state%xn  = snew(i,spec_comp:spec_comp+nspec-1)/eos_state%rho

             pt_index(:) = (/i, -1, -1/)

             ! (rho,P) --> T,h
             call eos(eos_input_rp, eos_state, pt_index)

             snew(i,rhoh_comp) = snew(i,rho_comp) * eos_state%h

          end if
       enddo
    end if

  end subroutine update_rhoh_1d
#endif

#if (AMREX_SPACEDIM == 2)
  subroutine update_rhoh_2d(lev, lo, hi, &
                              sold,   so_lo, so_hi, nc_so, &
                              snew,   sn_lo, sn_hi, nc_sn, &
                              sfluxx, x_lo, x_hi, nc_x, &
                              sfluxy, y_lo, y_hi, nc_y, &
                              force,  f_lo, f_hi, nc_f, &
                              p0_new, &
                              dx, dt, &
                              nspec) bind(C,name="update_rhoh_2d")

    integer         , intent(in   ) :: lev, lo(2), hi(2)
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
    double precision, intent(in   ) :: p0_new(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: dx(2), dt
    integer         , intent(in   ) :: nspec

    ! Local variables
    integer          :: i, j
    double precision :: divterm

    integer :: pt_index(3)
    type(eos_t) :: eos_state

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          divterm = (sfluxx(i+1,j,rhoh_comp) - sfluxx(i,j,rhoh_comp))/dx(1) &
               + (sfluxy(i,j+1,rhoh_comp) - sfluxy(i,j,rhoh_comp))/dx(2)
          snew(i,j,rhoh_comp) = sold(i,j,rhoh_comp) + dt*(-divterm + force(i,j,rhoh_comp))

       end do
    end do

    if (do_eos_h_above_cutoff) then
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (snew(i,j,rho_comp) .le. base_cutoff_density) then
                eos_state%rho = snew(i,j,rho_comp)
                eos_state%T   = sold(i,j,temp_comp)
                eos_state%p   = p0_new(lev,j)
                eos_state%xn  = snew(i,j,spec_comp:spec_comp+nspec-1)/eos_state%rho

                pt_index(:) = (/i, j, -1/)

                ! (rho,P) --> T,h
                call eos(eos_input_rp, eos_state, pt_index)

                snew(i,j,rhoh_comp) = snew(i,j,rho_comp) * eos_state%h
             end if

          enddo
       enddo
    end if

  end subroutine update_rhoh_2d
#endif


#if (AMREX_SPACEDIM == 3)
  subroutine update_rhoh_3d(lev, lo, hi, &
                              sold,   so_lo, so_hi, nc_so, &
                              snew,   sn_lo, sn_hi, nc_sn, &
                              sfluxx, x_lo, x_hi, nc_x, &
                              sfluxy, y_lo, y_hi, nc_y, &
                              sfluxz, z_lo, z_hi, nc_z, &
                              force,  f_lo, f_hi, nc_f, &
                              p0_new, &
                              dx, dt, &
                              nspec) bind(C,name="update_rhoh_3d")

    integer         , intent(in   ) :: lev, lo(3), hi(3)
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
    double precision, intent(in   ) :: p0_new(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: dx(3), dt
    integer         , intent(in   ) :: nspec

    ! Local variables
    integer          :: i, j, k
    double precision :: divterm

    integer :: pt_index(3)
    type(eos_t) :: eos_state

    !$OMP PARALLEL PRIVATE(i,j,k,divterm)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             divterm = (sfluxx(i+1,j,k,rhoh_comp) - sfluxx(i,j,k,rhoh_comp))/dx(1) &
                  + (sfluxy(i,j+1,k,rhoh_comp) - sfluxy(i,j,k,rhoh_comp))/dx(2) &
                  + (sfluxz(i,j,k+1,rhoh_comp) - sfluxz(i,j,k,rhoh_comp))/dx(3)
             snew(i,j,k,rhoh_comp) = sold(i,j,k,rhoh_comp) + dt * (-divterm + force(i,j,k,rhoh_comp))

          enddo
       enddo
    enddo
    !$OMP END PARALLEL

    if ( do_eos_h_above_cutoff ) then
       !$OMP PARALLEL DO PRIVATE(i,j,k,eos_state,pt_index)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                if (snew(i,j,k,rho_comp) .le. base_cutoff_density) then

                   eos_state%rho = snew(i,j,k,rho_comp)
                   eos_state%T   = sold(i,j,k,temp_comp)
                   eos_state%p   = p0_new(lev,k)
                   eos_state%xn  = snew(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

                   pt_index(:) = (/i, j, k/)

                   ! (rho,P) --> T,h
                   call eos(eos_input_rp, eos_state, pt_index)

                   snew(i,j,k,rhoh_comp) = snew(i,j,k,rho_comp) * eos_state%h

                end if

             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
    end if

  end subroutine update_rhoh_3d

  subroutine update_rhoh_3d_sphr(lo, hi, &
                                   sold,   so_lo, so_hi, nc_so, &
                                   snew,   sn_lo, sn_hi, nc_sn, &
                                   sfluxx, x_lo, x_hi, nc_x, &
                                   sfluxy, y_lo, y_hi, nc_y, &
                                   sfluxz, z_lo, z_hi, nc_z, &
                                   force,  f_lo, f_hi, nc_f, &
                                   p0_new_cart, p_lo, p_hi, &
                                   dx, dt, &
                                   nspec) bind(C,name="update_rhoh_3d_sphr")

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
    integer         , intent(in   ) :: p_lo(3), p_hi(3)
    double precision, intent(in   ) :: p0_new_cart(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
    double precision, intent(in   ) :: dx(3), dt
    integer         , intent(in   ) :: nspec

    ! Local variables
    integer          :: i, j, k
    double precision :: divterm

    integer :: pt_index(3)
    type(eos_t) :: eos_state

    !$OMP PARALLEL PRIVATE(i,j,k,divterm)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             divterm = (sfluxx(i+1,j,k,rhoh_comp) - sfluxx(i,j,k,rhoh_comp))/dx(1) &
                  + (sfluxy(i,j+1,k,rhoh_comp) - sfluxy(i,j,k,rhoh_comp))/dx(2) &
                  + (sfluxz(i,j,k+1,rhoh_comp) - sfluxz(i,j,k,rhoh_comp))/dx(3)
             snew(i,j,k,rhoh_comp) = sold(i,j,k,rhoh_comp) + dt * (-divterm + force(i,j,k,rhoh_comp))

          enddo
       enddo
    enddo
    !$OMP END PARALLEL

    if ( do_eos_h_above_cutoff ) then
       !$OMP PARALLEL DO PRIVATE(i,j,k,eos_state,pt_index)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                if (snew(i,j,k,rho_comp) .le. base_cutoff_density) then

                   eos_state%rho = snew(i,j,k,rho_comp)
                   eos_state%T   = sold(i,j,k,temp_comp)
                   eos_state%p   = p0_new_cart(i,j,k)
                   eos_state%xn  = snew(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

                   pt_index(:) = (/i, j, k/)

                   ! (rho,P) --> T,h
                   call eos(eos_input_rp, eos_state, pt_index)

                   snew(i,j,k,rhoh_comp) = snew(i,j,k,rho_comp) * eos_state%h

                end if

             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
    end if

  end subroutine update_rhoh_3d_sphr
#endif

end module update_scal_module

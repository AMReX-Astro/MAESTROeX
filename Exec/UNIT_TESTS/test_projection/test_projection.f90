#include "AMReX_BC_TYPES.H"

module test_projection_module

  use bl_types
  use bl_constants_module
  use bl_error
  use probin_module, only: project_type, prob_lo

  implicit none

  private

contains

  !===========================================================================
  subroutine init_mac_velocity(umac, dx, mla, the_bc_level)

    integer :: n, i, ng, dm, nlevs

    type(multifab) , intent(inout) :: umac(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(ml_layout)   , intent(inout) :: mla
    type(bc_level)    , intent(in   ) :: the_bc_level(:)

    integer :: lo(get_dim(umac(1,1))), hi(get_dim(umac(1,1)))

    real(kind=dp_t), pointer :: ump(:,:,:,:), vmp(:,:,:,:), wmp(:,:,:,:)

    nlevs = size(umac(:,1))
    dm = get_dim(umac(1,1))

    ng = nghost(umac(1,1))

    do n=1,nlevs
       do i = 1, nfabs(umac(n,1))
          ump => dataptr(umac(n,1), i)
          vmp => dataptr(umac(n,2), i)

          lo = lwb(get_box(umac(n,1), i))
          hi = upb(get_box(umac(n,1), i))

          select case (dm)
          case (2)
             call init_mac_velocity_2d(ump(:,:,1,1), vmp(:,:,1,1), ng, &
                  lo, hi, dx(n,:))

          case (3)
             wmp => dataptr(umac(n,3), i)
             call init_mac_velocity_3d(ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), ng, &
                  lo, hi, dx(n,:))

          end select
       end do
    end do

    ! make edge states consistent across levels
    do n = nlevs,2,-1
       do i = 1, dm
          call ml_edge_restriction_c(umac(n-1,i),1,umac(n,i),1,mla%mba%rr(n-1,:),i,1)
       enddo
    enddo

  end subroutine init_mac_velocity

  subroutine init_mac_velocity_2d(umac, vmac, ng, lo, hi, dx)

    ! initialize the velocity field to a divergence-free field.  This
    ! velocity field comes from Almgren, Bell, and Szymczak 1996.

    use probin_module, only: prob_lo

    integer         , intent(in   ) :: lo(:), hi(:), ng
    real (kind=dp_t), intent(inout) :: umac(lo(1)-ng:,lo(2)-ng:)
    real (kind=dp_t), intent(inout) :: vmac(lo(1)-ng:,lo(2)-ng:)
    real (kind=dp_t), intent(in   ) :: dx(:)

    ! Local variables
    integer :: i, j
    real (kind=dp_t) :: x, y

    ! x-velocity  (x are edges, y are centers)
    do j = lo(2), hi(2)
       y = (dble(j)+HALF)*dx(2) + prob_lo(2)

       do i = lo(1), hi(1)+1
          x = (dble(i))*dx(1) + prob_lo(1)

          umac(i,j) = -sin(M_PI*x)**2 * sin(TWO*M_PI*y)

       enddo
    enddo

    ! y-velocity  (x are centers, y are edges)
    do j = lo(2), hi(2)+1
       y = (dble(j))*dx(2) + prob_lo(2)

       do i = lo(1), hi(1)
          x = (dble(i)+HALF)*dx(1) + prob_lo(1)

          vmac(i,j) =  sin(M_PI*y)**2 * sin(TWO*M_PI*x)

       enddo
    enddo

  end subroutine init_mac_velocity_2d


  subroutine init_mac_velocity_3d(umac, vmac, wmac, ng, lo, hi, dx)

    ! initialize the velocity field to a divergence-free field.  This
    ! velocity field comes from the idea that the curl of any vector
    ! is divergence-free.
    !
    ! we take: Phi = (alpha, beta, gamma) as our initial vector, with:
    !
    !   alpha = sin 4pi y  sin 2pi z
    !   beta  = sin 2pi x  sin 4pi z
    !   gamma = sin 4pi x  sin 2pi y
    !
    ! (this is Interior and even in all directions)
    !
    ! then U = curl{Phi} gives our field

    use probin_module, only: prob_lo

    integer         , intent(in   ) :: lo(:), hi(:), ng
    real (kind=dp_t), intent(inout) :: umac(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real (kind=dp_t), intent(inout) :: vmac(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real (kind=dp_t), intent(inout) :: wmac(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real (kind=dp_t), intent(in   ) :: dx(:)

    ! Local variables
    integer :: i, j, k
    real (kind=dp_t) :: x, y, z

    ! x-velocity  (x are edges, y and z are centers)
    do k = lo(3), hi(3)
       z = (dble(k)+HALF)*dx(3) + prob_lo(3)

       do j = lo(2), hi(2)
          y = (dble(j)+HALF)*dx(2) + prob_lo(2)

          do i = lo(1), hi(1)+1
             x = (dble(i))*dx(1) + prob_lo(1)

             umac(i,j,k) = TWO*M_PI*sin(FOUR*M_PI*x)*cos( TWO*M_PI*y) - &
                  FOUR*M_PI*sin( TWO*M_PI*x)*cos(FOUR*M_PI*z)

          enddo
       enddo
    enddo

    ! y-velocity  (x and z are centers, y are edges)
    do k = lo(3), hi(3)
       z = (dble(k)+HALF)*dx(3) + prob_lo(3)

       do j = lo(2), hi(2)+1
          y = (dble(j))*dx(2) + prob_lo(2)

          do i = lo(1), hi(1)
             x = (dble(i)+HALF)*dx(1) + prob_lo(1)

             vmac(i,j,k) = TWO*M_PI*sin(FOUR*M_PI*y)*cos( TWO*M_PI*z) - &
                  FOUR*M_PI*cos(FOUR*M_PI*x)*sin( TWO*M_PI*y)

          enddo
       enddo
    enddo

    ! z-velocity  (x and y are centers, z are edges)
    do k = lo(3), hi(3)+1
       z = (dble(k))*dx(3) + prob_lo(3)

       do j = lo(2), hi(2)
          y = (dble(j)+HALF)*dx(2) + prob_lo(2)

          do i = lo(1), hi(1)
             x = (dble(i)+HALF)*dx(1) + prob_lo(1)

             wmac(i,j,k) = TWO*M_PI*cos( TWO*M_PI*x)*sin(FOUR*M_PI*z) - &
                  FOUR*M_PI*cos(FOUR*M_PI*y)*sin( TWO*M_PI*z)


          enddo
       enddo
    enddo

  end subroutine init_mac_velocity_3d


  !===========================================================================
  subroutine add_grad_scalar(lo, hi, gphi, g_lo, g_hi, nc_g, &
       U, u_lo, u_hi, nc_v, domain_phys_bc, dx)  bind(C, name="add_grad_scalar")

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: g_lo(3), g_hi(3), nc_g
    integer, intent(in) :: v_lo(3), v_hi(3), nc_v
    integer, intent(in) :: domain_phys_bc(amrex_spacedim,2)
    double precision, intent(inout) :: gphi(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3),1:nc_g)
    double precision, intent(inout) :: U(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),1:nc_u)
    double precision, intent(in) :: dx(3)

    integer :: i, j, k
    double precision :: x, y, z
    double precision :: phi(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

    if (domain_phys_bc(1,1) .eq. SlipWall .and. &
         domain_phys_bc(1,2) .eq. SlipWall .and. &
         domain_phys_bc(2,1) .eq. SlipWall .and. &
         domain_phys_bc(2,2) .eq. SlipWall .and. &
         domain_phys_bc(3,1) .eq. SlipWall .and. &
         domain_phys_bc(3,2) .eq. SlipWall) then

       do k = lo(3), hi(3)
          z = (dble(k)+0.5d0)*dx(3) + prob_lo(3)

          do j = lo(2), hi(2)
             y = (dble(j)+0.5d0)*dx(2) + prob_lo(2)

             do i = lo(1), hi(1)
                x = (dble(i)+0.5d0)*dx(1) + prob_lo(1)

                gphi(i,j,k,1) =  160.0d0*x*(ONE - x)
                gphi(i,j,k,2) =  160.0d0*y*(ONE - y)
                gphi(i,j,k,3) =  160.0d0*z*(ONE - z)

                U(i,j,k,1) = U(i,j,k,1) + gphi(i,j,k,1)
                U(i,j,k,2) = U(i,j,k,2) + gphi(i,j,k,2)
                U(i,j,k,3) = U(i,j,k,3) + gphi(i,j,k,3)

             enddo
          enddo
       enddo

    else if (domain_phys_bc(1,1) .eq. Interior .and. &
         domain_phys_bc(1,2) .eq. Interior .and. &
         domain_phys_bc(2,1) .eq. Interior .and. &
         domain_phys_bc(2,2) .eq. Interior .and. &
         domain_phys_bc(3,1) .eq. Interior .and. &
         domain_phys_bc(3,2) .eq. Interior) then

       do k = lo(3), hi(3)
          z = (dble(k)+0.5d0)*dx(3) + prob_lo(3)

          do j = lo(2), hi(2)
             y = (dble(j)+0.5d0)*dx(2) + prob_lo(2)

             do i = lo(1), hi(1)
                x = (dble(i)+0.5d0)*dx(1) + prob_lo(1)

                phi(i,j,k) = 5.0d0 * cos(2.d0*M_PI*y)*cos(2.d0*M_PI*x)*cos(2.d0*M_PI*z)

             enddo
          enddo
       enddo

       do k = lo(3), hi(3)
          z = (dble(k)+0.5d0)*dx(3) + prob_lo(3)

          do j = lo(2), hi(2)
             y = (dble(j)+0.5d0)*dx(2) + prob_lo(2)

             do i = lo(1), hi(1)
                x = (dble(i)+0.5d0)*dx(1) + prob_lo(1)

                gphi(i,j,k,1) =  (phi(i+1,j,k)-phi(i-1,j,k))/(2.d0*dx(1))
                gphi(i,j,k,2) =  (phi(i,j+1,k)-phi(i,j-1,k))/(2.d0*dx(2))
                gphi(i,j,k,3) =  (phi(i,j,k+1)-phi(i,j,k-1))/(2.d0*dx(3))

                U(i,j,k,1) = U(i,j,k,1) + gphi(i,j,k,1)
                U(i,j,k,2) = U(i,j,k,2) + gphi(i,j,k,2)
                U(i,j,k,3) = U(i,j,k,3) + gphi(i,j,k,3)

             enddo
          enddo
       enddo

    else
       call bl_error('Not set up for these boundary conditions')

    endif

  end subroutine add_grad_scalar


  !===========================================================================
  subroutine add_grad_scalar_mac(lo, hi, gphi_mac, g_lo, g_hi, nc_g, &
       umac, u_lo, u_hi, nc_u, domain_phys_bc, box_phys_bc, dx)

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: g_lo(3), g_hi(3), nc_g
    integer, intent(in) :: v_lo(3), v_hi(3), nc_v
    integer, intent(in) :: domain_phys_bc(amrex_spacedim,2)
    integer, intent(in) :: box_phys_bc(amrex_spacedim,2)
    double precision, intent(inout) :: gphi_mac(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3),1:nc_g)
    double precision, intent(inout) :: umac(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),1:nc_u)
    double precision, intent(in) :: dx(3)

    integer :: i, j, k
    double precision :: x, y, z
    double precision :: phi(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

    if (domain_phys_bc(1,1) .eq. SlipWall .and. &
         domain_phys_bc(1,2) .eq. SlipWall .and. &
         domain_phys_bc(2,1) .eq. SlipWall .and. &
         domain_phys_bc(2,2) .eq. SlipWall .and. &
         domain_phys_bc(3,1) .eq. SlipWall .and. &
         domain_phys_bc(3,2) .eq. SlipWall) then

       ! Add on the gradient of a scalar (phi) that satisfies
       ! grad(phi).n = 0.

       ! x-velocity  (x are edges, y and z are centers)
       do k = lo(3), hi(3)
          z = (dble(k)+HALF)*dx(3) + prob_lo(3)

          do j = lo(2), hi(2)
             y = (dble(j)+HALF)*dx(2) + prob_lo(2)

             do i = lo(1), hi(1)+1
                x = (dble(i))*dx(1) + prob_lo(1)

                gphi(i,j,k,1) =  160.0d0*x*(ONE - x)
                umac(i,j,k,1) = umac(i,j,k,1) + gphi(i,j,k,1)

             enddo
          enddo
       enddo

       ! y-velocity  (x and z are centers, y are edges)
       do k = lo(3), hi(3)
          z = (dble(k)+HALF)*dx(3) + prob_lo(3)

          do j = lo(2), hi(2)+1
             y = (dble(j))*dx(2) + prob_lo(2)

             do i = lo(1), hi(1)
                x = (dble(i)+HALF)*dx(1) + prob_lo(1)

                gphi(i,j,k,2) =  160.0d0*y*(ONE - y)
                umac(i,j,k,2) = umac(i,j,k,2) + gphi(i,j,k,2)

             enddo
          enddo
       enddo

       ! z-velocity  (x and y are centers, z are edges)
       do k = lo(3), hi(3)+1
          z = (dble(k))*dx(3) + prob_lo(3)

          do j = lo(2), hi(2)
             y = (dble(j)+HALF)*dx(2) + prob_lo(2)

             do i = lo(1), hi(1)
                x = (dble(i)+HALF)*dx(1) + prob_lo(1)

                gphi(i,j,k,3) =  160.0d0*z*(ONE - z)
                umac(i,j,k,3) = umac(i,j,k,3) + gphi(i,j,k,3)

             enddo
          enddo
       enddo


    else if (domain_phys_bc(1,1) .eq. Interior .and. &
         domain_phys_bc(1,2) .eq. Interior .and. &
         domain_phys_bc(2,1) .eq. Interior .and. &
         domain_phys_bc(2,2) .eq. Interior .and. &
         domain_phys_bc(3,1) .eq. Interior .and. &
         domain_phys_bc(3,2) .eq. Interior) then


       do k = lo(3)-1, hi(3)+1
          z = (dble(k)+HALF)*dx(3) + prob_lo(3)

          do j = lo(2)-1, hi(2)+1
             y = (dble(j)+HALF)*dx(2) + prob_lo(2)

             do i = lo(1)-1, hi(1)+1
                x = (dble(i)+HALF)*dx(1) + prob_lo(1)

                phi(i,j,k) = 5.0d0 * cos(2.d0*M_PI*y)*cos(2.d0*M_PI*x)*cos(2.d0*M_PI*z)
             enddo
          enddo
       enddo

       ! x-velocity  (x are edges, y and z are centers)
       do k = lo(3), hi(3)
          z = (dble(k)+HALF)*dx(3) + prob_lo(3)

          do j = lo(2), hi(2)
             y = (dble(j)+HALF)*dx(2) + prob_lo(2)

             do i = lo(1), hi(1)+1
                x = (dble(i))*dx(1) + prob_lo(1)

                gphi(i,j,k,1) = (phi(i,j,k) - phi(i-1,j,k))/dx(1)
                umac(i,j,k,1) = umac(i,j,k,1) + gphi(i,j,k,1)

             enddo
          enddo
       enddo

       ! y-velocity  (x and z are centers, y are edges)
       do k = lo(3), hi(3)
          z = (dble(k)+HALF)*dx(3) + prob_lo(3)

          do j = lo(2), hi(2)+1
             y = (dble(j))*dx(2) + prob_lo(2)

             do i = lo(1), hi(1)
                x = (dble(i)+HALF)*dx(1) + prob_lo(1)

                gphi(i,j,k,2) = (phi(i,j,k) - phi(i,j-1,k))/dx(2)
                umac(i,j,k,2) = umac(i,j,k,2) + gphi(i,j,k,2)

             enddo
          enddo
       enddo

       ! z-velocity  (x and y are centers, z are edges)
       do k = lo(3), hi(3)+1
          z = (dble(k))*dx(3) + prob_lo(3)

          do j = lo(2), hi(2)
             y = (dble(j)+HALF)*dx(2) + prob_lo(2)

             do i = lo(1), hi(1)
                x = (dble(i)+HALF)*dx(1) + prob_lo(1)

                gphi(i,j,k,3) = (phi(i,j,k) - phi(i,j,k-1))/dx(3)
                umac(i,j,k,3) = umac(i,j,k,3) + gphi(i,j,k,3)

             enddo
          enddo
       enddo

    else
       call bl_error('Not set up for these boundary conditions')
    endif

    ! impose BCs

    ! x lo
    select case (box_phys_bc(1,1))
    case (SlipWall)
       umac(lo(1),lo(2):hi(2),lo(3):hi(3),1) = ZERO
    case (Interior, INTERIOR)
    case default
       call bl_error("invalid x lo BC")
    end select

    ! x hi
    select case(box_phys_bc(1,2))
    case (SlipWall)
       umac(hi(1)+1,lo(2):hi(2),lo(3):hi(3),1) = ZERO
    case (Interior, INTERIOR)
    case default
       call bl_error("invalid x hi BC")
    end select

    ! y lo
    select case (box_phys_bc(2,1))
    case (SlipWall)
       umac(lo(1):hi(1),lo(2),lo(3):hi(3),2) = ZERO
    case (Interior, INTERIOR)
    case default
       call bl_error("invalid y lo BC")
    end select

    ! y hi
    select case(box_phys_bc(2,2))
    case (SlipWall)
       umac(lo(1):hi(1),hi(2)+1,lo(3):hi(3),2) = ZERO
    case (Interior, INTERIOR)
    case default
       call bl_error("invalid y hi BC")
    end select

    ! z lo
    select case (box_phys_bc(3,1))
    case (SlipWall)
       umac(lo(1):hi(1),lo(2):hi(2),lo(3),3) = ZERO
    case (Interior, INTERIOR)
    case default
       call bl_error("invalid z lo BC")
    end select

    ! z hi
    select case(box_phys_bc(3,2))
    case (SlipWall)
       umac(lo(1):hi(1),lo(2):hi(2),hi(3)+1,3) = ZERO
    case (Interior, INTERIOR)
    case default
       call bl_error("invalid z hi BC")
    end select


  end subroutine add_grad_scalar_mac

  subroutine add_grad_scalar_2d_mac(umac, vmac, ng_um, &
       gphi_x, gphi_y, ng_gp, &
       lo, hi, dx, &
       domain_phys_bc, box_phys_bc)

    use     bc_module
    use probin_module, only: prob_lo

    integer         , intent(in   ) :: lo(:), hi(:), ng_um, ng_gp
    real (kind=dp_t), intent(inout) ::   umac(lo(1)-ng_um:,lo(2)-ng_um:)
    real (kind=dp_t), intent(inout) ::   vmac(lo(1)-ng_um:,lo(2)-ng_um:)
    real (kind=dp_t), intent(inout) :: gphi_x(lo(1)-ng_gp:,lo(2)-ng_gp:)
    real (kind=dp_t), intent(inout) :: gphi_y(lo(1)-ng_gp:,lo(2)-ng_gp:)
    real (kind=dp_t), intent(in   ) :: dx(:)
    integer         , intent(in   ) :: domain_phys_bc(:,:)
    integer         , intent(in   ) :: box_phys_bc(:,:)


    ! Local variables
    integer :: i, j
    real (kind=dp_t) :: x, y
    real (kind=dp_t), allocatable :: phi(:,:)

    allocate(phi(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1))

    if (domain_phys_bc(1,1) .eq. SlipWall .and. &
         domain_phys_bc(1,2) .eq. SlipWall .and. &
         domain_phys_bc(2,1) .eq. SlipWall .and. &
         domain_phys_bc(2,2) .eq. SlipWall) then

       ! Add on the gradient of a scalar (phi) that satisfies
       ! grad(phi).n = 0.

       ! x-velocity  (x are edges, y are centers)
       do j = lo(2), hi(2)
          y = (dble(j)+HALF)*dx(2) + prob_lo(2)

          do i = lo(1), hi(1)+1
             x = (dble(i))*dx(1) + prob_lo(1)

             gphi_x(i,j) =  FOUR*x*(ONE - x)
             umac(i,j) = umac(i,j) + gphi_x(i,j)

          enddo
       enddo

       ! y-velocity  (x are centers, y are edges)
       do j = lo(2), hi(2)+1
          y = (dble(j))*dx(2) + prob_lo(2)

          do i = lo(1), hi(1)
             x = (dble(i)+HALF)*dx(1) + prob_lo(1)

             gphi_y(i,j) =  FOUR*y*(ONE - y)
             vmac(i,j) = vmac(i,j) + gphi_y(i,j)

          enddo
       enddo


    else if (domain_phys_bc(1,1) .eq. Interior .and. &
         domain_phys_bc(1,2) .eq. Interior .and. &
         domain_phys_bc(2,1) .eq. Interior .and. &
         domain_phys_bc(2,2) .eq. Interior) then

       do j = lo(2)-1, hi(2)+1
          y = (dble(j)+0.5d0)*dx(2) + prob_lo(2)

          do i = lo(1)-1, hi(1)+1
             x = (dble(i)+0.5d0)*dx(1) + prob_lo(1)

             phi(i,j) = 0.1d0 * cos(2.d0*M_PI*y)*cos(2.d0*M_PI*x)

          enddo
       enddo

       ! x-velocity  (x are edges, y are centers)
       do j = lo(2), hi(2)
          y = (dble(j)+HALF)*dx(2) + prob_lo(2)

          do i = lo(1), hi(1)+1
             x = (dble(i))*dx(1) + prob_lo(1)

             gphi_x(i,j) = (phi(i,j) - phi(i-1,j))/dx(1)
             umac(i,j) = umac(i,j) + gphi_x(i,j)

          enddo
       enddo

       ! y-velocity  (x are centers, y are edges)
       do j = lo(2), hi(2)+1
          y = (dble(j))*dx(2) + prob_lo(2)

          do i = lo(1), hi(1)
             x = (dble(i)+HALF)*dx(1) + prob_lo(1)

             gphi_y(i,j) = (phi(i,j) - phi(i,j-1))/dx(2)
             vmac(i,j) = vmac(i,j) + gphi_y(i,j)

          enddo
       enddo

    else
       call bl_error('Not set up for these boundary conditions')
    end if


    ! impose BCs

    ! x lo
    select case (box_phys_bc(1,1))
    case (SlipWall)
       umac(lo(1),lo(2):hi(2)) = ZERO
    case (Interior, INTERIOR)
    case default
       call bl_error("invalid x lo BC")
    end select

    ! x hi
    select case(box_phys_bc(1,2))
    case (SlipWall)
       umac(hi(1)+1,lo(2):hi(2)) = ZERO
    case (Interior, INTERIOR)
    case default
       call bl_error("invalid x hi BC")
    end select

    ! y lo
    select case (box_phys_bc(2,1))
    case (SlipWall)
       vmac(lo(1):hi(1),lo(2)) = ZERO
    case (Interior, INTERIOR)
    case default
       call bl_error("invalid y lo BC")
    end select

    ! y hi
    select case(box_phys_bc(2,2))
    case (SlipWall)
       vmac(lo(1):hi(1),hi(2)+1) = ZERO
    case (Interior, INTERIOR)
    case default
       call bl_error("invalid y hi BC")
    end select

    deallocate(phi)

  end subroutine add_grad_scalar_2d_mac


  subroutine add_grad_scalar_3d_mac(umac, vmac, wmac, ng_um, &
       gphi_x, gphi_y, gphi_z, ng_gp, &
       lo, hi, dx, &
       domain_phys_bc, box_phys_bc)

    use     bc_module
    use probin_module, only: prob_lo

    integer         , intent(in   ) :: lo(:), hi(:), ng_um, ng_gp
    real (kind=dp_t), intent(inout) ::   umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real (kind=dp_t), intent(inout) ::   vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real (kind=dp_t), intent(inout) ::   wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real (kind=dp_t), intent(inout) :: gphi_x(lo(1)-ng_gp:,lo(2)-ng_gp:,lo(3)-ng_gp:)
    real (kind=dp_t), intent(inout) :: gphi_y(lo(1)-ng_gp:,lo(2)-ng_gp:,lo(3)-ng_gp:)
    real (kind=dp_t), intent(inout) :: gphi_z(lo(1)-ng_gp:,lo(2)-ng_gp:,lo(3)-ng_gp:)
    real (kind=dp_t), intent(in   ) :: dx(:)
    integer         , intent(in   ) :: domain_phys_bc(:,:)
    integer         , intent(in   ) :: box_phys_bc(:,:)


    ! Local variables
    integer :: i, j, k
    real (kind=dp_t) :: x, y, z
    real (kind=dp_t), allocatable :: phi(:,:,:)

    allocate(phi(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))

    if (domain_phys_bc(1,1) .eq. SlipWall .and. &
         domain_phys_bc(1,2) .eq. SlipWall .and. &
         domain_phys_bc(2,1) .eq. SlipWall .and. &
         domain_phys_bc(2,2) .eq. SlipWall .and. &
         domain_phys_bc(3,1) .eq. SlipWall .and. &
         domain_phys_bc(3,2) .eq. SlipWall) then

       ! Add on the gradient of a scalar (phi) that satisfies
       ! grad(phi).n = 0.

       ! x-velocity  (x are edges, y and z are centers)
       do k = lo(3), hi(3)
          z = (dble(k)+HALF)*dx(3) + prob_lo(3)

          do j = lo(2), hi(2)
             y = (dble(j)+HALF)*dx(2) + prob_lo(2)

             do i = lo(1), hi(1)+1
                x = (dble(i))*dx(1) + prob_lo(1)

                gphi_x(i,j,k) =  160.0d0*x*(ONE - x)
                umac(i,j,k) = umac(i,j,k) + gphi_x(i,j,k)

             enddo
          enddo
       enddo

       ! y-velocity  (x and z are centers, y are edges)
       do k = lo(3), hi(3)
          z = (dble(k)+HALF)*dx(3) + prob_lo(3)

          do j = lo(2), hi(2)+1
             y = (dble(j))*dx(2) + prob_lo(2)

             do i = lo(1), hi(1)
                x = (dble(i)+HALF)*dx(1) + prob_lo(1)

                gphi_y(i,j,k) =  160.0d0*y*(ONE - y)
                vmac(i,j,k) = vmac(i,j,k) + gphi_y(i,j,k)

             enddo
          enddo
       enddo

       ! z-velocity  (x and y are centers, z are edges)
       do k = lo(3), hi(3)+1
          z = (dble(k))*dx(3) + prob_lo(3)

          do j = lo(2), hi(2)
             y = (dble(j)+HALF)*dx(2) + prob_lo(2)

             do i = lo(1), hi(1)
                x = (dble(i)+HALF)*dx(1) + prob_lo(1)

                gphi_z(i,j,k) =  160.0d0*z*(ONE - z)
                wmac(i,j,k) = wmac(i,j,k) + gphi_z(i,j,k)

             enddo
          enddo
       enddo


    else if (domain_phys_bc(1,1) .eq. Interior .and. &
         domain_phys_bc(1,2) .eq. Interior .and. &
         domain_phys_bc(2,1) .eq. Interior .and. &
         domain_phys_bc(2,2) .eq. Interior .and. &
         domain_phys_bc(3,1) .eq. Interior .and. &
         domain_phys_bc(3,2) .eq. Interior) then


       do k = lo(3)-1, hi(3)+1
          z = (dble(k)+HALF)*dx(3) + prob_lo(3)

          do j = lo(2)-1, hi(2)+1
             y = (dble(j)+HALF)*dx(2) + prob_lo(2)

             do i = lo(1)-1, hi(1)+1
                x = (dble(i)+HALF)*dx(1) + prob_lo(1)

                phi(i,j,k) = 5.0d0 * cos(2.d0*M_PI*y)*cos(2.d0*M_PI*x)*cos(2.d0*M_PI*z)
             enddo
          enddo
       enddo

       ! x-velocity  (x are edges, y and z are centers)
       do k = lo(3), hi(3)
          z = (dble(k)+HALF)*dx(3) + prob_lo(3)

          do j = lo(2), hi(2)
             y = (dble(j)+HALF)*dx(2) + prob_lo(2)

             do i = lo(1), hi(1)+1
                x = (dble(i))*dx(1) + prob_lo(1)

                gphi_x(i,j,k) = (phi(i,j,k) - phi(i-1,j,k))/dx(1)
                umac(i,j,k) = umac(i,j,k) + gphi_x(i,j,k)

             enddo
          enddo
       enddo

       ! y-velocity  (x and z are centers, y are edges)
       do k = lo(3), hi(3)
          z = (dble(k)+HALF)*dx(3) + prob_lo(3)

          do j = lo(2), hi(2)+1
             y = (dble(j))*dx(2) + prob_lo(2)

             do i = lo(1), hi(1)
                x = (dble(i)+HALF)*dx(1) + prob_lo(1)

                gphi_y(i,j,k) = (phi(i,j,k) - phi(i,j-1,k))/dx(2)
                vmac(i,j,k) = vmac(i,j,k) + gphi_y(i,j,k)

             enddo
          enddo
       enddo

       ! z-velocity  (x and y are centers, z are edges)
       do k = lo(3), hi(3)+1
          z = (dble(k))*dx(3) + prob_lo(3)

          do j = lo(2), hi(2)
             y = (dble(j)+HALF)*dx(2) + prob_lo(2)

             do i = lo(1), hi(1)
                x = (dble(i)+HALF)*dx(1) + prob_lo(1)

                gphi_z(i,j,k) = (phi(i,j,k) - phi(i,j,k-1))/dx(3)
                wmac(i,j,k) = wmac(i,j,k) + gphi_z(i,j,k)

             enddo
          enddo
       enddo

    else
       call bl_error('Not set up for these boundary conditions')
    end if


    ! impose BCs

    ! x lo
    select case (box_phys_bc(1,1))
    case (SlipWall)
       umac(lo(1),lo(2):hi(2),lo(3):hi(3)) = ZERO
    case (Interior, INTERIOR)
    case default
       call bl_error("invalid x lo BC")
    end select

    ! x hi
    select case(box_phys_bc(1,2))
    case (SlipWall)
       umac(hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = ZERO
    case (Interior, INTERIOR)
    case default
       call bl_error("invalid x hi BC")
    end select

    ! y lo
    select case (box_phys_bc(2,1))
    case (SlipWall)
       vmac(lo(1):hi(1),lo(2),lo(3):hi(3)) = ZERO
    case (Interior, INTERIOR)
    case default
       call bl_error("invalid y lo BC")
    end select

    ! y hi
    select case(box_phys_bc(2,2))
    case (SlipWall)
       vmac(lo(1):hi(1),hi(2)+1,lo(3):hi(3)) = ZERO
    case (Interior, INTERIOR)
    case default
       call bl_error("invalid y hi BC")
    end select

    ! z lo
    select case (box_phys_bc(3,1))
    case (SlipWall)
       wmac(lo(1):hi(1),lo(2):hi(2),lo(3)) = ZERO
    case (Interior, INTERIOR)
    case default
       call bl_error("invalid z lo BC")
    end select

    ! z hi
    select case(box_phys_bc(3,2))
    case (SlipWall)
       wmac(lo(1):hi(1),lo(2):hi(2),hi(3)+1) = ZERO
    case (Interior, INTERIOR)
    case default
       call bl_error("invalid z hi BC")
    end select

    deallocate(phi)

  end subroutine add_grad_scalar_3d_mac


  !===========================================================================
  subroutine convert_MAC_to_cc(lo, hi, umac, m_lo, m_hi, nc_m, u, u_lo, u_hi, nc_u) bind(C, name="convert_MAC_to_cc")

    ! convert a MAC velocity field to a cell-centered one -- no ghost
    ! cell filling is done here.

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: m_lo(3), m_hi(3)
    integer, intent(in) :: u_lo(3), u_hi(3)
    double precision, intent(in) :: umac(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3),1:nc_m)
    double precision, intent(in) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),1:nc_u)

    integer :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             u(i,j,k,1) = HALF*(umac(i,j,k) + umac(i+1,j,k))
             u(i,j,k,2) = HALF*(vmac(i,j,k) + vmac(i,j+1,k))
             u(i,j,k,3) = HALF*(wmac(i,j,k) + wmac(i,j,k+1))

          enddo
       enddo
    enddo


  end subroutine convert_MAC_to_cc


  subroutine get_project_type(ptype) bind(C, name="get_project_type")

    integer, intent(inout) :: ptype

    ptype = project_type

  end subroutine get_project_type

end module test_projection_module

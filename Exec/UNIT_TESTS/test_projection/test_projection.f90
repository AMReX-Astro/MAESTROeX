module test_projection_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use ml_cc_restriction_module
  use ml_layout_module
  use define_bc_module
  use multifab_fill_ghost_module
  use multifab_physbc_module

  implicit none

  private
  public :: init_velocity, init_mac_velocity, &
       add_grad_scalar, add_grad_scalar_mac, convert_MAC_to_cc

contains

  !===========================================================================
  subroutine init_velocity(U, dx, mla, the_bc_level)

    integer :: n, i, ng, dm, nlevs

    type(multifab) , intent(inout) :: U(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(ml_layout)   , intent(inout) :: mla
    type(bc_level)    , intent(in   ) :: the_bc_level(:)

    integer :: lo(get_dim(U(1))), hi(get_dim(U(1)))

    real(kind=dp_t), pointer :: up(:,:,:,:)

    nlevs = size(U)
    dm = get_dim(U(1))

    ng = nghost(U(1))

    do n=1,nlevs
       do i = 1, nfabs(U(n))
          up => dataptr(U(n), i)
          lo = lwb(get_box(U(n), i))
          hi = upb(get_box(U(n), i))

          select case (dm)
          case (2)
             call init_velocity_2d(up(:,:,1,:), ng, lo, hi, dx(n,:))

          case (3)
             call init_velocity_3d(up(:,:,:,:), ng, lo, hi, dx(n,:))

          end select
       end do
    end do

    ! fill ghostcells
    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(U(nlevs))

       ! fill non-periodic domain boundary ghost cells                         
       call multifab_physbc(U(nlevs),1,1,dm,the_bc_level(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the             
       ! finer grids are done first                                            
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n                
          ! data covering it                                                   
          call ml_cc_restriction(U(n-1),U(n),mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from                  
          ! level n-1 data note that multifab_fill_boundary and                
          ! multifab_physbc are called for both levels n-1 and n               
          call multifab_fill_ghost_cells(U(n),U(n-1),nghost(U(n)), &
                                         mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1), &
                                         the_bc_level(n), &
                                         1,1,dm)
          
       enddo

    end if

  end subroutine init_velocity

  subroutine init_velocity_2d(U, ng, lo, hi, dx)

    ! initialize the velocity field to a divergence-free field.  This
    ! velocity field comes from Almgren, Bell, and Szymczak 1996.

    use probin_module, only: prob_lo

    integer         , intent(in   ) :: lo(:), hi(:), ng
    real (kind=dp_t), intent(inout) :: U(lo(1)-ng:,lo(2)-ng:,:)
    real (kind=dp_t), intent(in   ) :: dx(:)

    ! Local variables
    integer :: i, j
    real (kind=dp_t) :: x, y

    do j = lo(2), hi(2)
       y = (dble(j)+0.5d0)*dx(2) + prob_lo(2)

       do i = lo(1), hi(1)
          x = (dble(i)+0.5d0)*dx(1) + prob_lo(1)
    
          U(i,j,1) = -sin(M_PI*x)**2 * sin(TWO*M_PI*y)
          U(i,j,2) =  sin(M_PI*y)**2 * sin(TWO*M_PI*x)  

       enddo
    enddo

  end subroutine init_velocity_2d

  subroutine init_velocity_3d(U, ng, lo, hi, dx)

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
    ! (this is periodic and even in all directions)
    !
    ! then U = curl{Phi} gives our field

    use probin_module, only: prob_lo

    integer         , intent(in   ) :: lo(:), hi(:), ng
    real (kind=dp_t), intent(inout) :: U(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(in   ) :: dx(:)

    ! Local variables
    integer :: i, j, k
    real (kind=dp_t) :: x, y, z

    do k = lo(3), hi(3)
       z = (dble(k)+0.5d0)*dx(3) + prob_lo(3)

       do j = lo(2), hi(2)
          y = (dble(j)+0.5d0)*dx(2) + prob_lo(2)

          do i = lo(1), hi(1)
             x = (dble(i)+0.5d0)*dx(1) + prob_lo(1)
    
             U(i,j,k,1) = TWO*M_PI*sin(FOUR*M_PI*x)*cos( TWO*M_PI*y) - &
                         FOUR*M_PI*sin( TWO*M_PI*x)*cos(FOUR*M_PI*z)

             U(i,j,k,2) = TWO*M_PI*sin(FOUR*M_PI*y)*cos( TWO*M_PI*z) - &
                         FOUR*M_PI*cos(FOUR*M_PI*x)*sin( TWO*M_PI*y)

             U(i,j,k,3) = TWO*M_PI*cos( TWO*M_PI*x)*sin(FOUR*M_PI*z) - &
                         FOUR*M_PI*cos(FOUR*M_PI*y)*sin( TWO*M_PI*z)

          enddo
       enddo
    enddo

  end subroutine init_velocity_3d


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
    ! (this is periodic and even in all directions)
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
  subroutine add_grad_scalar(U, gphi, dx, mla, the_bc_level)

    integer :: n, i, ng, dm, nlevs

    type(multifab) , intent(inout) :: U(:)
    type(multifab) , intent(inout) :: gphi(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(ml_layout)   , intent(inout) :: mla
    type(bc_level)    , intent(in   ) :: the_bc_level(:)

    integer :: lo(get_dim(U(1))), hi(get_dim(U(1)))

    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: gp(:,:,:,:)

    nlevs = size(U)
    dm = get_dim(U(1))

    ng = nghost(U(1))

    do n=1,nlevs
       do i = 1, nfabs(U(n))
          up => dataptr(U(n), i)
          gp => dataptr(gphi(n), i)
          lo = lwb(get_box(U(n), i))
          hi = upb(get_box(U(n), i))

          select case (dm)
          case (2)
             call add_grad_scalar_2d(up(:,:,1,:), gp(:,:,1,:), ng, lo, hi, dx(n,:), &
                                     the_bc_level(n)%phys_bc_level_array(0,:,:), &
                                     the_bc_level(n)%phys_bc_level_array(i,:,:))

          case (3)
             call add_grad_scalar_3d(up(:,:,:,:), gp(:,:,:,:), ng, lo, hi, dx(n,:), &
                                     the_bc_level(n)%phys_bc_level_array(0,:,:), &
                                     the_bc_level(n)%phys_bc_level_array(i,:,:))

          end select
       end do
    end do

    ! fill ghostcells
    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(U(nlevs))

       ! fill non-periodic domain boundary ghost cells                         
       call multifab_physbc(U(nlevs),1,1,dm,the_bc_level(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the             
       ! finer grids are done first                                            
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n                
          ! data covering it                                                   
          call ml_cc_restriction(U(n-1),U(n),mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from                  
          ! level n-1 data note that multifab_fill_boundary and                
          ! multifab_physbc are called for both levels n-1 and n               
          call multifab_fill_ghost_cells(U(n),U(n-1),nghost(U(n)), &
                                         mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1), &
                                         the_bc_level(n), &
                                         1,1,dm)
          
       enddo

    end if

  end subroutine add_grad_scalar

  subroutine add_grad_scalar_2d(U, gphi, ng, lo, hi, dx, &
                                domain_phys_bc, box_phys_bc)

    use     bc_module
    use probin_module, only: prob_lo

    integer         , intent(in   ) :: lo(:), hi(:), ng
    real (kind=dp_t), intent(inout) ::    U(lo(1)-ng:,lo(2)-ng:,:)
    real (kind=dp_t), intent(inout) :: gphi(lo(1)-ng:,lo(2)-ng:,:)
    real (kind=dp_t), intent(in   ) :: dx(:)
    integer         , intent(in   ) :: domain_phys_bc(:,:)
    integer         , intent(in   ) :: box_phys_bc(:,:)


    ! Local variables
    integer :: i, j
    real (kind=dp_t) :: x, y
    real (kind=dp_t), allocatable :: phi(:,:)

    allocate(phi(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1)) 

    if (domain_phys_bc(1,1) .eq. SLIP_WALL .and. &
        domain_phys_bc(1,2) .eq. SLIP_WALL .and. &
        domain_phys_bc(2,1) .eq. SLIP_WALL .and. &
        domain_phys_bc(2,2) .eq. SLIP_WALL) then

       ! Add on the gradient of a scalar (phi) that satisfies
       ! grad(phi).n = 0.
       do j = lo(2), hi(2)
          y = (dble(j)+0.5d0)*dx(2) + prob_lo(2)

          do i = lo(1), hi(1)
             x = (dble(i)+0.5d0)*dx(1) + prob_lo(1)
      
             gphi(i,j,1) =  FOUR*x*(ONE - x)
             gphi(i,j,2) =  FOUR*y*(ONE - y)

             U(i,j,1) = U(i,j,1) + gphi(i,j,1)
             U(i,j,2) = U(i,j,2) + gphi(i,j,2)
   
          enddo
       enddo

    else if (domain_phys_bc(1,1) .eq. PERIODIC .and. &
             domain_phys_bc(1,2) .eq. PERIODIC .and. &
             domain_phys_bc(2,1) .eq. PERIODIC .and. &
             domain_phys_bc(2,2) .eq. PERIODIC) then

       do j = lo(2)-1, hi(2)+1
          y = (dble(j)+0.5d0)*dx(2) + prob_lo(2)
    
          do i = lo(1)-1, hi(1)+1
             x = (dble(i)+0.5d0)*dx(1) + prob_lo(1)
   
               phi(i,j) = 0.1d0 * cos(2.d0*M_PI*y)*cos(2.d0*M_PI*x)
   
          enddo
       enddo

       do j = lo(2), hi(2)
          y = (dble(j)+0.5d0)*dx(2) + prob_lo(2)

          do i = lo(1), hi(1)
             x = (dble(i)+0.5d0)*dx(1) + prob_lo(1)

             gphi(i,j,1) =  (phi(i+1,j)-phi(i-1,j))/(2.d0*dx(1))
             gphi(i,j,2) =  (phi(i,j+1)-phi(i,j-1))/(2.d0*dx(2))
   
             U(i,j,1) = U(i,j,1) + gphi(i,j,1)
             U(i,j,2) = U(i,j,2) + gphi(i,j,2)
   
          enddo
       enddo

    else
       call bl_error('Not set up for these boundary conditions')
    end if

    deallocate(phi)

  end subroutine add_grad_scalar_2d


  subroutine add_grad_scalar_3d(U, gphi, ng, lo, hi, dx, &
                                domain_phys_bc, box_phys_bc)

    use     bc_module
    use probin_module, only: prob_lo

    integer         , intent(in   ) :: lo(:), hi(:), ng
    real (kind=dp_t), intent(inout) ::    U(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(inout) :: gphi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(in   ) :: dx(:)
    integer         , intent(in   ) :: domain_phys_bc(:,:)
    integer         , intent(in   ) :: box_phys_bc(:,:)


    ! Local variables
    integer :: i, j, k
    real (kind=dp_t) :: x, y, z
    real (kind=dp_t), allocatable :: phi(:,:,:)

    allocate(phi(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)) 

    if (domain_phys_bc(1,1) .eq. SLIP_WALL .and. &
        domain_phys_bc(1,2) .eq. SLIP_WALL .and. &
        domain_phys_bc(2,1) .eq. SLIP_WALL .and. &
        domain_phys_bc(2,2) .eq. SLIP_WALL .and. &
        domain_phys_bc(3,1) .eq. SLIP_WALL .and. &
        domain_phys_bc(3,2) .eq. SLIP_WALL) then

       ! Add on the gradient of a scalar (phi) that satisfies
       ! grad(phi).n = 0.
       do k = lo(3), hi(3)
          z = (dble(k)+0.5d0)*dx(3) + prob_lo(3)

          do j = lo(2), hi(2)
             y = (dble(j)+0.5d0)*dx(2) + prob_lo(2)

             do i = lo(1), hi(1)
                x = (dble(i)+0.5d0)*dx(1) + prob_lo(1)
      
                gphi(i,j,k,1) =  160.0_dp_t*x*(ONE - x)
                gphi(i,j,k,2) =  160.0_dp_t*y*(ONE - y)
                gphi(i,j,k,3) =  160.0_dp_t*z*(ONE - z)

                U(i,j,k,1) = U(i,j,k,1) + gphi(i,j,k,1)
                U(i,j,k,2) = U(i,j,k,2) + gphi(i,j,k,2)
                U(i,j,k,3) = U(i,j,k,3) + gphi(i,j,k,3)
   
             enddo
          enddo
       enddo

    else if (domain_phys_bc(1,1) .eq. PERIODIC .and. &
             domain_phys_bc(1,2) .eq. PERIODIC .and. &
             domain_phys_bc(2,1) .eq. PERIODIC .and. &
             domain_phys_bc(2,2) .eq. PERIODIC .and. &
             domain_phys_bc(3,1) .eq. PERIODIC .and. &
             domain_phys_bc(3,2) .eq. PERIODIC) then

       do k = lo(3)-1, hi(3)+1
          z = (dble(k)+0.5d0)*dx(3) + prob_lo(3)

          do j = lo(2)-1, hi(2)+1
             y = (dble(j)+0.5d0)*dx(2) + prob_lo(2)
    
             do i = lo(1)-1, hi(1)+1
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
    end if

    deallocate(phi)

  end subroutine add_grad_scalar_3d


  !===========================================================================
  subroutine add_grad_scalar_mac(umac, gphi_mac, dx, mla, the_bc_level)

    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(inout) :: gphi_mac(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(ml_layout)   , intent(inout) :: mla
    type(bc_level)    , intent(in   ) :: the_bc_level(:)

    integer :: n, i, ng_um, ng_gp, dm, nlevs

    integer :: lo(get_dim(umac(1,1))), hi(get_dim(umac(1,1)))

    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)

    real(kind=dp_t), pointer :: gxp(:,:,:,:)
    real(kind=dp_t), pointer :: gyp(:,:,:,:)
    real(kind=dp_t), pointer :: gzp(:,:,:,:)

    nlevs = size(umac(:,1))
    dm = get_dim(umac(1,1))

    ng_um = nghost(umac(1,1))
    ng_gp = nghost(gphi_mac(1,1))

    do n=1,nlevs
       do i = 1, nfabs(umac(n,1))
          ump => dataptr(umac(n,1), i)
          vmp => dataptr(umac(n,2), i)

          gxp => dataptr(gphi_mac(n,1), i)
          gyp => dataptr(gphi_mac(n,2), i)

          lo = lwb(get_box(umac(n,1), i))
          hi = upb(get_box(umac(n,1), i))

          select case (dm)
          case (2)
             call add_grad_scalar_2d_mac(ump(:,:,1,1), vmp(:,:,1,1), ng_um, &
                                         gxp(:,:,1,1), gyp(:,:,1,1), ng_gp, &
                                         lo, hi, dx(n,:), &
                                         the_bc_level(n)%phys_bc_level_array(0,:,:), &
                                         the_bc_level(n)%phys_bc_level_array(i,:,:))

          case (3)
             wmp => dataptr(umac(n,3), i)
             gzp => dataptr(gphi_mac(n,3), i)
             call add_grad_scalar_3d_mac(ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), ng_um, &
                                         gxp(:,:,:,1), gyp(:,:,:,1), gzp(:,:,:,1), ng_gp, &
                                         lo, hi, dx(n,:), &
                                         the_bc_level(n)%phys_bc_level_array(0,:,:), &
                                         the_bc_level(n)%phys_bc_level_array(i,:,:))

          end select
       end do
    end do

    ! make edge states consistent across levels
    do n = nlevs,2,-1
       do i = 1, dm
          call ml_edge_restriction_c(umac(n-1,i),1,umac(n,i),1,mla%mba%rr(n-1,:),i,1)
          call ml_edge_restriction_c(gphi_mac(n-1,i),1,gphi_mac(n,i),1,mla%mba%rr(n-1,:),i,1)
       enddo
    enddo

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

    if (domain_phys_bc(1,1) .eq. SLIP_WALL .and. &
        domain_phys_bc(1,2) .eq. SLIP_WALL .and. &
        domain_phys_bc(2,1) .eq. SLIP_WALL .and. &
        domain_phys_bc(2,2) .eq. SLIP_WALL) then

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


    else if (domain_phys_bc(1,1) .eq. PERIODIC .and. &
             domain_phys_bc(1,2) .eq. PERIODIC .and. &
             domain_phys_bc(2,1) .eq. PERIODIC .and. &
             domain_phys_bc(2,2) .eq. PERIODIC) then

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
    case (SLIP_WALL)
       umac(lo(1),lo(2):hi(2)) = ZERO
    case (PERIODIC, INTERIOR)
    case default
       call bl_error("invalid x lo BC")
    end select

    ! x hi
    select case(box_phys_bc(1,2))
    case (SLIP_WALL)
       umac(hi(1)+1,lo(2):hi(2)) = ZERO
    case (PERIODIC, INTERIOR)
    case default
       call bl_error("invalid x hi BC")
    end select

    ! y lo
    select case (box_phys_bc(2,1))
    case (SLIP_WALL)
       vmac(lo(1):hi(1),lo(2)) = ZERO
    case (PERIODIC, INTERIOR)
    case default
       call bl_error("invalid y lo BC")
    end select

    ! y hi
    select case(box_phys_bc(2,2))
    case (SLIP_WALL)
       vmac(lo(1):hi(1),hi(2)+1) = ZERO
    case (PERIODIC, INTERIOR)
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

    if (domain_phys_bc(1,1) .eq. SLIP_WALL .and. &
        domain_phys_bc(1,2) .eq. SLIP_WALL .and. &
        domain_phys_bc(2,1) .eq. SLIP_WALL .and. &
        domain_phys_bc(2,2) .eq. SLIP_WALL .and. &
        domain_phys_bc(3,1) .eq. SLIP_WALL .and. &
        domain_phys_bc(3,2) .eq. SLIP_WALL) then

       ! Add on the gradient of a scalar (phi) that satisfies
       ! grad(phi).n = 0.

       ! x-velocity  (x are edges, y and z are centers)
       do k = lo(3), hi(3)
          z = (dble(k)+HALF)*dx(3) + prob_lo(3)

          do j = lo(2), hi(2)
             y = (dble(j)+HALF)*dx(2) + prob_lo(2)

             do i = lo(1), hi(1)+1
                x = (dble(i))*dx(1) + prob_lo(1)
      
                gphi_x(i,j,k) =  160.0_dp_t*x*(ONE - x)
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
      
                gphi_y(i,j,k) =  160.0_dp_t*y*(ONE - y)
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
      
                gphi_z(i,j,k) =  160.0_dp_t*z*(ONE - z)
                wmac(i,j,k) = wmac(i,j,k) + gphi_z(i,j,k)
   
             enddo
          enddo
       enddo


    else if (domain_phys_bc(1,1) .eq. PERIODIC .and. &
             domain_phys_bc(1,2) .eq. PERIODIC .and. &
             domain_phys_bc(2,1) .eq. PERIODIC .and. &
             domain_phys_bc(2,2) .eq. PERIODIC .and. &
             domain_phys_bc(3,1) .eq. PERIODIC .and. &
             domain_phys_bc(3,2) .eq. PERIODIC) then


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
    case (SLIP_WALL)
       umac(lo(1),lo(2):hi(2),lo(3):hi(3)) = ZERO
    case (PERIODIC, INTERIOR)
    case default
       call bl_error("invalid x lo BC")
    end select

    ! x hi
    select case(box_phys_bc(1,2))
    case (SLIP_WALL)
       umac(hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = ZERO
    case (PERIODIC, INTERIOR)
    case default
       call bl_error("invalid x hi BC")
    end select

    ! y lo
    select case (box_phys_bc(2,1))
    case (SLIP_WALL)
       vmac(lo(1):hi(1),lo(2),lo(3):hi(3)) = ZERO
    case (PERIODIC, INTERIOR)
    case default
       call bl_error("invalid y lo BC")
    end select

    ! y hi
    select case(box_phys_bc(2,2))
    case (SLIP_WALL)
       vmac(lo(1):hi(1),hi(2)+1,lo(3):hi(3)) = ZERO
    case (PERIODIC, INTERIOR)
    case default
       call bl_error("invalid y hi BC")
    end select

    ! z lo
    select case (box_phys_bc(3,1))
    case (SLIP_WALL)
       wmac(lo(1):hi(1),lo(2):hi(2),lo(3)) = ZERO
    case (PERIODIC, INTERIOR)
    case default
       call bl_error("invalid z lo BC")
    end select

    ! z hi
    select case(box_phys_bc(3,2))
    case (SLIP_WALL)
       wmac(lo(1):hi(1),lo(2):hi(2),hi(3)+1) = ZERO
    case (PERIODIC, INTERIOR)
    case default
       call bl_error("invalid z hi BC")
    end select

    deallocate(phi)

  end subroutine add_grad_scalar_3d_mac


  !===========================================================================  
  subroutine convert_MAC_to_cc(umac, u)

    ! convert a MAC velocity field to a cell-centered one -- no ghost
    ! cell filling is done here.

    integer :: n, i, ng_um, ng_u, dm, nlevs

    type(multifab) , intent(in   ) :: umac(:,:)
    type(multifab) , intent(inout) :: u(:)

    integer :: lo(get_dim(u(1))), hi(get_dim(u(1)))

    real(kind=dp_t), pointer :: ump(:,:,:,:), vmp(:,:,:,:), wmp(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)

    nlevs = size(u)
    dm = get_dim(u(1))

    ng_u  = nghost(u(1))
    ng_um = nghost(umac(1,1))

    do n=1,nlevs
       do i = 1, nfabs(u(n))
          up => dataptr(u(n), i)

          ump => dataptr(umac(n,1), i)
          vmp => dataptr(umac(n,2), i)

          lo = lwb(get_box(u(n), i))
          hi = upb(get_box(u(n), i))

          select case (dm)
          case (2)
             call convert_MAC_to_cc_2d(ump(:,:,1,1), vmp(:,:,1,1), ng_um, &
                                       up(:,:,1,:), ng_u, lo, hi)

          case (3)
             wmp => dataptr(umac(n,3), i)
             call convert_MAC_to_cc_3d(ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), ng_um, &
                                       up(:,:,:,:), ng_u, lo, hi)

          end select
       end do
    end do
    
  end subroutine convert_MAC_to_cc

  subroutine convert_MAC_to_cc_2d(umac, vmac, ng_um, u, ng_u, lo, hi)

    integer         , intent(in   ) :: lo(:), hi(:), ng_um, ng_u
    real (kind=dp_t), intent(in   ) :: umac(lo(1)-ng_um:,lo(2)-ng_um:)
    real (kind=dp_t), intent(in   ) :: vmac(lo(1)-ng_um:,lo(2)-ng_um:)
    real (kind=dp_t), intent(inout) ::    u(lo(1)-ng_u :,lo(2)-ng_u :,:)

    integer :: i, j

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          u(i,j,1) = HALF*(umac(i,j) + umac(i+1,j))
          u(i,j,2) = HALF*(vmac(i,j) + vmac(i,j+1))

       enddo
    enddo
       
  end subroutine convert_MAC_to_cc_2d

  subroutine convert_MAC_to_cc_3d(umac, vmac, wmac, ng_um, u, ng_u, lo, hi)

    integer         , intent(in   ) :: lo(:), hi(:), ng_um, ng_u
    real (kind=dp_t), intent(in   ) :: umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real (kind=dp_t), intent(in   ) :: vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real (kind=dp_t), intent(in   ) :: wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real (kind=dp_t), intent(inout) ::    u(lo(1)-ng_u :,lo(2)-ng_u :,lo(3)-ng_u :,:)

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
       
  end subroutine convert_MAC_to_cc_3d

end module test_projection_module

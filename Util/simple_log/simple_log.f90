! the simple_log_module provides a simple buffer that stores strings
! and writes to the screen at the same time.
!
! at the moment, the log size is fixed to MAX_LINES -- in the future
! we can update it reallocate the buffer as needed.

module simple_log_module

  use bl_types
  use bl_error_module

  implicit none

  integer, parameter :: MAX_COLS = 128
  integer, parameter :: MAX_LINES = 128

  character (len=MAX_COLS), allocatable, save :: log_data(:)
  integer, save :: log_lines = 0

  logical, save :: initialized = .false.

  interface log
     module procedure sd_log
     module procedure si_log
     module procedure s_log
  end interface log

contains

  subroutine simple_log_init()
    allocate(log_data(MAX_LINES))
    initialized = .true.
  end subroutine simple_log_init

  subroutine simple_log_finalize()
    deallocate(log_data)
  end subroutine simple_log_finalize
  
  subroutine sd_log(str, d)
    character (len=*), intent(in) :: str    
    real (kind=dp_t), intent(in) :: d

    character (len=20) :: dstr

    write (dstr,'(g20.10)') d
    call s_log(str // dstr)
  end subroutine sd_log

  subroutine si_log(str, i)
    character (len=*), intent(in) :: str    
    integer, intent(in) :: i

    character (len=8) :: istr

    write (istr,'(i8)') i
    call s_log(str // istr)
  end subroutine si_log

  subroutine s_log(str)
    character (len=*), intent(in) :: str

    if (.not. initialized) call simple_log_init()
    
    if (log_lines > MAX_LINES-1) then
       call bl_error("ERROR: log buffer exceeded in module simple_log")
    endif

    ! output to the screen
    write (*,*) trim(str)

    ! and store in the log
    log_lines = log_lines + 1
    if (len(str) > MAX_COLS) then
       log_data(log_lines) = trim(str(1:MAX_COLS))
    else
       log_data(log_lines) = trim(str)
    endif
  end subroutine s_log

  subroutine log_break()
    if (log_lines > MAX_LINES-1) then
       call bl_error("ERROR: log buffer exceeded in module simple_log")
    endif

    log_lines = log_lines + 1
    log_data(log_lines) = "--------------------------------------------------------------------------------"

  end subroutine log_break

end module simple_log_module
  

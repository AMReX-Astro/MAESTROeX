module problem_functions_module

  use amrex_constants_module
  use probin_module, only: min_time_step, react_its, run_prefix
  ! use meth_params_module, only: prob_lo

  implicit none

  private

contains

  subroutine get_min_timestep(dt) bind(C, name="get_min_timestep")

    double precision, intent(inout) :: dt

    dt = min_time_step

  end subroutine get_min_timestep

  subroutine get_react_its(its) bind(C, name="get_react_its")

    integer, intent(inout) :: its

    its = react_its

  end subroutine get_react_its

  subroutine get_run_prefix(prefix, max_len) bind(C, name="get_run_prefix")

    integer, value, intent(in) :: max_len
    character, intent(inout) :: prefix(max_len)

    integer :: i, length

    length = len_trim(run_prefix)

    if (length > max_len-1) then
       call bl_error("length of run_prefix is greater than the length of the character array you've given me")
    endif

    do i = 1, length
       prefix(i) = run_prefix(i:i)
    enddo

    ! c++ char arrays are null-terminated
    prefix(length+1) = char(0)

  end subroutine get_run_prefix

end module problem_functions_module



subroutine ca_network_init() bind(C, name="ca_network_init")

  use network, only: network_init

  call network_init()

end subroutine ca_network_init


subroutine ca_get_num_spec(nspec_out) &
     bind(C, name="ca_get_num_spec")

  use network, only: nspec
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer, intent(out) :: nspec_out

  nspec_out = nspec

end subroutine ca_get_num_spec


subroutine ca_get_spec_names(spec_names,ispec,len) &
     bind(C, name="ca_get_spec_names")

  use network, only: nspec, short_spec_names
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer, intent(in   ) :: ispec
  integer, intent(inout) :: len
  integer, intent(inout) :: spec_names(len)

  ! Local variables
  integer   :: i

  len = len_trim(short_spec_names(ispec+1))

  do i = 1,len
     spec_names(i) = ichar(short_spec_names(ispec+1)(i:i))
  end do

end subroutine ca_get_spec_names

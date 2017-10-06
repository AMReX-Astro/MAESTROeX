! :::
! ::: ----------------------------------------------------------------
! :::

subroutine ca_get_num_spec(nspec_out) bind(C, name="ca_get_num_spec")

  use network, only: nspec
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer, intent(out) :: nspec_out

  nspec_out = nspec

end subroutine ca_get_num_spec

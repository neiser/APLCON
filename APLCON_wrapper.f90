module APLCON_wrapper
  use iso_c_binding
  implicit none
contains

  ! helper function to convert c string to fortran string
  function c_to_f_string(s) result(str)
    character(kind=c_char,len=1), intent(in) :: s(*)
    character(len=:), allocatable :: str
    integer i, nchars
    i = 1
    do
       if (s(i) == c_null_char) exit
       i = i + 1
    end do
    nchars = i - 1  ! Exclude null character from Fortran string
    allocate(character(len=nchars) :: str)
    str = transfer(s(1:nchars), str)
  end function c_to_f_string

  subroutine C_APLCON_APLCON(NVAR,MCST) bind(c)
    integer(c_int), value :: NVAR, MCST
    CALL APLCON(NVAR,MCST)
  end subroutine C_APLCON_APLCON

  subroutine C_APLCON_APNAME(I,NAME) bind(c)
    integer(c_int), value :: I
    character(kind=c_char,len=1), intent(in) :: NAME(*)
    CALL APNAME(I,NAME)
  end subroutine C_APLCON_APNAME

end module APLCON_wrapper

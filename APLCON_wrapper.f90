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
    integer(c_int), value, intent(in) :: NVAR, MCST
    CALL APLCON(NVAR,MCST)
  end subroutine C_APLCON_APLCON

  subroutine C_APLCON_APRINT(LUNP,IPR) bind(c)
    integer(c_int), value, intent(in) :: LUNP, IPR
    CALL APRINT(LUNP,IPR)
  end subroutine C_APLCON_APRINT
  
  subroutine C_APLCON_APNAME(I,NAME) bind(c)
    integer(c_int), value, intent(in) :: I
    character(kind=c_char,len=1), intent(in) :: NAME(*)
    CALL APNAME(I,c_to_f_string(NAME))
  end subroutine C_APLCON_APNAME

  subroutine C_APLCON_CIPRV(LUP,X,VX,N) bind(c)
    integer(c_int), value, intent(in) :: LUP, N
    real(c_double), dimension(*), intent(in) :: X, VX    
    CALL CIPRV(LUP,X,VX,N)
  end subroutine C_APLCON_CIPRV

  subroutine C_APLCON_CFCORR(LUP,V,N) bind(c)
    integer(c_int), value, intent(in) :: LUP, N
    real(c_double), dimension(*), intent(in) :: V
    CALL CFCORR(LUP,V,N)
  end subroutine C_APLCON_CFCORR

  subroutine C_APLCON_APROFL(I1,I2) bind(c)
    integer(c_int), value, intent(in) :: I1, I2
    CALL APROFL(I1,I2)
  end subroutine C_APLCON_APROFL

  subroutine C_APLCON_APLOOP(X,VX,F,IRET) bind(c)
    real(c_double), dimension(*), intent(inout) :: X,VX,F
    integer(c_int), intent(out) :: IRET
    CALL APLOOP(X,VX,F,IRET)
  end subroutine C_APLCON_APLOOP

  
end module APLCON_wrapper

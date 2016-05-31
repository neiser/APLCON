module APLCON_wrapper
  use iso_c_binding
  implicit none
contains

  ! main routines
  subroutine C_APLCON_APLCON(NVAR,MCST) bind(c)
    integer(c_int), value, intent(in) :: NVAR, MCST
    CALL APLCON(NVAR,MCST)
    CALL FLUSH
  end subroutine C_APLCON_APLCON

  subroutine C_APLCON_APLOOP(X,VX,F,IRET) bind(c)
    real(c_double), dimension(*), intent(inout) :: X,VX,F
    integer(c_int), intent(out) :: IRET
    CALL APLOOP(X,VX,F,IRET)
    CALL FLUSH
  end subroutine C_APLCON_APLOOP

  ! routines to obtain results
  subroutine C_APLCON_CHNDPV(CHI2,ND,PVAL) bind(c)
    real(c_float), intent(out) :: CHI2,PVAL
    integer(c_int), intent(out) :: ND
    CALL CHNDPV(CHI2,ND,PVAL)
  end subroutine C_APLCON_CHNDPV

  subroutine C_APLCON_APSTAT(FOPT,NFUN,NITER) bind(c)
    real(c_double), intent(out) :: FOPT
    integer(c_int), intent(out) :: NFUN, NITER
    CALL APSTAT(FOPT,NFUN,NITER)
  end subroutine C_APLCON_APSTAT

  subroutine C_APLCON_APPULL(PULLS) bind(c)
    real(c_double), dimension(*), intent(out) :: PULLS
    CALL APPULL(PULLS)
  end subroutine C_APLCON_APPULL

  ! variable reduction
  subroutine C_APLCON_SIMSEL(X,VX,NY,LIST,Y,VY) bind(c)
    real(c_double), dimension(*), intent(in) :: X,VX,LIST
    integer(c_int), value, intent(in) :: NY
    real(c_double), dimension(*), intent(out) :: Y,VY
    CALL SIMSEL(X,VX,NY,LIST,Y,VY)
  end subroutine C_APLCON_SIMSEL

  subroutine C_APLCON_SIMTRN(X,VX,NX) bind(c)
    real(c_double), dimension(*), intent(inout) :: X,VX
    integer(c_int), value, intent(in) :: NX
    CALL SIMTRN(X,VX,NX)
  end subroutine C_APLCON_SIMTRN

  ! printout related routines (useful for debugging)
  subroutine C_APLCON_CIPRV(LUP,X,VX,N) bind(c)
    integer(c_int), value, intent(in) :: LUP, N
    real(c_double), dimension(*), intent(in) :: X, VX
    CALL CIPRV(LUP,X,VX,N)
  end subroutine C_APLCON_CIPRV

  subroutine C_APLCON_CFPRV(LUP,X,VX,N) bind(c)
    integer(c_int), value, intent(in) :: LUP, N
    real(c_double), dimension(*), intent(in) :: X, VX
    CALL CFPRV(LUP,X,VX,N)
  end subroutine C_APLCON_CFPRV

  subroutine C_APLCON_CFCORR(LUP,V,N) bind(c)
    integer(c_int), value, intent(in) :: LUP, N
    real(c_double), dimension(*), intent(in) :: V
    CALL CFCORR(LUP,V,N)
  end subroutine C_APLCON_CFCORR

  ! ENTRY points in aplusopt.F
  ! mainly setup/config routines
  ! not all are mentioned in the README
  subroutine C_APLCON_APRINT(LUNP,IPR) bind(c)
    integer(c_int), value, intent(in) :: LUNP, IPR
    CALL APRINT(LUNP,IPR)
  end subroutine C_APLCON_APRINT

  subroutine C_APLCON_APDEPS(EPSF) bind(c)
    real(c_double), value, intent(in) :: EPSF
    CALL APDEPS(EPSF)
  end subroutine C_APLCON_APDEPS

  subroutine C_APLCON_APEPSCHI(EPSCHI) bind(c)
    real(c_double), value, intent(in) :: EPSCHI
    CALL APEPSCHI(EPSCHI)
  end subroutine C_APLCON_APEPSCHI

  subroutine C_APLCON_APDERF(DERFAC) bind(c)
    real(c_double), value, intent(in) :: DERFAC
    CALL APDERF(DERFAC)
  end subroutine C_APLCON_APDERF

  subroutine C_APLCON_APDERU(DERUFC) bind(c)
    real(c_double), value, intent(in) :: DERUFC
    CALL APDERU(DERUFC)
  end subroutine C_APLCON_APDERU

  subroutine C_APLCON_APDLOW(DERLOW) bind(c)
    real(c_double), value, intent(in) :: DERLOW
    CALL APDLOW(DERLOW)
  end subroutine C_APLCON_APDLOW

  subroutine C_APLCON_APITER(ITERMX) bind(c)
    integer(c_int), value, intent(in) :: ITERMX
    CALL APITER(ITERMX)
  end subroutine C_APLCON_APITER

  subroutine C_APLCON_APROFL(I1,I2) bind(c)
    integer(c_int), value, intent(in) :: I1, I2
    CALL APROFL(I1,I2)
  end subroutine C_APLCON_APROFL

  subroutine C_APLCON_APSTEP(I,STEP) bind(c)
    integer(c_int), value, intent(in) :: I
    real(c_double), value, intent(in) :: STEP
    CALL APSTEP(I,STEP)
  end subroutine C_APLCON_APSTEP

  subroutine C_APLCON_APFIX(I) bind(c)
    integer(c_int), value, intent(in) :: I
    CALL APFIX(I)
  end subroutine C_APLCON_APFIX

  subroutine C_APLCON_APLIMT(I,XLOW,XHIG) bind(c)
    integer(c_int), value, intent(in) :: I
    real(c_double), value, intent(in) :: XLOW, XHIG
    CALL APLIMT(I,XLOW,XHIG)
  end subroutine C_APLCON_APLIMT

  subroutine C_APLCON_APTRIN(I) bind(c)
    integer(c_int), value, intent(in) :: I
    CALL APTRIN(I)
  end subroutine C_APLCON_APTRIN

  subroutine C_APLCON_APOISS(I) bind(c)
    integer(c_int), value, intent(in) :: I
    CALL APOISS(I)
  end subroutine C_APLCON_APOISS

  subroutine C_APLCON_ABINOM(I) bind(c)
    integer(c_int), value, intent(in) :: I
    CALL ABINOM(I)
  end subroutine C_APLCON_ABINOM

  subroutine C_APLCON_APLOGN(I) bind(c)
    integer(c_int), value, intent(in) :: I
    CALL APLOGN(I)
  end subroutine C_APLCON_APLOGN

  subroutine C_APLCON_APSQRT(I) bind(c)
    integer(c_int), value, intent(in) :: I
    CALL APSQRT(I)
  end subroutine C_APLCON_APSQRT

  subroutine C_APLCON_APOWER(I) bind(c)
    integer(c_int), value, intent(in) :: I
    CALL APOWER(I)
  end subroutine C_APLCON_APOWER

  subroutine C_APLCON_APOSIT(I) bind(c)
    integer(c_int), value, intent(in) :: I
    CALL APOSIT(I)
  end subroutine C_APLCON_APOSIT

end module APLCON_wrapper

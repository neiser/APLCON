*     DCSPLN     calculate spline coefficients
*     DCSPR      print spline function 
*     DCSPF      spline function evaluation 
*     DCSPD      spline derivative evaluation
*     DCSPX      next argument with given function/derivative value
*     IDCSP      index for X argument

*     ******************************************************************

      SUBROUTINE DCSPLN(X,Y,C,N)
*
*     calculation of spline coefficients from pairs x,y
*
*                 - -   -
*     CALL DCSPLN(X,Y,C,N)
*                     -
*
*     end conditions are: not-a-knot
*
*     ARRAY C(1...5,N) after return:
*
*     C(1,I) =  F(X(I)) = Y(I)
*     C(2,I) =  first derivative  F'(X(I))
*     C(3,I) =  second derivative F''(X(I))/2.0
*     C(4,I) =  third derivative  F'''(X(I))/6.0
*     C(5,I) =  X(I)
*
*     calculation of interpolated value for given x :
*        with X(I) =< X =< X(I+1):
*     H=X-C(5,I)
*     F=C(1,I)+H*(C(2,I)+H*(C(3,I)+H*C(4,I)/3.0)/2.0)
*
*     ------------------------------------------------------------------
      IMPLICIT NONE

      DOUBLE PRECISION C(5,*),X(*),Y(*)
      DOUBLE PRECISION G,DX,DVD1,DVD3
      INTEGER N,I 
*     ...
      IF(N.LT.3) STOP 'N < 3 in DCSPLIN'
      DO I=1,N
       C(1,I)=Y(I)
       C(5,I)=X(I)
      END DO

      ENTRY ECSPLN(C,N)

      DO I=2,N
       C(3,I)= C(5,I)-C(5,I-1)          ! first differences
       C(4,I)=(C(1,I)-C(1,I-1))/C(3,I)  ! first divided differences
      END DO

      C(4,1)=C(3,3)         ! not-a-knot at left side
      C(3,1)=C(3,2)+C(3,3)
      C(2,1)=((C(3,2)+2.0D0*C(3,1))*C(4,2)*C(3,3)+C(3,2)**2*C(4,3))
     +        /C(3,1)

      DO I=2,N-1            ! generate equations and carry out 
       G=-C(3,I+1)/C(4,I-1) ! forward pass of Gauss elimination
       C(2,I)=G*C(2,I-1)+3.0D0*(C(3,I)*C(4,I+1)+C(3,I+1)*C(4,I))
       C(4,I)=G*C(3,I-1)+2.0D0*(C(3,I)+C(3,I+1))
      END DO

      IF(N.EQ.3) THEN       
         C(2,N)=2.0D0*C(4,N)  ! not-a-knot at right side (N=3)
         C(4,N)=1.0D0
         G=-1.0D0/C(4,N-1)
      ELSE
         G=C(3,N-1)+C(3,N)  ! not-a-knot at right side
         C(2,N)=((C(3,N)+2.0D0*G)*C(4,N)*C(3,N-1)
     +           +C(3,N)**2*(C(1,N-1)-C(1,N-2))/C(3,N-1))/G
         G=-G/C(4,N-1)
         C(4,N)=C(3,N-1)
      END IF
      C(4,N)=C(4,N)+G*C(3,N-1)  ! complete forward pass       
      C(2,N)=(G*C(2,N-1)+C(2,N))/C(4,N)

      DO I=N-1,1,-1             ! backward substitution
       C(2,I)=(C(2,I)-C(3,I)*C(2,I+1))/C(4,I)
      END DO

      DO I=2,N                  ! generate coefficients
       DX=C(3,I)
       DVD1=(C(1,I)-C(1,I-1))/DX
       DVD3=C(2,I-1)+C(2,I)-2.0D0*DVD1
       C(3,I-1)=(DVD1-C(2,I-1)-DVD3)/DX
       C(4,I-1)=(DVD3/DX)/DX
      END DO

      DX=C(5,N)-C(5,N-1)        ! coefficients at x(n)
      C(2,N)=C(2,N-1)+DX*(2.0D0*C(3,N-1)+DX*3.0D0*C(4,N-1))
      C(3,N)=C(3,N-1)+DX*3.0D0*C(4,N-1)
      C(4,N)=C(4,N-1)
      END

      SUBROUTINE DCSPR(C,N)
*
*     Print  table  of  abscissas,  function  values,  derivatives   and
*     integrals for 2*N-1 X-values
*                - -
*     CALL DCSPR(C,N)
*
*
      IMPLICIT NONE
      INTEGER N,I,NDER 
      DOUBLE PRECISION C(5,N),PR(5),X,DCSPF,DCSPD
*     ...
      WRITE(*,101) N
      DO I=1,2*N-1
       IF(MOD(I,2).EQ.1) X=C(5,(I+1)/2)
       IF(MOD(I,2).EQ.0) X=0.5D0*(C(5,I/2)+C(5,I/2+1))
       PR(1)=X
       DO NDER=1,4
        PR(NDER+1)=DCSPD(X,C,N,NDER-1)
       END DO
       PR(2)=DCSPF(X,C,N)
c       PR(6)=DCSPD(X,C,N,-1)
       IF(MOD(I,2).EQ.1) WRITE(*,102) PR
      END DO
      RETURN
  101 FORMAT(5X,'Printout of spline function with dimension N =',
     1   I4,' (CSP)'/
     2 5X,'X',14X,'f(X)',11X,5Hf'(X),10X,6Hf''(X),9X,7Hf'''(X))

  102 FORMAT(3G15.5,G14.4,G13.3)
      END

      SUBROUTINE DCSFPL(C,N) ! plot
      IMPLICIT NONE
      INTEGER N,I,NDER 
      DOUBLE PRECISION C(5,N),X,DCSPF !,DCSPD
      REAL XPL(256),YPL(256)
*     ...
C      CALL FPS ! reset
      NDER=1
      DO I=1,2*N-1
       IF(MOD(I,2).EQ.1) X=C(5,(I+1)/2)
       IF(MOD(I,2).EQ.0) X=0.5D0*(C(5,I/2)+C(5,I/2+1))
       XPL(I)=X 
       YPL(I)=DCSPF(X,C,N)
      END DO
      CALL FPNXY('*',2*N-1,XPL,YPL)
      CALL FPL 
      END


      SUBROUTINE DCSPG(C,N,NDER,XR,YR,M)
*
      IMPLICIT NONE
      INTEGER N,M,I,NDER
      DOUBLE PRECISION C(5,N),X,DCSPD,XR(*),YR(*)
*     ...
      M=2*N-1 
      DO I=1,M
       IF(MOD(I,2).EQ.1) X=C(5,(I+1)/2)
       IF(MOD(I,2).EQ.0) X=0.5D0*(C(5,I/2)+C(5,I/2+1))
       XR(I)=X
       YR(I)=DCSPD(X,C,N,NDER)
      END DO
      END


      DOUBLE PRECISION FUNCTION DCSPF(X,C,N)
*
*     function value for given argument x
*                - - -
*     F  = DCSPF(X,C,N) = function value of spline function at x
*
      IMPLICIT NONE
      INTEGER N,IDCSP,I
      DOUBLE PRECISION C(5,N),H,X
*     ...
      I=IDCSP(X,C,N)
      H=X-C(5,I)
      DCSPF=C(1,I)+H*(C(2,I)+H*(C(3,I)+H*C(4,I)))
      END

      DOUBLE PRECISION FUNCTION DCSPD(X,C,N,NDER)
*
*     calculate function value, derivative or integral value  for  given
*     argument X
*
*                - - - ----
*     FD = DCSPD(X,C,N,NDER)
*
*     FD = function value, derivative or integrated spline function
*     at X
*
*     NDER = -1   integration
*          =  0   function value
*          =  1   first derivative
*          =  2   second derivative
*          =  3   third derivative
*
*     For x outside the range of the spline function,  extrapolation  is
*     done for the  function  value  and  derivatives.  for  integration
*     however, the function is assumed to be zero outside the range. the
*     lower limit of integration is equal to the left edge  (C(5,1))  of
*     the spline function.
*
      IMPLICIT NONE
      INTEGER N,I,J,IDCSP,NDER
      DOUBLE PRECISION C(5,N)
      DOUBLE PRECISION F,H,X
*     ...
      DCSPD=0.0D0
      IF(NDER.GT.3.OR.NDER.LT.(-1)) RETURN

      I=IDCSP(X,C,N)  ! index
      H=X-C(5,I)
      IF(NDER.EQ.3) THEN       ! third derivative
         DCSPD=6.0D0*C(4,I)
      ELSE IF(NDER.EQ.2) THEN  ! second derivative
         DCSPD=2.0D0*C(3,I)+6.0D0*H*C(4,I)
      ELSE IF(NDER.EQ.1) THEN  ! first derivative
         DCSPD=C(2,I)+H*(2.0D0*C(3,I)+3.0D0*H*C(4,I))
      ELSE IF(NDER.EQ.0) THEN  ! function value
         DCSPD=C(1,I)+H*(C(2,I)+H*(C(3,I)+H*C(4,I)))
      ELSE                     !  integration
         IF(X.LT.C(5,1)) RETURN
         DO J=1,I
          H=C(5,J+1)-C(5,J)
          IF(J.EQ.I.AND.X.LT.C(5,N)) H=X-C(5,J)
          F=C(1,J)+H*(C(2,J)/2.0D0+H*(C(3,J)/3.0D0+H*C(4,J)/4.0D0))
          DCSPD=DCSPD+F*H
         END DO
      END IF
      END


      DOUBLE PRECISION FUNCTION DCSPN(XL,FD,C,N,NDER,IDIR)
*
*     calculate next x-value with given function value or derivative

*                -- -- - - ---- ----
*     X  = DCSPN(XL,FD,C,N,NDER,IDIR)      NDER = -1, 0, 1, 2
*          next abscissa value with dcspd(x,c,n,nder) = fd
 
      IMPLICIT NONE
      INTEGER N,NDER,IDIR,I,IT,IDCSP
      DOUBLE PRECISION C(5,N),XL,FD,DCSPD,XN,XA,XB,FA,FB,DB
      DOUBLE PRECISION FN,FP,DX
*     ...
      DCSPN=XL                  ! default return = start value
      I=IDCSP(XL,C,N)           ! initial value
      DX=0.01D0*(C(5,I+1)-C(5,I))
c      WRITE(*,*) 'Initial I N ',I,N,C(5,I),C(5,I+1)         
      IF(IDIR.LT.0) THEN
         XA=C(5,I)-DX
         FA=DCSPD(XA,C,N,NDER)-FD
         XB=XL+DX
         FB=DCSPD(XB,C,N,NDER)-FD
c         WRITE(*,*) '--- include ',XA,XB,FA,FB 
      ELSE
         XA=XL-DX
         FA=DCSPD(XA,C,N,NDER)-FD
         XB=C(5,I+1)+DX
         FB=DCSPD(XB,C,N,NDER)-FD
c         WRITE(*,*) '+++ include ',XA,XB,FA,FB
      END IF

      IF(FA*FB.GT.0.0D0) THEN
 10      I=I+SIGN(1,IDIR)       ! +1 or -1 
         IF(I.LT.1.OR.I.GE.N) RETURN
         XA=C(5,I)
         FA=DCSPD(XA,C,N,NDER)-FD
         XB=C(5,I+1)
         FB=DCSPD(XB,C,N,NDER)-FD
         IF(FA*FB.GT.0.0D0) GOTO 10
      END IF

      XN=0.5D0*(XA+XB)          ! center
      DB=C(5,I+1)-C(5,I)        ! interval length
      DO IT=1,10                ! reduce interval by 0.001
       FN=DCSPD(XN,C,N,NDER)-FD
       IF(FN*FA.GT.0.0D0) THEN
          FA=FN
          XA=XN
       ELSE 
          FB=FN
          XB=XN
       END IF
       XN=0.5D0*(XA+XB)
      END DO

      IT=0                      ! improve value by Newton steps
 20   IT=IT+1
      FN=DCSPD(XN,C,N,NDER)-FD ! function value
      FP=DCSPD(XN,C,N,NDER+1)  ! derivative
      XN=XN-FN/FP
      IF(IT.LT.10.AND.ABS(FN).GT.ABS(1.0D-12*DB*FP)) GOTO 20
      DCSPN=XN ! convergence
      END 

      INTEGER FUNCTION IDCSP(X,C,N) ! return index for X
*             C(5,I) <= X < C(5,I+1), if inside
      IMPLICIT NONE
      INTEGER N,I,IL,IH
      DOUBLE PRECISION X,C(5,N)
      DATA  I/1/
      SAVE I                    ! use evtl. previous index
*     ...
      IF(I.GE.N) I=(N+1)/2      ! correct start index
      IF(X.LT.C(5,I)) THEN
         IL=MAX(I-2,1)
         IF(X.LT.C(5,IL)) IL=1 
         IH=I
      ELSE 
         IL=I
         IH=MIN(I+2,N-1)
         IF(X.GE.C(5,IH)) IH=N-1
      END IF
 10   I=(IL+IH)/2              ! binary search
      IF(I.NE.IL) THEN
         IF(X.GE.C(5,I)) THEN    
            IL=I
         ELSE
            IH=I
         END IF
         GOTO 10
       END IF
       IDCSP=I
       END


 
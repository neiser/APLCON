

      SUBROUTINE A1PROF(KRET)        ! 1-parameter profile analysis
c      SUBROUTINE A1PROF(X,VX,KRET)        ! 1-parameter profile analysis
      IMPLICIT NONE
      INTEGER KRET !,J
#include "comcfit.inc"
#include "cprofil.inc"
#include "declarefl.inc"
#include "nauxfit.inc" 
      DOUBLE PRECISION XCD(100),YCD(100)     
*     KRET=-1 start next fit
*     KRET=
*     KRET= 
C     DOUBLE PRECISION X(*),VX(*)
*     SIGMAS = nr of parabolic sigmas (chi^2 limit is SIGMAS**2)
*     NKNOTS = nr of knot points in one direction
      DOUBLE PRECISION CHSUP 
      INTEGER    NKNOTS,KNOTS,ILRM,ILRP,NTLIMP,ILR1,ILR2
      REAL       SIGMAS
      PARAMETER (NKNOTS=6)
      PARAMETER (SIGMAS=4.4)
*     ...
c      WRITE(*,*) 'Entered A1PROF'
*
*     __________________________________________________________________
*     prepare 1-parameter points
*
      NADFS=1                  ! 1-parameter profile analysis
      IPF=NPSEC(1,ISECA)       ! profile index
      CENTER=AUX(INDAS+IPF)    ! fitted value of variable
c      WRITE(*,*) 'a12prof INDVS ',INDVS
c      WRITE(*,*) (AUX(INDVS+J),J=1,6)
      SIGMAX=SQRT(AUX(INDVS+(IPF*IPF+IPF)/2)) ! 

      XLR(0)=CENTER            ! original point
      FLR(0)=0.0D0             ! function difference
      IPAK=IPF                 ! index for unpacking            
#include "unpackfl.inc" 
      NTLIMP=NTLIM             ! property of variable IPF/IPAK  
      KNOTS=(2*NKNOTS+1)/2     ! make KNOTS even number
      CHSUP=SIGMAS**2          ! chi^2 limit
c      WRITE(*,*) 'KNOTS CHSUP',KNOTS,CHSUP
      REDL=(1.1*SIGMAS)/FLOAT(KNOTS) ! one step in units of sigma
      ILRM=MLR 
      REDR=REDL
      IF(NTLIMP.EQ.1) THEN     ! variable > 0 (positive)
         IF(CENTER.LT.SIGMAX*SIGMAS) REDL=CENTER/SIGMAS/FLOAT(KNOTS)
         ILRM=CENTER/SIGMAS/REDL
         ILRM=MIN(MLR,ILRM)    ! max number of steps to zero   
      END IF 
      DO I=1,MLR ! prepare all points
       CLR(+I)=+REDR*DFLOAT(I) ! unit of parab standard deviations
       CLR(-I)=-REDL*DFLOAT(I) ! unit of parab standard deviations
       XLR(+I)=XLR(0)+CLR(+I)*SIGMAX ! value of variable
       XLR(-I)=XLR(0)+CLR(-I)*SIGMAX ! value of variable
      END DO
c      WRITE(*,123) 'CLR',CLR
c      WRITE(*,123) 'XLR',XLR
c123  FORMAT(A/(4F15.4)) 
      ISTAR=0
      ILR = 0                  ! reset index of profile point
      ILRP=-2                  ! direction and step size
*
*     __________________________________________________________________
*     next profile point     
*
 10   IF(ILR*ILRP.GT.3.AND.
     +  ((FLR(ILR).GT.CHSUP.OR.ABS(ILR).GE.ILRM).OR.
     +   (FLR(ILR).LT.FLR(ILR-ILRP)))) THEN
         IF(ILRP.LT.0) THEN
            ILR1=ILR             ! end-point negative side
            ILRM=MLR             ! max number of steps
         END IF 
         IF(ILRP.GT.0) ILR2=ILR  ! end-point positive side
         ILR=ILR+ILRP/2          ! even -> odd
         ILRP=-ILRP              ! reverse direction
      END IF  
      ILR=ILR+ILRP               ! next value: change index by step
      CONSTR=XLR(ILR)            ! fixed parameter value
c      WRITE(*,*) '... using CONSTR=',CONSTR,ILR
c321  FORMAT(A/(4F15.5))
      NF=NFPRIM+1                ! number of constraints
      CALL APRINT(6,0)           ! suppress printout during profile fit       
*
*     __________________________________________________________________
*     return to fit ....
c      WRITE(*,*) 'return to fit'
c      WRITE(*,*) 'return to fit',INDST,(AUX(INDST+J),J=1,NX)
      KRET=-1  
*
*      ...
      RETURN
*
*     __________________________________________________________________
*     back from fit: store chi^2
*
C     ENTRY B1PROF(X,VX,KRET)    ! entry after fit with add. constraint
      ENTRY B1PROF(KRET)    ! entry after fit with add. constraint
c      WRITE(*,*) 'in B1PROF',INDST,(AUX(INDST+J),J=1,NX)
      FLR(ILR)=CHISQ-CHIREF      !  20 store chi^2 diff. of secondary fit
      IF(ILR+ILRP.NE.0) GOTO 10  ! continue 
*
*     __________________________________________________________________
*     search finished: spline, store, print 
*
      DO I=-1,ILR1,-1
       IF(FLR(I).LE.FLR(I+1)) GOTO 21
      END DO
      I=ILR1-1
 21   ILR1=I+1
      DO I=1,ILR2
       IF(FLR(I).LE.FLR(I-1)) GOTO 22
      END DO
      I=ILR2+1
 22   ILR2=I-1
      NLR=0                      ! reset number of points
c      WRITE(*,*) 'ILR1,ILR2',ILR1,ILR2
      DO I=ILR1,ILR2
c       WRITE(*,*) 'I xlr flr',I,XLR(I),FLR(I)
       IF(FLR(I).GE.0.0) THEN              ! fit converged
          NLR=NLR+1              ! nr of points
          XCD(NLR)=XLR(I)         ! value of variable
          YCD(NLR)=FLR(I)         ! function difference
          XC(NLR)=XLR(I)
          YC(NLR)=FLR(I)  
       END IF
      END DO
      CALL APROFJ(XCD,YCD,NLR)   ! internal
c      WRITE(*,*) 'in B1PROF vend',INDST,(AUX(INDST+J),J=1,NX) 
      CALL DCSPLN(XCD,YCD,CSP,NLR) ! spline definition
      CALL ASPROF(IPF,0,0,0.0,NLR,XC,YC)       ! store 
      CALL PPROFL                ! print profile info
      WRITE(*,*) 'PPROF1 '
      CALL PPROF1
      KRET=0                     ! this profile analysis finished
      END 


      SUBROUTINE A2PROF(KRET)      ! 2-parameter profile analysis
C      SUBROUTINE A2PROF(X,VX,KRET)      ! 2-parameter profile analysis
      IMPLICIT NONE
#include "comcfit.inc"
#include "cprofil.inc"
#include "declarefl.inc"
#include "nauxfit.inc"
      DOUBLE PRECISION XCD(100),YCD(100)
C      DOUBLE PRECISION X(*),VX(*)
      DOUBLE PRECISION CHSUP,FAC,DCSPN
      INTEGER KRET
      INTEGER KNOTS,NKNOTS,IJSYM,ILR1,ILR2,NTLIMX,NTLIMY,ILRM,ILRP,J
      DOUBLE PRECISION PI
      PARAMETER (PI=3.14159265358979323846D0) 
      REAL       SIGMAS
      PARAMETER (NKNOTS=6)
      PARAMETER (SIGMAS=4.4)
*
*     __________________________________________________________________
*     prepare 2-parameter points
*
      NADFS=2                  ! 2-parameter profile analysis
      IPFX=NPSEC(1,ISECA)      ! profile index
      IPFY=NPSEC(2,ISECA)      ! profile index
      CENTEX=AUX(INDAS+IPFX)
      CENTEY=AUX(INDAS+IPFY)
c      SIGMAX=SQRT(AUX(INDVS+(IPF*IPF+IPF)/2))

      CALL DVELLP(AUX(INDVS+IJSYM(IPFX,IPFX)),
     +            AUX(INDVS+IJSYM(IPFX,IPFY)),
     +            AUX(INDVS+IJSYM(IPFY,IPFY)),
     +            SIGMAX,SIGMAY,CPHI,SPHI)      ! cos and sin phi 
c      WRITE(*,*) 'DEVELL',SIGMAX,SIGMAY,CPHI,SPHI
      XLR( 0)=CENTEX           ! original point
      YLR( 0)=CENTEY           ! original point
      FLR( 0)=0.0D0            ! function difference
      IPAK=IPFX
#include "unpackfl.inc"      
      NTLIMX=NTLIM             ! property of X
      IPAK=IPFY
#include "unpackfl.inc"      
      NTLIMY=NTLIM             ! property of Y
      KNOTS=(2*NKNOTS+1)/2     ! make KNOTS even number
      CHSUP=SIGMAS**2          ! chi^2 limit
      NSTAR=12                 ! 15 degree steps
      CRT= COS(PI/DFLOAT(NSTAR))  ! rotation cos
      SRT=-SIN(PI/DFLOAT(NSTAR))  ! rotation sin
      CTH=1.0D0
      STH=0.0D0 
      ISTAR=0

 01   CONTINUE
      FAC=1.1*SIGMAS          ! default factor
      REDL=FAC/FLOAT(KNOTS)   ! one step in units of sigma
      REDR=REDL
c      WRITE(*,*) '(1) REDL REDR',REDR,REDL 
      ILRM=MLR 
      IF(NTLIMX.EQ.1) THEN
         FAC=MIN(FAC,-XLR(0)/(CPHI*CTH*SIGMAX+SPHI*STH*SIGMAY))
      END IF
      IF(NTLIMY.EQ.1) THEN
         FAC=MIN(FAC,-YLR(0)/(-SPHI*CTH*SIGMAX+CPHI*STH*SIGMAY))
      END IF
c      WRITE(*,*) 'FAC ',FAC
      IF(FAC.NE.1.1*SIGMAS) THEN 
         REDL=FAC/FLOAT(KNOTS)
         ILRM=FAC/REDL+0.01       
         ILRM=MIN(MLR,ILRM)
c         WRITE(*,*) 'FAC,KNOTS,REDL,ILRM',FAC,KNOTS,REDL,ILRM
      END IF
C      WRITE(*,*) 'REDR REDL',REDR,REDL
      DO I=1,MLR ! prepare all points
       CLR(+I)=+REDR*DFLOAT(I)  ! unit of parab standard deviations
       CLR(-I)=-REDL*DFLOAT(I)  ! unit of parab standard deviations
       XLR(+I)=XLR(0)+CLR(+I)*(CPHI*CTH*SIGMAX+SPHI*STH*SIGMAY)
       XLR(-I)=XLR(0)+CLR(-I)*(CPHI*CTH*SIGMAX+SPHI*STH*SIGMAY)
       YLR(+I)=YLR(0)+CLR(+I)*(-SPHI*CTH*SIGMAX+CPHI*STH*SIGMAY)
       YLR(-I)=YLR(0)+CLR(-I)*(-SPHI*CTH*SIGMAX+CPHI*STH*SIGMAY)
      END DO
c      WRITE(*,321) 'CLR',CLR
c      WRITE(*,321) 'XLR',XLR
c      WRITE(*,321) 'YLR',YLR
c321  FORMAT(A/(4F15.5)) 
      CHT=CTH                  ! rotation for next vector search
      CTH= CRT*CHT+SRT*STH     ! rotate
      STH=-SRT*CHT+CRT*STH     ! rotate

      ISTAR=ISTAR+1            ! increase rotation index
      ILR = 0                  ! reset index of profile point
      ILRP=-2                  ! direction and step size
*
*     __________________________________________________________________
*     next profile point     
*
 10   IF(ILR*ILRP.GT.0.AND.(FLR(ILR).GT.CHSUP.OR.ABS(ILR).GE.ILRM)) THEN
         IF(ILRP.LT.0) THEN
            ILR1=ILR             ! end-point negative side
            ILRM=MLR             ! max number of steps
         END IF
         IF(ILRP.GT.0) ILR2=ILR  ! end-point positive side
         ILR=ILR+ILRP/2          ! even -> odd
         ILRP=-ILRP              ! reverse direction
      END IF
      ILR=ILR+ILRP               ! next value: change index by step
      CONSTX=XLR(ILR)            ! fixed parameter value
      CONSTY=YLR(ILR)            ! fixed parameter value
c      WRITE(*,*) 'CONSTX CONSTY ',CONSTX,CONSTY
*
*     __________________________________________________________________
*     prepare X(.) and VX(.) 
*
C     CALL XVCOPY(X,VX,AUX(INDXS+1),AUX(INDWM+1)) ! copy to X, VX
c      DO I=1,NX 
c       IF(AUX(INDWM+(I*I+I)/2).NE.0.0D0) THEN ! variance > 0: "measured"
c          X(I)=AUX(INDXS+I)      ! store measured original X back
c       END IF
c      END DO
c      DO I=1,(NX*NX+NX)/2
c       VX(I)=AUX(INDWM+I)        ! store back original VX
c      END DO
      NF=NFPRIM+2                ! number of constraints
c      CALL BPLCON(NX,NF)         ! re-init fit
      CALL APRINT(6,0)           ! suppress printout during profile fit       
*
*     __________________________________________________________________
*     return to fit ....
*
*      ...
      KRET=-1
      RETURN
*
*     __________________________________________________________________
*     back from fit - store chi^2
*
C     ENTRY B2PROF(X,VX,KRET)
      ENTRY B2PROF(KRET)
      FLR(ILR)=CHISQ-CHIREF   ! 20 store chi^2 diff. of secondary fit
c      WRITE(*,*) 'chi^2 calue ',FLR(ILR),' stored ',ILR,ILRP 
      IF(ILR+ILRP.NE.0) GOTO 10  ! continue
*
*     __________________________________________________________________
*     search finished 
*
      NLR=0                      ! reset number of points
      DO I=ILR1,ILR2
       IF(FLR(I).GE.0.0) THEN              ! fit converged
          NLR=NLR+1              ! nr of points
          XC(NLR)=CLR(I)         ! value of variable
          YC(NLR)=FLR(I)         ! function difference
          XCD(NLR)=CLR(I)         ! value of variable
          YCD(NLR)=FLR(I)         ! function difference 
c          WRITE(*,*) NLR,XCD(NLR),YCD(NLR)
       END IF
      END DO
      CALL DCSPLN(XCD,YCD,CSP,NLR) ! spline definition
c      CALL DCSPR (CSP,NLR)         ! print
      DO J=1,3                     ! find intersection points and store
       FAC=DCSPN(0.0D0,DFLOAT(J*J),CSP,NLR,0,+1)
c       WRITE(*,*) 'J FAC ',J,FAC,REDR
       XCONT(ISTAR,J)=XLR(0)+FAC*(XLR(1)-XLR(0))/REDR
       YCONT(ISTAR,J)=YLR(0)+FAC*(YLR(1)-YLR(0))/REDR 
       FAC=DCSPN(0.0D0,DFLOAT(J*J),CSP,NLR,0,-1)
c       WRITE(*,*) 'J FAC ',J,FAC,REDL
       XCONT(ISTAR+NSTAR,J)=XLR(0)-FAC*(XLR(0)-XLR(1))/REDL
       YCONT(ISTAR+NSTAR,J)=YLR(0)-FAC*(YLR(0)-YLR(1))/REDL 
      END DO
      IF(ISTAR.LT.NSTAR) GOTO 01

*     finished: store data 
      DO J=1,3
       CALL ASPROF(IPFX,IPFY,J,FLOAT(J*J),NSTAR+NSTAR,
     +             XCONT(1,J),YCONT(1,J))
      END DO 
      CALL PPROFL                ! print profile  
      KRET=0 
      END 


      SUBROUTINE ASPROF(I1,I2,IV,VALUE,NARG,XD,YD) ! store 
*     store result/get result
#include "nauxfit.inc"
#include "cprofil.inc"
#include "comcfit.inc"
*     I1,I2 = indices of variables
*     IV    = index of coordinate set
*     VALUE = delta F 
c      DOUBLE PRECISION
      REAL XD(*),YD(*)
      REAL             XS(*),YS(*)
*     packinfo-----------------
      DOUBLE PRECISION  DAUX
      REAL              SAUX(2)
      EQUIVALENCE (DAUX,SAUX(1))
*     ...
      N=NARG
c      WRITE(*,*) 'ASPROF ',N,NDENDA,NDENDE,NDACTL
      IND=NDENDE       ! get index
      SAUX(1)=N*128+IV ! pack
      SAUX(2)=VALUE
      AUX(IND+1)=DAUX  ! first control word
      SAUX(1)=I1
      SAUX(2)=I2
      AUX(IND+2)=DAUX  ! second control word
      IND=IND+2
c      WRITE(*,*) 'AUX(IND)=DAUX',AUX(IND),DAUX
      DO I=1,N
       SAUX(1)=XD(I)
       SAUX(2)=YD(I)
c       WRITE(*,*) I,XD(I),YD(I)
       AUX(IND+I)=DAUX ! coordinates
      END DO
      NDENDE=IND+N     ! N coordinate pairs
      SAUX(1)=0
      SAUX(2)=0
      AUX(NDENDE+1)=DAUX
      RETURN

      ENTRY AGPROF(I1,I2,IV,VALUE,NARG,XS,YS) ! get 
      NARG=0
      IND=NDENDA
 10   IF(IND.GE.NDENDE) RETURN 
      DAUX=AUX(IND+1)
      N    =SAUX(1) ! number of xy-pairs
      IVI=MOD(N,128) ! index of value
      N=N/128       ! number of words
      VALUEI=SAUX(2) ! value
      IF(N.EQ.0) RETURN 
      DAUX=AUX(IND+2)
      II1=SAUX(1)    ! index of variable
      II2=SAUX(2)    ! index of variable
      IND=IND+2
      IF(II1.NE.I1.OR.II2.NE.I2.OR.IVI.NE.IV) THEN
         IND=IND+N
         GOTO 10
      END IF
      VALUE=VALUEI
      DO I=1,N
       DAUX=AUX(IND+I)
       XS(I)=SAUX(1)
       YS(I)=SAUX(2)
      END DO 
      NARG=N 
      RETURN

      ENTRY AGPRAL ! print all
      WRITE(*,*) '_____________________________________________'
      WRITE(*,*) 'AGPRAL = print all '
      IND=NDENDA
 01   IF(IND.GE.NDENDE) THEN
         WRITE(*,*) '_____________________________________________'
         RETURN
      END IF
      DAUX=AUX(IND+1)
      N  =SAUX(1)   ! number of xy-pairs
      IVI=MOD(N,128) ! index of value
      N=N/128       ! number of words
      IF(N.EQ.0) RETURN
      VALUEI=SAUX(2) ! value
      DAUX=AUX(IND+2)
      II1=SAUX(1)    ! index of variable 
      II2=SAUX(2)    ! index of variable
      IND=IND+2
*
      WRITE(*,*) ' '
      WRITE(*,*) 'Profile ',N,IVI,VALUEI,II1,II2
      WRITE(*,101)
      DO I=1,N
       DAUX=AUX(IND+I)
       WRITE(*,102) I,SAUX(1),SAUX(2)
      END DO 
      IND=IND+N
      GOTO 01
 101  FORMAT(' index       X-profile     Y-profile')
 102  FORMAT(I6,4X,2G14.5)
      END

      SUBROUTINE PPROFL     ! print profile results
#include "comcfit.inc"
#include "cprofil.inc"
      CHARACTER*1 C, PANAME*16  
      DOUBLE PRECISION PCD(6),SCD(6),FDL,FDR,FD,DCSPN ! DFD(6)
      DOUBLE PRECISION XFDL,XFDR 
c      DATA DFD/ 1.0D0,2.706D0,3.84D0,6.63D0,7.88D0,10.82D0/
      DATA PCD/68.27D0, 90.0D0,95.0D0, 99.D0,99.5D0,99.9D0/
      DATA SCD/ 1.0D0,1.645D0,1.96D0,2.575D0,2.807D0,3.29D0/ 
*     ...
      WRITE(*,*) ' '
      WRITE(*,*) ('_',J=1,71)
      IF(NADFS.EQ.1) THEN ! 1-dim profile 
         CALL APNAME(-IPF,PANAME)
         WRITE(*,102) 'Result of 1-dim profile analysis for parameter',
     +               IPF,':',PANAME
         WRITE(*,*) ' '
         WRITE(*,*) '    f(X) = chi^2 as a function',
     +                 ' of X = parameter value' 
         WRITE(*,*) ' ' 
         CALL DCSPR (CSP,NLR) ! print
         WRITE(*,*) ' '
         CALL DCSFPL(CSP,NLR) ! plot

         WRITE(*,*) ' '
         WRITE(*,102) 'Confidence intervals for parameter',
     +                 IPF,':',PANAME
         WRITE(*,*) ' '
         WRITE(*,*) '      Prob   sigmas ',
     +      '_____________________________________',
     +      '   sigma units'
         WRITE(*,105)  0,CENTER-SIGMAX,CENTER,CENTER+SIGMAX
 105     FORMAT(I4,15X,G13.5,G13.5,G13.5)
c106     FORMAT(A,1X,A16,':',I3) 
         DO I=1,6    ! find interpolated points
c          FD=DFD(I)
          FD=SCD(I)**2
          XFDL=DCSPN(XLR(0),FD,CSP,NLR,0,-1)
          XFDR=DCSPN(XLR(0),FD,CSP,NLR,0,+1)
          IF(XFDL.LT.XLR(0).AND.XFDR.GT.XLR(0)) THEN
             FDL=(XLR(0)-XFDL)/SIGMAX
             FDR=(XFDR-XLR(0))/SIGMAX 
             WRITE(*,104) I,PCD(I),SCD(I),XFDL,XFDR,FDL,FDR
          END IF 
 104     FORMAT(I4,F7.1,' %',F6.2,G13.5,5X,'...',5X,G13.5,2F7.2)
         END DO
      ELSE IF(NADFS.EQ.2) THEN ! 2-dim profile
         CALL APNAME(-IPFX,PANAME)
         WRITE(*,102) 'Result of 2-dim profile analysis for parameter',
     +               IPFX,':',PANAME
         CALL APNAME(-IPFY,PANAME)
         WRITE(*,102) '                                 and parameter', 
     +               IPFY,':',PANAME

         WRITE(*,*) ' '
         WRITE(*,*) 'Plot of profile:'
         WRITE(*,*) ' '
         CALL FPS
         DO J=1,3
          IF(J.EQ.1) C='1'
          IF(J.EQ.2) C='2'
          IF(J.EQ.3) C='3' 
          DO I=1,NSTAR+NSTAR
           CALL FPNXY(C,1,XCONT(I,J),YCONT(I,J))
          END DO
         END DO
         CALL FPL
         WRITE(*,*) ' '
      END IF
      WRITE(*,*) ('_',J=1,71)
 102  FORMAT(1X,A,I3,A,A)
      END
 

      SUBROUTINE APROF1(I1,XA,XB)
*     prepare profile data for variable I1  
      REAL XS(101),YS(101)
      DOUBLE PRECISION DIERFC,XCR  ! CLS
      DOUBLE PRECISION CSP(5,73),XD,DCSPF,XZERO,CH2,XA,XB,X,Y,XCL,CL
      DOUBLE PRECISION XCD(*),YCD(*),XCEN(3,3),CL2
      DOUBLE PRECISION CLA,CLB,CLA1,CLB1,CLA3,CLB3
      LOGICAL READY
      DATA READY/.FALSE./ 
*     ...
      READY=.FALSE.
      XA=0.0D0
      XB=0.0D0
      CALL AGPROF(I1,0,0,VALUE,NLR,XS,YS) ! get data from memory
      IF(NLR.EQ.0) RETURN
      READY=.TRUE.
      XA=XS(  1)
      XB=XS(NLR)
      WRITE(*,*) 'XA XB',XA,XB 
      IM=1                            ! index of ...
      DO I=1,NLR
       CSP(1,I)=YS(I)
       CSP(5,I)=XS(I)
       IF(CSP(1,I).LT.CSP(1,IM)) IM=I ! ... smallest Y value 
      END DO
      XZERO=XS(IM)                    ! X at smallest Y value
      CALL ECSPLN(CSP,NLR)            ! spline definition
      RETURN

      ENTRY APROFJ(XCD,YCD,NCD)       ! intern 
      NLR=NCD  
      CALL DCSPLN(XCD,YCD,CSP,NLR)    ! spline definition
      IM=1                            ! index of ...
      DO I=1,NLR
       IF(CSP(1,I).LT.CSP(1,IM)) IM=I ! ... smallest Y value
      END DO
      XZERO=CSP(5,IM)                 ! X at smallest Y value 
      READY=.TRUE.
      RETURN

      ENTRY APRFUN(X,Y)               ! return Y for given X
      Y=0.0D0
      IF(.NOT.READY) RETURN
      IF(X.LT.XS(1).OR.X.GT.XS(NLR)) RETURN
      XD=X
      Y=DCSPF(XD,CSP,NLR)
      RETURN 
      
      ENTRY APROCL(XCL,CL)            ! return X for given conlev 
*     return x for confidence level CL
*      CL < 0         left limit  
*      0  < CL < 0.5  right limit
      XCL=0.0D0
      IF(.NOT.READY) RETURN
      IF(CL.LT.0.0) THEN
         CH2=2.0D0*DIERFC(2.0D0*ABS(CL))**2
         XCL=DCSPN(XZERO,CH2,CSP,NLR,0,-1)
         IF(XCL.EQ.XZERO)   XCL=0.0D0
      ELSE
         CH2=2.0D0*DIERFC(2.0D0*CL)**2 
         XCL=DCSPN(XZERO,CH2,CSP,NLR,0,+1)
         IF(XCL.EQ.XZERO) XCL=0.0D0
      END IF
      RETURN

      ENTRY APRCEN(XCL,XCR,CL)          ! return XL,XR for given CL
*     0.5 < CL        central limit 
      XCL=0.0D0
      XCR=0.0D0
      IF(.NOT.READY) RETURN

c      CL=0.6827 ! force 1 sigma

      CL2=2.0*(1.0D0-CL) 
       
      CLA=0.5*CL2 ! left  C.L.
      CLB=CL2-CLA ! right C.L.
      FACT=0.9D0
c      WRITE(*,*) 'CLA DIERFC(CLA)',CLA,DIERFC(CLA)
c      WRITE(*,*) 'CLB DIERFC(CLB)',CLB,DIERFC(CLB)
     
c      WRITE(*,*) 'CLA Chi^2 value =',CLA,2.0D0*DIERFC(CLA)**2
c      WRITE(*,*) 'CLA Chi^2 value =',CLB,2.0D0*DIERFC(CLB)**2


      DO ITER=1,100

      XCL=DCSPN(XZERO,2.0D0*DIERFC(CLA)**2,CSP,NLR,0,-1)
      XCR=DCSPN(XZERO,2.0D0*DIERFC(CLB)**2,CSP,NLR,0,+1)
c      WRITE(*,*) 'APRCEN: XZERO',CL,XCL,XCR,' <<'
      XCEN(2,1)=XCL
      XCEN(2,2)=XCR
      XCEN(2,3)=XCEN(2,2)-XCEN(2,1)

        
      CLA1=CLA*FACT
      CLB1=CL2-CLA1  

      XCEN(1,1)=DCSPN(XZERO,2.0D0*DIERFC(CLA1)**2,CSP,NLR,0,-1)
      XCEN(1,2)=DCSPN(XZERO,2.0D0*DIERFC(CLB1)**2,CSP,NLR,0,+1)   
      XCEN(1,3)=XCEN(1,2)-XCEN(1,1)
      CLA3=CLA/FACT
      CLB3=CL2-CLA3
      XCEN(3,1)=DCSPN(XZERO,2.0D0*DIERFC(CLA3)**2,CSP,NLR,0,-1)
      XCEN(3,2)=DCSPN(XZERO,2.0D0*DIERFC(CLB3)**2,CSP,NLR,0,+1)   
      XCEN(3,3)=XCEN(3,2)-XCEN(3,1)
c      WRITE(*,*) 'XCEN',XCEN

      IMIN=1
      DO I=1,3
       IF(XCEN(I,3).LT.XCEN(IMIN,3)) IMIN=I
      END DO
c      WRITE(*,*) ' minimum is ',XCEN(IMIN,3)
      XCL=XCEN(IMIN,1)
      XCR=XCEN(IMIN,2)

      IF(IMIN.EQ.1) THEN         ! to minus side
         CLA=CLA*FACT
      ELSE IF(IMIN.EQ.2) THEN    ! reduce factor
         FACT=0.5D0*(FACT+1.0D0) 
      ELSE IF(IMIN.EQ.3) THEN    ! to plus side 
         CLA=CLA/FACT
      END IF
      CLB=CL2-CLA                ! new second limit
c     WRITE(*,*) 'CLAB is',CLA,CLB,FACT
      IF(1.0D0-FACT.LT.0.001) GOTO 10
      END DO
 10   CONTINUE
c      WRITE(*,*) 'CLAB is',CLA,CLB,FACT,XCL,XCR,ITER
      END 

      SUBROUTINE PPROF1     ! print 1-dim profile results
#include "comcfit.inc"
#include "cprofil.inc"
      CHARACTER PANAME*16 ! C*1
      DOUBLE PRECISION PCD(8),SCD(8),FDL,FDR,FD,DCSPN
      DOUBLE PRECISION XFDL,XFDR,CL
      DOUBLE PRECISION CLIM(8)
*     percentage values 
      DATA PCD/68.27D0,90.0D0,95.00D0,95.45D0,
     +         99.00D0,99.5D0,99.73D0,99.90D0/
*     standard deviations
      DATA SCD/ 1.000D0,1.645D0,1.96D0,2.00D0,
     +          2.575D0,2.807D0,3.00D0,3.29D0/
*     C.L. for limits
      DATA CLIM/0.10D0,0.05D0,0.02D0,0.01D0,0.005D0,0.002D0,0.001D0,
     +          0.0005D0/ 
*     ...
      WRITE(*,*) ' '
      WRITE(*,*) ('_',J=1,71)

      CALL APNAME(-IPF,PANAME)
      WRITE(*,102) 'Result of 1-dim profile analysis for parameter',
     +            IPF,':',PANAME
      WRITE(*,*) ' '
      WRITE(*,*) '    f(X) = chi^2 as a function',
     +           ' of X = parameter value'
      WRITE(*,*) ' '
      CALL DCSPR (CSP,NLR) ! print
      WRITE(*,*) ' '
      CALL DCSFPL(CSP,NLR) ! plot

      WRITE(*,*) ' '
      WRITE(*,102) 'Confidence intervals for parameter',
     +              IPF,':',PANAME
      WRITE(*,*) ' '
      WRITE(*,*) '      Prob   sigmas ',
     +   '_____________________________________',
     +   '   sigma units'
      WRITE(*,105)  0,CENTER-SIGMAX,CENTER,CENTER+SIGMAX
      DO I=1,8    ! find interpolated points
       FD=SCD(I)**2
       XFDL=DCSPN(XLR(0),FD,CSP,NLR,0,-1)
       XFDR=DCSPN(XLR(0),FD,CSP,NLR,0,+1)
       IF(XFDL.LT.XLR(0).AND.XFDR.GT.XLR(0)) THEN
          FDL=(XLR(0)-XFDL)/SIGMAX
          FDR=(XFDR-XLR(0))/SIGMAX
          WRITE(*,104) I,PCD(I),SCD(I),XFDL,XFDR,FDL,FDR
         
c          CL=1.0D0-0.01D0*PCD(I)
c          WRITE(*,*) 'Try CL=',CL 
c          XFDL=0.0
c          XFDR=0.0
c          CALL APRCEN(XFDL,XFDR,0.5D0*CL)
c          WRITE(*,104) I,PCD(I),SCD(I),XFDL,XFDR,FDL,FDR

       END IF
      END DO

      WRITE(*,*) ' '
      WRITE(*,*) '          C.L.    lower limit          ',
     +           '      C.L.    upper limit'
      DO I=1,8
       CL=CLIM(I)
       CALL APROCL(XFDL,-CL)
       CALL APROCL(XFDR,+CL)
       FD=100.0D0*(1.0D0-CL)
       IF(XFDL.NE.0.0D0) THEN
          IF(XFDR.NE.0.0D0) THEN
             WRITE(*,114) I,FD,'<',XFDL,FD,'>',XFDR 
          ELSE
             WRITE(*,114) I,FD,'<',XFDL
          END IF
       ELSE IF(XFDR.NE.0.0D0) THEN
          WRITE(*,214) I,FD,'>',XFDR
       END IF 
c       WRITE(*,114) I,FD,'<',XFDL,FD,'>',XFDR
      END DO 
 114  FORMAT(I4,5X,F7.2,' %',4X,A1,G11.5,10X,F7.2,' %',4X,A1,G11.5) 
 214  FORMAT(I4,40X,F7.2,' %',4X,A1,G11.5)

      WRITE(*,*) ' '
      WRITE(*,102) 'Minimum length confidence intervals for parameter',
     +              IPF,':',PANAME
      WRITE(*,*) ' '
      WRITE(*,*) '      Prob   sigmas ',
     +   '_____________________________________',
     +   '   sigma units'
      WRITE(*,105)  0,CENTER-SIGMAX,CENTER,CENTER+SIGMAX
      DO I=1,8    ! find interpolated points
          CL=1.0D0-0.01D0*PCD(I)
          XFDL=0.0
          XFDR=0.0
c          CALL APRCEN(XFDL,XFDR,0.5D0*CL)
c          WRITE(*,104) I,PCD(I),SCD(I),XFDL,XFDR !,FDL,FDR
          CALL APRCEN(XFDL,XFDR,0.01D0*PCD(I)) 
          WRITE(*,104) I,PCD(I),SCD(I),XFDL,XFDR,(XFDR-XFDL)/SIGMAX 
      END DO

      WRITE(*,*) ('_',J=1,71)
 102  FORMAT(1X,A,I3,A,A)
 104  FORMAT(I4,F7.1,' %',F6.2,G13.5,5X,'...',5X,G13.5,2F7.2)
 105  FORMAT(I4,15X,G13.5,G13.5,G13.5)   
c106  FORMAT(A,1X,A16,':',I3)
      END

      FUNCTION DIERFC(Y) ! INVERSE OF ERROR FUNCTION IN DOUBLE PRECISION
      IMPLICIT DOUBLE PRECISION (A - H, O - Z)
      PARAMETER (
     +    QA = 9.16461398268964D-01, 
     +    QB = 2.31729200323405D-01, 
     +    QC = 4.88826640273108D-01, 
     +    QD = 1.24610454613712D-01, 
     +    Q0 = 4.99999303439796D-01, 
     +    Q1 = 1.16065025341614D-01, 
     +    Q2 = 1.50689047360223D-01, 
     +    Q3 = 2.69999308670029D-01, 
     +    Q4 = -7.28846765585675D-02)
      PARAMETER (
     +    PA = 3.97886080735226000D+00, 
     +    PB = 1.20782237635245222D-01, 
     +    P0 = 2.44044510593190935D-01, 
     +    P1 = 4.34397492331430115D-01, 
     +    P2 = 6.86265948274097816D-01, 
     +    P3 = 9.56464974744799006D-01, 
     +    P4 = 1.16374581931560831D+00, 
     +    P5 = 1.21448730779995237D+00, 
     +    P6 = 1.05375024970847138D+00, 
     +    P7 = 7.13657635868730364D-01, 
     +    P8 = 3.16847638520135944D-01, 
     +    P9 = 1.47297938331485121D-02, 
     +    P10 = -1.05872177941595488D-01, 
     +    P11 = -7.43424357241784861D-02)
      PARAMETER 	(
     +    P12 = 2.20995927012179067D-03, 
     +    P13 = 3.46494207789099922D-02, 
     +    P14 = 1.42961988697898018D-02, 
     +    P15 = -1.18598117047771104D-02, 
     +    P16 = -1.12749169332504870D-02, 
     +    P17 = 3.39721910367775861D-03, 
     +    P18 = 6.85649426074558612D-03, 
     +    P19 = -7.71708358954120939D-04, 
     +    P20 = -3.51287146129100025D-03, 
     +    P21 = 1.05739299623423047D-04, 
     +    P22 = 1.12648096188977922D-03)
      Z = Y
      IF (Y .GT. 1) Z = 2 - Y
      W = QA - LOG(Z)
      U = SQRT(W)
      S = (QC + LOG(U)) / W
      T = 1 / (U + QB)
      X = U * (1 - S * (0.5D0 + S * QD)) - 
     +    ((((Q4 * T + Q3) * T + Q2) * T + Q1) * T + Q0) * T
      T = PA / (PA + X)
      U = T - 0.5D0
      S = (((((((((P22 * U + P21) * U + P20) * U + 
     +    P19) * U + P18) * U + P17) * U + P16) * U + 
     +    P15) * U + P14) * U + P13) * U + P12
      S = ((((((((((((S * U + P11) * U + P10) * U + 
     +    P9) * U + P8) * U + P7) * U + P6) * U + P5) * U + 
     +    P4) * U + P3) * U + P2) * U + P1) * U + P0) * T - 
     +    Z * EXP(X * X - PB)
      X = X + S * (1 + X * S)
      IF (Y .GT. 1) X = -X
      DIERFC = X
      END



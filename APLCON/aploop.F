 
      SUBROUTINE APLCON(NVAR,MCST)                ! dimension parameters
*     ==================================================================
*     initialize and define dimension  NVAR/MCST
*     set default parameters
*     define pointer to arrays within aux array
*     clear aux array 
*     ==================================================================
      IMPLICIT NONE 
      INTEGER NVAR,MCST,IJ  ! ,I1,I2,NFF,I
#include "comcfit.inc"
#include "nauxfit.inc"
#include "cprofil.inc"
      LOGICAL START
      DATA START/.TRUE./
*     ...
      IF(START) THEN
         START=.FALSE.
         NCASE=0          ! reset counter
         IPR=5
      END IF 
      NCASE=NCASE+1       ! count cases

      NSECA=0             ! reset number of profile searches
      NFADD=0             ! reset number of additional constraints
      NADFS=0             ! primary fit initialization
                          ! using space for NF+2 instead of NF

      NX    =NVAR         ! number of variables
      NF    =MCST         ! number of constraint equations
      NFPRIM=NF           ! primary value of NF 
      NDPDIM=NAUX         ! dimension of AUX array
      LUNSIM=6            ! printout unit (default)
c     IPR=5               ! default print flag is 5
      DERFAC=1.0D-3       ! derivative factor
      DERUFC=1.0D-5       ! factor or unmeasured variable
      DERLOW=1.0D-2       ! factor for lower limit
      EPSF  =1.0D-6       ! accuracy limit
      EPSCHI=0.1D-4       ! chi2 accuracy limit      
      ITERMX=100          ! max number of iterations
      ITERMX=50
      ITERMX=20
      ITERMX=10 
      NAUXC=NAUX          ! copy AUX dimension

      INIT  =0
      ISTAT =0            ! init phase

      NXF=NX+NF           ! total number of fit equations
      MXF=(NXF*NXF+NXF)/2 ! number elements symmetric matrix

* ______________________________________________________________________
*     Indices for steps, flags and limits are defined here 
                          ! reserve NX * (NF + 2) for Jacobian
      INDST=NX*(NF+2)     ! steps
      INDTR=INDST+NX      ! transformation flags
      INDLM=INDTR+NX      ! 2*limits for variables
      NDTOT=INDLM+2*NX    ! space used so far
      IF(NDTOT.GT.NAUX) THEN
         CALL APRINI(1)   ! error print
         STOP
      END IF
*     __________________________________________________________________
*     storage of initial sub-arrays
*           1 ...    NX * (NF+2)   Jacobian, derivative matrix   A(.)
*     INDST+1 ...    NX            steps                         ST(.)
*     INDTR+1 ...    NX            properties of variables       
*     INDLM+1 ...    2 * NX        limits for variables          XL(2,.)
*                    ----------- 
*     NDTOT =        NX * (NF+6)   initial memory area
      DO IJ=1,NDTOT       ! NX*(NF+6)
       AUX(IJ)=0.0D0      ! clear  A(.), ST(.),...,XL(2,.)
      END DO
      NDF   =NF           ! reset n d f
      NCALLS=0            ! reset number of calls

      CALL APRINI(0)      ! initial print without X, VX  
c      CALL APNAME(0,' ')  ! reset parameter names
      END 

      SUBROUTINE APLOOP(X,VX,F, IRET)        ! steering routine for loop
      INTEGER IRET
      DOUBLE PRECISION X(*),VX(*),F(*)
#include "nauxfit.inc"
#include "comcfit.inc"
#include "cprofil.inc"

      IF(NCALLS.NE.0) GOTO 10
*     __________________________________________________________________
*     indices/pointer etc at first APLOOP entry 
      NFF=NF+NFADD               ! max. number of constraints
      NXF=NX+NFF                 ! total number of fit equations
C      WRITE(*,*) 'NFF,NXF=',NFF,NXF
      MXF=(NXF*NXF+NXF)/2        ! number elements symmetric matrix
      INDFC=INDLM+2*NX           ! pointer to FC(NF) = copy of F(NF)
      INDCF=INDFC+NFF            ! pointer to FCOPY 
      INDHH=INDCF+NFF            ! pointer to HH(NF) = copy of F(.)
      INDXS=INDHH+NFF            ! save X(.)        pointer
      INDDX=INDXS+NX             ! step
      INDXP=INDDX+NX             ! previos step
      INDRH=INDXP+NX             ! right-hand side
      INDWM=INDRH+NXF            ! weight matrix
      INDIA=INDWM+MXF            ! matrix diagonal "DIAG"
      INDQN=INDIA+NXF            ! next pointer    "QNEXT"
      INDAS=INDQN+NXF            ! X result
      INDPU=INDAS+NX             ! pulls, solution X and Vx
      NDTOT=INDPU+NX             ! total number of words (so far)  
      NDTOTL=NDTOT
C      WRITE(*,*) 'NDTOTL,NAUX=',NDTOTL,NAUX
      IF(NDTOTL.GT.NAUX) THEN
C         WRITE(*,*) 'NDTOTL,NAUX=',NDTOTL,NAUX
         CALL AIPRIN(X,VX,1,IRET) ! error printout
         STOP
      END IF
      DO I=INDLM+2*NX+1,NDTOT
       AUX(I)=0.0D0          ! reset part of aux, unused so far
      END DO
      CALL ASTEPS(X,VX, AUX(1+INDST)) ! initial steps ST(.)
      IF(IPR.NE.0) THEN
         WRITE(*,*) 'Printout of initial values and correlation coeffs:' 
         CALL CIPRV (LUNSIM,X,VX,NX) 
         CALL CFCORR(LUNSIM,VX,NX)
      END IF 
*     __________________________________________________________________
*     internal APLOOP    
 10   IRET=-1                ! default status is -1 = continue
      CALL IPLOOP(X,VX,F, AUX(INDXS+1),AUX(INDDX+1),
     +                    AUX(INDCF+1),AUX(INDXP+1),AUX(INDRH+1),  IRET)
      END 

      SUBROUTINE IPLOOP(X,VX,F, XS,DX,FCOPY,XP,RH, IRET)
     +                                       ! steering routine for loop
*     ==================================================================
*     initial print
*        add 1/2 constraints for profile analysis
*     call IPLDER
*     call IPLCON
*        profile analysis   
*     ==================================================================
      IMPLICIT NONE
      INTEGER J,IRET,JRET,NFIT,KRET,IPRSAV
c      INTEGER NITER,NFIT
c      INTEGER J,IRET,JRET,NSECAS,IJSYM,ILRP,ILR1,ILR2,NFUN ,NN,NTLIMP
      DOUBLE PRECISION X(*),VX(*),F(*) ! ,FOPT,FAC
      DOUBLE PRECISION XS(*),DX(*),FCOPY(*),XP(*),RH(*),FEX(100)
*     local variables
      DOUBLE PRECISION FJ,FADD(2)        ! constraint values
#include "comcfit.inc" 
#include "nauxfit.inc"
#include "cprofil.inc"
C#include "declarefl.inc"
c      REAL XRP,YRP
c      DOUBLE PRECISION DCSPN
*     __________________________________________________________________
      INTEGER ISTATU
*                    = -1  numerical derivatives
*                    =  0  constraint function evaluation
*                    =  1  constraint function with test afterwards
*                    =  2  end-of-fit
*     __________________________________________________________________
*     ... 
      NCALLS=NCALLS+1                    ! count calls

*     __________________________________________________________________
*     initialization 
      IF(NCALLS.NE.1) GOTO 20
      ISTATU=0                                                       !!!
      NFIT=0                             ! reset fit count
      ITER=0
      NADFS=0                            ! no additional constraints
      ISECA=0                            ! index of profile analysis 
      IF(IPR.GE.3) CALL AIPRIN(X,VX,1,IRET)   ! initial printout
      DO J=1,NX
       XS(J)=X(J)                        ! save initial X values
       DX(J)=0.0D0                       ! reset correction DX
      END DO
*     __________________________________________________________________
*     start/restart 
 10   NFIT=NFIT+1                        ! count fits
c      IF(NFIT.GE.10) STOP
      ITER=0                        
      NCST=0
      CHISQ=0.0D0 
c      WRITE(*,*) 'Start/restart fit NFIT ISTATU ',NFIT,ISTATU
c      WRITE(*,*) 'NF NX NADFS ',NF,NX,NADFS
c      WRITE(*,321) 'Vx',(VX(J),J=1,6)
c      WRITE(*,321) 'X ',(X(J),J=1,3) 
c      WRITE(*,321) 'XS',(XS(J),J=1,3)
c      WRITE(*,321) 'DX',(DX(J),J=1,3)  
c      WRITE(*,*)   'CONSTR IPF',CONSTR,IPF
c      WRITE(*,321) 'F',(F(J),J=1,NF) 
c      WRITE(*,*) 'INDST ',INDST
c      WRITE(*,321) 'ST',(AUX(INDST+J),J=1,3)
*     __________________________________________________________________
*     add constraints for profile analysis
 20   IF(NADFS.EQ.1) THEN                ! add constraints for profiles
         FADD(1)=CONSTR-X(IPF)             ! 1-parameter constraint
c         WRITE(*,*) 'Extra constraint',FADD(1),CONSTR,X(IPF),IPF
      ELSE IF(NADFS.EQ.2) THEN 
         FADD(1)=CONSTX-X(IPFX)          ! 2-parameter constraint
         FADD(2)=CONSTY-X(IPFY)
      END IF
      DO J=1,NF
       IF(NADFS.EQ.1.AND.J.EQ.NF) THEN
          FJ=FADD(1)
       ELSE IF(NADFS.EQ.2.AND.J+1.GE.NF) THEN
          FJ=FADD(J-NF+2)
       ELSE
          FJ=F(J)
       END IF
       FEX(J)=FJ   ! extended F
      END DO 
*     __________________________________________________________________
*     constraint test summary
c      IF(NFIT.GT.1) WRITE(*,*) 'ISTATU=',ISTATU
      IF(ISTATU.LT.0) GOTO 40
      FTESTP=FTEST                       ! save previous value
      FTEST=0.0D0                        ! reset constraint tests 
      FRMSP=FRMS
      FRMS =0.0D0                                          
      DO J=1,NF
       IF(NADFS.EQ.1.AND.J.EQ.NF) THEN
          FJ=FADD(1)
       ELSE IF(NADFS.EQ.2.AND.J+1.GE.NF) THEN
          FJ=FADD(J-NF+2)
       ELSE
          FJ=F(J)
       END IF
       FCOPY(J)=FJ                       ! copy constraint vector
       FTEST=FTEST+ABS(FJ)               ! sum absolute values
       FRMS =FRMS +FJ**2                 ! sum squares
      END DO
c      WRITE(*,*) 'ISTATU=',ISTATU
c      WRITE(*,321) 'X    ',(X    (J),J=1,NX)
c      WRITE(*,321) 'FCOPY',(FCOPY(J),J=1,NF)
c321  FORMAT(A/(4F15.6)) 
      FTEST=MAX(1.0D-16,FTEST/FLOAT(NF)) ! average |F|
      FRMS=SQRT(FRMS/FLOAT(NF)+1.0D-32)  ! LS mean 
      IF(ITER.EQ.0) CALL AIPRIN(X,VX,3,IRET)
c      IF(NFIT.GT.1) WRITE(*,*) 'ISTATU=',ISTATU,FTEST
      IF(ISTATU.EQ.1) GOTO 60
*     __________________________________________________________________
*     start numerical derivatives
 30   ISTATU=-1                                                      !!!
c      IF(NFIT.GT.1) WRITE(*,*) 'start num'
*     __________________________________________________________________
*     derivative calculation
 40   IF(ISTATU+1.NE.0) GOTO 50          ! ISTATU=-1 
      CALL ANUMDE(X,FEX, AUX,              ! derivative matrix A
     +                    AUX(INDST+1),     ! steps  ST(.)
     +                    AUX(INDLM+1),     ! limits XL(2,.)
     +                    AUX(1+INDFC),     ! copy FC(.) central F(.)
     +                    AUX(1+INDHH),     ! copy HH(.) shifted F(.)
     +                    JRET)
c      IF(NFIT.GT.1) WRITE(*,*) 'do num',JRET
      IRET=-1
      IF(JRET.LT.0) RETURN            !...for constraint calculation
*     __________________________________________________________________
*     add elements of Jacobian for profile analysis
 50   IF(NADFS.EQ.1) THEN                ! 1-parameter derivative matrix
         DO J=1,NX
          AUX(J+NX*(NF-1))=0.0D0
         END DO
         AUX(IPF+NX*(NF-1))=-1.0D0
      ELSE IF(NADFS.EQ.2) THEN           ! 2-parameter derivative matrix
         DO J=1,NX
          AUX(J+NX*(NF-2))=0.0D0
          AUX(J+NX*(NF-1))=0.0D0
         END DO
         AUX(IPFX+NX*(NF-2))=-1.0D0
         AUX(IPFY+NX*(NF-1))=-1.0D0
      END IF
*     __________________________________________________________________
*     next iteration
      CALL ANITER(X,VX,FCOPY,
     +            AUX,XP,RH,AUX(INDWM+1),DX)
c      WRITE(*,*) 'ANITER ',(AUX(INDST+J),J=1,NX)
      GOTO 70  
*     __________________________________________________________________
*     test cutsteps
 60   CALL ANTEST(IRET)
      IF(IRET+1.EQ.0) GOTO 30   ! numerical derivative:   ISTATU=-1
      IF(IRET  .GE.0) GOTO 80   ! convergence or failure: ISTATU= 2
*     __________________________________________________________________
*     apply corrections DX(.) to X(.) with transformations
 70   CALL ADDTOX(X,XS,DX,XP)  
      ISTATU=1                            ! test at next entry
      RETURN
*     __________________________________________________________________
*     end-of-primary-fit (NFIT=1)
 80   ISTATU=2
c      WRITE(*,*) 'statement 80 in IPLOOP, NFIT=',NFIT
      IF(NFIT.NE.1) GOTO 82
c      WRITE(*,*) 'INDDX,INDAS,INDWM,INDPU',INDDX,INDAS,INDWM,INDPU
      CALL ACOPXV(X,VX,AUX(INDDX+1),AUX(INDAS+1),
     +                 AUX(INDWM+1),AUX(INDPU+1))
c      WRITE(*,*) 'vor AIPRIN'
      CALL AIPRIN(X,VX,4,IRET)                        ! final print
c      WRITE(*,*) '>>>NSECA = ',NSECA
      IF(NSECA.EQ.0) RETURN
*     __________________________________________________________________
*     start of profile analysis
      INDVS=INDPU+NX                ! index of save Vx result
      NDTOT=INDVS+(NX*NX+NX)/2      ! total number of words (so far)
      NDENDE=NDTOT
      NDENDA=NDENDE
      DO J=1,(NX*NX+NX)/2
       AUX(INDVS+J)=VX(J)    ! saved primary cov. matrix
       VX(J)=AUX(INDWM+J)        ! store back original VX
      END DO
      CHIREF=CHISQ
      IPRSAV=IPR 
c      WRITE(*,*) 'Original VX '
c      WRITE(*,*) (VX(J),J=1,6)
c      WRITE(*,*) 'Copy INDVS ',INDVS,' fitted matrix'
c      WRITE(*,*) (AUX(INDVS+J),J=1,6)
*     prepare storage for secondary fits
*     __________________________________________________________________
*     start next profile analysis
 81   IF(ISECA.GE.NSECA) GOTO 99
      ISECA=ISECA+1           ! next profile analysis
      IF(NPSEC(2,ISECA).EQ.0) THEN
         NADFS=1 ! 1-parameter profile analysis
         NF=NFPRIM+NADFS
c         WRITE(*,*) 'Calling A1PROF'
c         CALL A1PROF(X,VX,KRET)
         CALL A1PROF(KRET)
      ELSE
         NADFS=2 ! 2-parameter profile analysis
         NF=NFPRIM+NADFS
c        CALL A2PROF(X,VX,KRET)
         CALL A2PROF(KRET)
      END IF 
      ISTATU=0
      GOTO 10 ! restart fit
*     __________________________________________________________________
*     secondary fit for profile analysis ended, store result
 82   IF(NADFS.EQ.1) THEN
c      WRITE(*,321) 'ST B1PROF V',(AUX(INDST+J),J=1,3)
c        CALL B1PROF(X,VX,KRET)          ! 1-parameter profile analysis
         CALL B1PROF(KRET)          ! 1-parameter profile analysis
c          WRITE(*,321) 'ST B1PROF N',(AUX(INDST+J),J=1,3) 
      ELSE 
c        CALL B2PROF(X,VX,KRET)         ! 2-parameter profile analysis
         CALL B2PROF(KRET)         ! 2-parameter profile analysis
      END IF
      ISTATU=0 
      IF(KRET.LT.0) GOTO 10              ! next fit
      GOTO 81 
 99   CHISQ=CHIREF                       ! store primary fit results
      IPR=IPRSAV                         ! printflag
      DO J=1,NX
       X(J)=AUX(INDAS+J)
      END DO
      DO J=1,(NX*NX+NX)/2
       VX(J)=AUX(INDVS+J)
      END DO  
      CALL AIPRIN(X,VX,4,IRET)                ! final print
      WRITE(*,*) ('_',J=1,71)
      WRITE(*,*) ('_',J=1,71)
      IRET=0
      END 


      SUBROUTINE ASTEPS(X,VX,ST)         ! define initial steps
*     ==================================================================
*     check transformed variables (e.g. > 0)
*     define initial step size for
*        o  measured variables
*        o  unmeasured variables
*     determine number of degrees of freedom    
*     transform steps for transformed variables
*     transform covariance matrix for transformed variables 
*     ==================================================================
*
      IMPLICIT NONE
      INTEGER J,II,IJSYM
      DOUBLE PRECISION VII 
#include "declarefl.inc"
#include "comcfit.inc"
#include "nauxfit.inc"
      DOUBLE PRECISION X(*),VX(*),ST(*) 
*     ...
      II=0
      DO I=1,NX              ! loop on all variables
       IPAK=I 
#include "unpackfl.inc"  
       II=II+I
       VII=ABS(VX(II))       ! original diagonal element
       IF(VII.EQ.0.0D0) THEN
          NTMES=1            ! set unmeasured flag
       END IF
*      _________________________________________________________________
*      check transformed variables
       IF(NTVAR.EQ.4) THEN   ! logarithmic transformation - check
          IF(X(I).LE.0.0) THEN
             IF(IPR.GE.2) THEN ! error condition
                WRITE(LUNSIM,*)
     +          'Variable',I,' is',X(I),' reset to normal'
             END IF 
             NTVAR=0        ! reset: X(i) has to be positive
          END IF
       ELSE IF(NTVAR.EQ.5) THEN ! sqrt transformation - check   
          IF(X(I).LE.0.0) THEN
             IF(IPR.GE.2) THEN ! error condition
                WRITE(LUNSIM,*)
     +          'Variable',I,' is',X(I),' reset to normal'
             END IF
             NTVAR=0        ! reset: X(i) has to be positive
          END IF
       END IF
*      _________________________________________________________________
*      define step size for derivative calculation
       IF(VII.NE.0.0D0) THEN         ! measured variable
          IF(ST(I).GT.0.0D0) THEN
             ST(I)=MIN(ST(I),DERFAC*SQRT(VII)) ! user step, if smaller 
          ELSE IF(ST(I).LE.0.0D0) THEN         ! step is undefined
             ST(I)=DERFAC*SQRT(VII)            ! step from cov matrix
             IF(NTVAR.NE.2) THEN
                ST(I)=MIN(ST(I),DERLOW*MAX(1.0D-6,ABS(X(I))))
             ELSE   ! Poisson
                ST(I)=MIN(ST(I),DERLOW*MAX(1.0D-6,ABS(X(I)+1.0D0)))
             END IF 
          END IF
       ELSE IF(VII.EQ.0.0D0) THEN    ! unmeasured variable
          NDF=NDF-1                  ! reduce degrees of freedom
          DO J=1,NX
           VX(IJSYM(I,J))=0.0D0      ! clear matrix elements 
          END DO
          ST(I)=DERUFC*MAX(1.0D0,ABS(X(I)))
       END IF
       IF(NTINE.EQ.1) ST(I)=0.0D0   ! fix
*      _________________________________________________________________
*      transform steps for transformed variables
       IF(ST(I).EQ.0.0D0) THEN      !
          NTINE=1                   ! fixed by user
       ELSE                         !
          IF(NTVAR.EQ.3) THEN       ! lognormal variable
             WRITE(*,*) 'Step from to ',ST(I),ST(I)/X(I)
             ST(I)=ST(I)/X(I)       ! change step to log step
          ELSE IF(NTVAR.EQ.4) THEN  ! sqrt variable
             ST(I)=0.5D0*ST(I)/SQRT(X(I)) ! change step to sqrt step 
          END IF 
       END IF
C       CALL ATETOI  ! transform external to internal variables

*     __________________________________________________________________
*     transform covariance matrix for transformed variables
       IF(NTVAR.EQ.4) THEN ! transform covariance matrix for logn
          WRITE(*,*) 'lognormal variable',I,X(I)
          DO J=1,NX
           VX(IJSYM(I,J))=VX(IJSYM(I,J))/X(I)
           IF(I.EQ.J) VX(IJSYM(I,J))=VX(IJSYM(I,J))/X(I)
          END DO
       ELSE IF(NTVAR.EQ.5) THEN ! ... and for sqrt
          WRITE(*,*) 'sqrt variable',I,X(I)
          DO J=1,NX
           VX(IJSYM(I,J))=VX(IJSYM(I,J))*0.5D0/SQRT(X(I))
           IF(I.EQ.J) VX(IJSYM(I,J))=VX(IJSYM(I,J))*0.5D0/SQRT(X(I))
          END DO
       END IF
#include "packfl.inc"
      END DO
      END


      SUBROUTINE ANUMDE(X,F, A,ST,XL,FC,HH, JRET)! numerical derivatives
*     ==================================================================
*     calculation of numerical derivatives, one variable at a time
*        check limits for variable
*        calculate displaced value of variable
*        calculate derivative in Jacobian matrix A 
*        classify derivative properties of variable
*        return with JRET=-1, unless finished
*
*     X(.),F(.) = variables and constraints
*     A(.)      = Jacobian 
*     ST(.)     = steps
*     XL(2,.)   = limits
*     FC(.)     = Constraints at central variable values
*     HH(.)     = Constraints at first step (variable + step)
*     __________________________________________________________________ 
*     Logic:
*     Entry  1: start with I=0 ...
*            2: continue with indices
*     Return 1: continue with constraint evaluation
*            2: Jacobian ready  
*     ==================================================================
*      
      IMPLICIT NONE
      INTEGER JRET,ILR,IJ,J
      INTEGER NZER,NONZ,NALZ 
#include "comcfit.inc"
#include "nauxfit.inc"
#include "declarefl.inc"
      DOUBLE PRECISION X(*),F(*),A(*),ST(*),XL(2,*),FC(*),HH(*)
      DOUBLE PRECISION XD(2),XT(2),XSAVE,DER,STM
      DOUBLE PRECISION RATDIF,RATMAX,DERZER
      LOGICAL LIMDEF,TINUE
      DATA TINUE/.FALSE./       ! entry flag
*     ...
c      WRITE(*,*) 'NUMDE entered',TINUE
      JRET=-1                   ! ...means continue at return
      IF(TINUE) GOTO 30         ! continue
      TINUE=.TRUE.
      I=0                       ! initialize derivative loop
 10   IF(I.GE.NX) THEN          ! finished
         JRET=0                 ! ... means differentation finished
         TINUE=.FALSE.          ! Jacobian ready
         RETURN
      END IF
      I=I+1                     ! next variable
c      WRITE(*,*) 'I ST(I)',I,ST(I)
      IF(ST(I).EQ.0.0) GOTO 10  ! skip fixed variable
      IPAK=I        
#include "unpackfl.inc"
c      WRITE(*,*) 'I NTDER',I,NTDER 
      IF(NTDER.GE.4) GOTO 10    ! skip repeated derivative calculation
      XSAVE=X(I)                ! save current value of variable
      ILR=0                     ! define steps
*     __________________________________________________________________
*     check limits for variable
      LIMDEF=XL(1,I).NE.XL(2,I) ! true if limits defined
      IF(LIMDEF) THEN
         IF(XSAVE+ST(I).GT.XL(2,I).OR.XSAVE-ST(I).GT.XL(1,I)) THEN
            STM=0.9999*MIN(XL(2,I)-XSAVE,XSAVE-XL(1,I)) ! minimal step size
            IF(3.0*STM.GT.ST(I)) THEN
               ST(I)=STM           ! use smaller symmetric step
            ELSE 
               STM=0.4999*MAX(XL(2,I)-XSAVE,XSAVE-XL(1,I)) ! minimal step size 
               IF(2.0*STM.LT.ST(I)) ST(I)=STM 
               IF(ST(I).LT.XL(2,I)-XSAVE) THEN
                  XD(1)=XSAVE+ST(I)
                  XD(2)=XSAVE+ST(I)*2.0 ! + one-sided steps
                  ILR=1
               ELSE
                  XD(1)=XSAVE-ST(I)
                  XD(2)=XSAVE-ST(I)*2.0 ! - one-sided steps
                  ILR=2  
               END IF
            END IF
         END IF
      END IF 
*     __________________________________________________________________
*     define displaced values for derivative calculation
      IF(ILR.EQ.0) THEN
 20      IF(NTVAR.EQ.0.OR.NTVAR.EQ.2.OR.NTVAR.EQ.3) THEN
            XT(1)=XSAVE+ST(I)   ! symmetric (two-sided) steps
            XT(2)=XSAVE-ST(I)
            XD(1)=XT(1)
            XD(2)=XT(2)
         ELSE IF(NTVAR.EQ.1) THEN   ! 1/x
            XT(1)=1.0D0/XSAVE+ST(I) ! internal
            XT(2)=1.0D0/XSAVE-ST(I) 
            XD(1)=1.0D0/XT(1)       ! external
            XD(2)=1.0D0/XT(2)
         ELSE IF(NTVAR.EQ.4) THEN   ! log-normal
            XT(1)=LOG(XSAVE)+ST(I)  ! internal
            XT(2)=LOG(XSAVE)-ST(I)
            XD(1)=EXP(XT(1))        ! external
            XD(2)=EXP(XT(2))
         ELSE IF(NTVAR.EQ.5) THEN ! sqrt
            IF(ST(I)**2.GE.XSAVE) THEN
               IF(XSAVE.LE.0.0D0) THEN
                  IF(IPR.GT.1) WRITE(LUNSIM,*)
     +            'Variable',I,' is',X(I),' reset to normal'
                  NTVAR=0 
#include "packfl.inc"
                  GOTO 20
               ELSE
                  ST(I)=0.9D0*SQRT(XSAVE)
               END IF
            END IF
            XT(1)=SQRT(XSAVE)+ST(I) ! internal 
            XT(2)=SQRT(XSAVE)-ST(I)
            XD(1)=XT(1)**2          ! external
            XD(2)=XT(2)**2
         ELSE IF(NTVAR.EQ.6) THEN ! x**power       
            XT(1)=XSAVE**XL(2,I)+ST(I)   ! internal
            XT(2)=XSAVE**XL(2,I)-ST(I)
            XD(1)=XT(1)**(1.0D0/XL(2,I)) ! external
            XD(2)=XT(2)**(1.0D0/XL(2,I))
         END IF
      END IF
*     __________________________________________________________________
*     set variable to displaced value and return for calculation
      X(I)=XD(1)                ! first step
      I=-I
      RETURN 
*     __________________________________________________________________
*     continue 
 30   IF(I.LT.0) THEN           ! calculation of first step done ... 
         DO J=1,NF
          HH(J)=F(J)            ! save constraint values
         END DO
         I=-I                   ! reverse flag
         X(I)=XD(2)             ! set next step ...
         RETURN                 ! ... and return for second step
      END IF
*     __________________________________________________________________
*     INIT ne 0: second step done - calculate derivative
      X(I)=XSAVE                ! restore variable I
      IJ=I                      ! derivative calculation
      NZER=0                    ! reset: number of zero derivatives
      NONZ=0                    ! reset: number of non-zero derivatives
      NALZ=0                    ! flag all derivatives are zero
      RATMAX=0.0D0              ! max of diff-ratio
      DERZER=0.0D0              ! abs of nonzero-derivative
      DO J=1,NF                 ! loop on all constraint functions
       IF(ILR.EQ.0) THEN        ! symmetric formula
c         DER=0.5D0*(HH(J)-F(J))/ST(I) ! numerical 1. derivative
          DER=(HH(J)-F(J))/(XT(1)-XT(2))  !!! internal variable
c          WRITE(*,246) I,J,HH(J),F(J),XT(1),XT(2),DER
c 246      FORMAT('der I J',2I3,5F10.4) 
       ELSE                     ! asymmetric formula
          DER=0.5D0*(3.0D0*FC(J)+F(J)-4.0D0*HH(J))/ST(I)
          IF(ILR.EQ.2) DER=-DER ! sign
       END IF
*      _________________________________________________________________
*      classify derivative properties of variable I
       IF(A(IJ).NE.0.0D0.OR.DER.NE.0.0D0) NALZ=1 ! non all zero
       IF(DER.EQ.0.0D0) THEN
          NZER=NZER+1           ! derivative zero - count 
       ELSE                     ! derivative non-zero
          RATDIF=ABS(A(IJ)-DER)/(ABS(A(IJ))+ABS(DER))
          RATMAX=MAX(RATMAX,RATDIF) 
          DERZER=ABS(DER)       ! abs value of this derivative
          NONZ=NONZ+1           ! count non-zero derivative
       END IF
       A(IJ)=DER                ! insert into Jacobian matrix A
       IJ=IJ+NX
      END DO
      IF(NTDER.EQ.0) THEN
         NTDER=1
      ELSE
         IF(NONZ.EQ.1.AND.ABS(DERZER-1.0D0).LT.1.0D-12) THEN
            NTDER=MIN(NTDER+1,7)
         ELSE IF(RATMAX.LT.1.0D-12) THEN
            NTDER=MIN(NTDER+1,7)
         ELSE
            NTDER=2             ! reset to 2
         END IF
      END IF 
      NTDER=0 
#include "packfl.inc"
      GOTO 10 
      END

      SUBROUTINE ANITER(X,VX,F,A,XP,RH,WM,DX) ! next iteration step
      IMPLICIT NONE
      DOUBLE PRECISION X(*),VX(*),F(*),A(*),XP(*),RH(*),WM(*),DX(*)
#include "comcfit.inc"
#include "nauxfit.inc"
#include "declarefl.inc"
      INTEGER IA,II,J,NRANK 
      DOUBLE PRECISION DIAG(1000),QNEXT(1000),SCALXY
*     ... 
      ITER=ITER+1                 ! start next iteration
      CHSQP=CHISQ                 ! save current chi^2

c      A(2)=1.0D0
c      A(3)=-1.0D0
c      WRITE(*,321) 'Jacobian',(A(I),I=1,NX*NF) 
c321  FORMAT(A/(4F15.6))
*     __________________________________________________________________
*     right-hand side of equation
      NCST=0
      DO  I=1,NX                  ! first NX components
       RH(I)=0.0D0                 ! define right hand side of equation
      END DO
      DO J=1,NF                   ! next NF components
       AUX(INDFC+J)=F(J)
       RH(NX+J)=-F(J)
      END DO
      IA=0
      DO J=1,NF                   ! "subtract" actual step 
       RH(NX+J)=RH(NX+J)+SCALXY(A(IA+1),DX,NX)
       IA=IA+NX
       AUX(INDHH+J)=RH(NX+J)       ! right hand side for chi**2
      END DO
c      WRITE(*,*) 'RHS',(RH(I),I=1,NX+NF) 
*     __________________________________________________________________
*     form matrix and solve
      DO I=1,(NX*NX+NX)/2
       WM(I)=-VX(I)                ! copy -VX(.) into W_11
      END DO
      II=0                        ! modify V for Poisson variables
      DO I=1,NX
       II=II+I
       IPAK=I
#include "unpackfl.inc"
c       WRITE(*,*) 'NTVAR=2 II',II,X(I),NTVAR
       IF(NTVAR.EQ.2) THEN ! Poisson
          WM(II)=-SQRT(1.0+X(I)**2)! -MAX(ABS(X(I)),1.0D0)
       END IF
      END DO
      CALL DUMINV(A, WM,RH,NX,NF, 1, NRANK, DIAG,QNEXT)
c      WRITE(*,*) NRANK,' is rank. NX NF ',NX,NF 
      CHISQ=-SCALXY(AUX(INDHH+1),RH(NX+1),NF)  ! next chi^2
      IF(CHISQ.LT.0.0D0) CHISQ=0.0D0
*     __________________________________________________________________
*     handle corrections and cutstep
      WEIGHT=1.0D0                ! default weight
      IF(ITER.GT.1.AND.CHISQ.GE.2.00D0*CHSQP) WEIGHT=0.1D0
      IF(ITER.GT.1.AND.CHISQ.GE.3.00D0*CHSQP) WEIGHT=0.05D0    

      DO I=1,NX
       XP(I)=DX(I)                ! save previous corrections       
       DX(I)=RH(I)                 ! store new corrections
      END DO
c      WRITE(*,*) 'Prev corrections',(XP(I),I=1,NX)
c      WRITE(*,*) 'Curr corrections',(DX(I),I=1,NX)
c      WRITE(*,*) 'Weight=',WEIGHT
      END 


      SUBROUTINE ADDTOX(X,XS,DX,XP)
      DOUBLE PRECISION X(*),XS(*),DX(*),XP(*)
#include "comcfit.inc"
#include "nauxfit.inc"
#include "declarefl.inc"
*     ...
      DO I=1,NX
       DX(I)=WEIGHT*DX(I)+(1.0D0-WEIGHT)*XP(I)   ! reduce step evtl.
       IPAK=I
#include "unpackfl.inc"
       IF(NTVAR.EQ.0.OR.NTVAR.EQ.2.OR.NTVAR.EQ.3) THEN 
          X(I)=XS(I)+DX(I) ! correct x and return to test constraints
       ELSE IF(NTVAR.EQ.4) THEN ! log-normal
          X(I)=EXP(LOG(XS(I))+DX(I))  
       ELSE IF(NTVAR.EQ.5) THEN ! sqrt
          X(I)=(SQRT(XS(I))+DX(I))**2
       END IF
      END DO
      END


      SUBROUTINE ANTEST(IRET)      ! test convergence
#include "comcfit.inc"
      DOUBLE PRECISION CM(14)
*     ....
      IRET=-1                      ! calculate new Jacobian  
      IUNPH=0   
*     __________________________________________________________________
*     combined measure?
c      WRITE(*,*) 'ITER,NCST',ITER,NCST,CHISQ,FRMS,CHSQP,FRMSP
      IF(ITER.EQ.1.AND.NCST.EQ.0) THEN
         FACTOR=CHISQ/(ABS(FTESTP-FTEST)+EPSF)
         DO I=1,14 
          CM(I)=0.0D0 
         END DO
         CM(9)=0.5D0                  ! damping factor 
         CM(1)=1.0D0          +1.0D0
         CM(2)=FRMS**2        +FRMSP**2
         CM(3)=FRMS**4        +FRMSP**4 
         CM(4)=CHISQ          +CHSQP
         CM(5)=CHISQ*FRMS**2  +CHSQP*FRMSP**2
         CALL LESFCM(CM)            ! fit
         CM(13)=CM(6)-CM(7)*FRMS**2 ! combined penalty
      ELSE IF(NCST.EQ.0) THEN    
         CM(1)=CM(9)*CM(1)+1.0D0
         CM(2)=CM(9)*CM(2)+FRMS**2
         CM(3)=CM(9)*CM(3)+FRMS**4
         CM(4)=CM(9)*CM(4)+CHISQ
         CM(5)=CM(9)*CM(5)+CHISQ*FRMS**2
         CALL LESFCM(CM)            ! fit
         CM(14)=CM(13) 
         CM(13)=CM(6)-CM(7)*FRMS**2 ! combined penalty
c         WRITE(*,*) 'Penalty',CM(6),CM(7),CM(13)
      END IF   
c      WRITE(*,*) 'ANTEST:',CHISQ,FTEST,CHISQ+FACTOR*FTEST
      PENALT=CM(13)
      CALL AIPRIN(X,VX,3,IRET)
*     __________________________________________________________________ 
*     cutstep
c      WRITE(*,*) 'Test cutstep:',FTEST,FTESTP,EPSF,2.0D0*FTESTP+EPSF
      IF(NCST.LT.2.AND.
     +  (IUNPH.NE.0.OR.(ITER.GT.1.AND.FTEST.GT.2.0D0*FTESTP+EPSF))) THEN
         NCST=NCST+1
         WEIGHT=0.25D0
         WEIGHT=0.50D0  
c         IF(FTEST/FTESTP.GT. 5.0D0) WEIGHT=0.10D0
c         IF(FTEST/FTESTP.GT.10.0D0) WEIGHT=0.05D0
         IRET=-2                  ! cutstep - add corrections 
         RETURN                   
      END IF   

c      IF(NCST.EQ.0.AND.PENALT.LT.0.0D0) THEN
c         WEIGHT=0.1D0
c         NCST=NCST+1
c         IRET=-2                  ! cutstep - add corrections
c         WRITE(*,*) 'Penalty!'
c         RETURN
c      END IF 

c      IF(NCST.EQ.0.AND.CM(13).GT.CM(14)) THEN
c         WEIGHT=0.1D0
c         NCST=NCST+1 
c         IRET=-2                  ! cutstep - add corrections
c         WRITE(*,*) 'Penalty increasing!'
c         RETURN
c      END IF

c      IF(ITER.EQ.1.AND.NCST.EQ.0) THEN
c         WEIGHT=0.1 
c         NCST=NCST+1
c         IRET=-2                  ! cutstep - add corrections
c         RETURN
c      END IF 

*     __________________________________________________________________
*     convergent
      IF(ITER.GE.2.AND.NCST.EQ.0) THEN
         DCHISQ=CHISQ-CHSQP
         IF(ABS(DCHISQ).LE.EPSCHI.AND.FTEST.LT.EPSF) THEN
            IRET=0                ! convergence

C            WRITE(*,*) '>>>> convergence'
 
            RETURN
         END IF 
      END IF                                                  
*     __________________________________________________________________
*     failure
      IF(ITER.GT.ITERMX) IRET=2   ! non-convergence
      END

      SUBROUTINE ACOPXV(X,VX,DX,AS,WM,PU)
      DOUBLE PRECISION X(*),VX(*),DX(*),AS(*),WM(*),PU(*),SCOPY
#include "comcfit.inc"
#include "nauxfit.inc"
*     ...
*     __________________________________________________________________
*     convergence: pull calculation
      II=0           
      DO I=1,NX
       AS(I)=X(I)
       II=II+I
       PU(I)=0.0  
       IF(VX(II).GT.0.0) THEN
          IF(VX(II)-WM(II).GT.0.0) PU(I)=DX(I)/SQRT(VX(II)-WM(II))
       END IF
      END DO
*     __________________________________________________________________
*     copy/exchange result/input covariance matrix
      DO I=1,(NX*NX+NX)/2
        SCOPY=VX(I)
        VX(I)=WM(I)              ! copy fitted covariance matrix
        WM(I)=SCOPY              ! ... and save input matrix
      END DO
*     __________________________________________________________________
*     apply transformation to covariance matrix
      CALL ATITOE(X,VX)
      END 

      SUBROUTINE ATETOI(X,VX)  ! transform external to internal variables
      IMPLICIT NONE
      DOUBLE PRECISION X(*),VX(*)
      INTEGER J,IJSYM 
#include "comcfit.inc"
#include "nauxfit.inc"
#include "declarefl.inc"
*     ...
      DO I=1,NX                 ! transformation back
       IPAK=I
#include "unpackfl.inc"
*      transform covariance matrix for transformed variables
       IF(NTVAR.EQ.4) THEN ! transform covariance matrix for logn
          WRITE(*,*) 'lognormal variable',I,X(I)
          DO J=1,NX
           VX(IJSYM(I,J))=VX(IJSYM(I,J))/X(I)
           IF(I.EQ.J) VX(IJSYM(I,J))=VX(IJSYM(I,J))/X(I)
          END DO 
       ELSE IF(NTVAR.EQ.5) THEN ! ... and for sqrt
          WRITE(*,*) 'sqrt variable',I,X(I)
          DO J=1,NX
           VX(IJSYM(I,J))=VX(IJSYM(I,J))*0.5D0/SQRT(X(I))
           IF(I.EQ.J) VX(IJSYM(I,J))=VX(IJSYM(I,J))*0.5D0/SQRT(X(I))
          END DO
       END IF
      END DO
      RETURN                                            

      ENTRY      ATITOE(X,VX)  ! transform internal to external variables
      DO I=1,NX                 ! transformation back
       IPAK=I
#include "unpackfl.inc"
       IF(NTVAR.EQ.4) THEN
          WRITE(*,*) 'cov matrix back with ',X(I)
          DO J=1,NX ! log-normal
           VX(IJSYM(I,J))=VX(IJSYM(I,J))*X(I)
           IF(I.EQ.J) VX(IJSYM(I,J))=VX(IJSYM(I,J))*X(I)
          END DO
       ELSE IF(NTVAR.EQ.5) THEN
          DO J=1,NX ! sqrt
           VX(IJSYM(I,J))=VX(IJSYM(I,J))*2.0D0*SQRT(X(I))
           IF(I.EQ.J) VX(IJSYM(I,J))=VX(IJSYM(I,J))*2.0D0*SQRT(X(I))
          END DO
       END IF
      END DO
      END

      SUBROUTINE AIPRIN(X,VX,IARG,IRET)
      DOUBLE PRECISION X(*),VX(*),CHP,CHPROB,DINGAU
      INTEGER IARG,J
C      INTEGER NDTOTL  
#include "comcfit.inc"
#include "nauxfit.inc"
#include "cprofil.inc"
      CHARACTER*19 TEXT(5)
      DATA         TEXT   /'Chisquare too high ',
     +                     'Too many iterations',
     +                     'Unphysical region  ',
     +                     'NDF less or equal 0',
     +                     'AUX dimension small'/
*     ...
C      NDTOTL=120
      IF(IARG.EQ.1) GOTO 10
      IF(IARG.EQ.2) GOTO 20
      IF(IARG.EQ.3) GOTO 30
      IF(IARG.EQ.4) GOTO 40
      IF(IARG.EQ.5) GOTO 50
*     __________________________________________________________________
*     APLCON Titel print
      IF(IPR.GE.3) THEN  ! print
         WRITE(LUNSIM,*) ' '
         WRITE(LUNSIM,*) 'APLCON - constrained least squares',
     +   '              Version  01/10/2011'
         WRITE(LUNSIM,*) ' '
         IF(IPR.GE.5.AND.NX.LE.128) THEN
            CALL CFPRV(LUNSIM,X,VX,NX)
         END IF
      END IF
      RETURN
*     __________________________________________________________________      
*     APLCON memory space error
 10   IF(NDTOTL.GT.NAUX) THEN
         WRITE(*,*) ' '
         WRITE(*,*) 'APLCON  - constrained least squares'
         WRITE(*,*) '          Case ',NCASE
         WRITE(*,*) 'Insufficient space in internal array AUX(',NAUX,')'
         WRITE(*,*) 'Required:',NDTOTL,' elements (at least)'
         WRITE(*,*) '-> stop'
         WRITE(*,*) ' '
         STOP
      END IF
      RETURN
*     __________________________________________________________________
*     initial printout normal start
 20   IF(IPR.GE.3) THEN
         WRITE(LUNSIM,*) ' '
         WRITE(LUNSIM,101) ' Constrained least squares fit:'
         WRITE(LUNSIM,101) '                         case',NCASE 
         WRITE(LUNSIM,101) '              nr of variables',NX
         WRITE(LUNSIM,101) '            nr of constraints',NF
         WRITE(LUNSIM,101) '           degrees of freedom',NDF     
         WRITE(LUNSIM,102) '                      epsilon',EPSF
         WRITE(LUNSIM,102) '     factors for numer. diff.',DERFAC,
     +                     DERUFC,DERLOW
         WRITE(LUNSIM,101) '          used array elements',NDTOT
         WRITE(LUNSIM,101) '         total array elements',NDPDIM     
      END IF

      IF(IPR.GE.6.AND.NX.LE.32) THEN ! print
         CALL CFPRV(LUNSIM,X,VX,NX)  ! initial parameter values
      END IF
 101  FORMAT(3X,A30,I7)
 102  FORMAT(3X,A30,5X,1P,3E8.1)
 109  FORMAT(3X,5G12.5)
 105  FORMAT(' No convergence (',I3,' ITER, IRET =',I2,')  ',A)
 119  FORMAT('Convergence after',I7,' iterations with chi^2=',G15.8,
     +       '  ndf=',I4)
 120  FORMAT(22X,'corresponding to P-value=',F11.6)
 121  FORMAT(16X,'1-,2-sided standard deviations=',F11.6,',',F13.6)
 122  FORMAT('P-value converted to number of standard deviations:')
 110  FORMAT(/'Iteration    calls        chi^2  Frms       |F| ',
     +     6X,'dchi^2',2X,'chi^2&Frms')
 111  FORMAT(1X,I5,2X,I10,13X,E10.2,E10.2)
      RETURN
*     __________________________________________________________________
*     iteration printout
 30   IF(IPR.GE.4.AND.ITER.EQ.0) THEN 
         WRITE(LUNSIM,110)
         WRITE(LUNSIM,111) ITER,NCALLS,FRMS,FTEST
      END IF 
      IF(IPR.GE.4.AND.ITER.GE.1) THEN  
         IF(NCST.EQ.0) THEN
            WRITE(LUNSIM,211) ITER,NCALLS,CHISQ,FRMS,FTEST,
     +      CHISQ-CHSQP,PENALT
         ELSE
            WRITE(LUNSIM,212) ITER,NCST,NCALLS,FRMS,FTEST 
         END IF 
      END IF
      IF(IPR.GE.7.AND.NCST.EQ.0) THEN
         CALL CFGMPR(LUNSIM,X,NX,1,'of X-values')
         CALL CFGMPR(LUNSIM,F,NF,1,' of constraint function values')
      END IF
  211 FORMAT(1X,I5,    2X,I10,F13.3,E10.2,E10.2,E11.2,F10.3) 
  212 FORMAT(1X,I5,'.',I1,I10,13X,E10.2,E10.2)
      RETURN
*     __________________________________________________________________
*     final return with IRET = 0,1,2,3,4,5
 40   IF(IPR.GE.1) THEN      ! print 
C         WRITE(*,*) 'Printing result IPR, IRET=',IPR,IRET
C         IF(IRET.GT.5) IRET=0 !!!!!!!!!!!!!!!
         IF(IRET.EQ.0) THEN             ! convergence (IRET=0)
C            WRITE(*,*) 'Printing result IPR, IRET=',IPR,IRET
            IF(IPR.GE.3) THEN
C               WRITE(*,*) 'before print ITER ... IPR=',IPR
               WRITE(LUNSIM,*) ' '
               WRITE(LUNSIM,119) ITER,CHISQ,NDF
               IF(NDF.GT.0) THEN
                  CHP=CHPROB(CHISQ,NDF)   ! p-value  
                  WRITE(LUNSIM,120) CHP
                  IF(CHP.LT.0.05D0) THEN 
                     WRITE(LUNSIM,122)
                     WRITE(LUNSIM,121) DINGAU(1.0D0-0.5D0*CHP),
     +                                 DINGAU(1.0D0-CHP)
                  END IF
               END IF
            END IF     
            IF(IPR.GE.4) THEN
               IF(NX.LE.128) THEN
                  CALL CFPRVP(LUNSIM,X,VX,AUX(INDPU+1),NX)
                  IF(IPR.GE.5) THEN
                     WRITE(LUNSIM,*) ' '
                     CALL CFCORR(LUNSIM,VX,NX)
                  END IF 
               ELSE 
                  WRITE(LUNSIM,*) '  x-vector (fitted):'
                  WRITE(LUNSIM,109) (X(I),I=1,NX)
               END IF
               WRITE(LUNSIM,*) ' '
            END IF              
         ELSE IF(IPR.GE.3) THEN         ! non-convergence (IRET > 0)
            WRITE(LUNSIM,105) ITER,IRET,TEXT(IRET)
            WRITE(LUNSIM,*) ' ' 
         END IF
      END IF 
      RETURN
*     __________________________________________________________________
*     
 50   CONTINUE
      RETURN
*     __________________________________________________________________
*     
      ENTRY APRINI(IARG)
      IF(IPR.GE.3) THEN  ! print
         WRITE(LUNSIM,*) ('_',J=1,71)
         WRITE(LUNSIM,*) ' '
         WRITE(LUNSIM,*) 'APLCON - constrained least squares',
     +   '                  Version  01/10/2011'
         WRITE(LUNSIM,*) ('_',J=1,71)
         WRITE(LUNSIM,*) ' '
      END IF
      END 

      SUBROUTINE LESFCM(C)
      DOUBLE PRECISION C(13)
*     ...
      C(8)=C(1)*C(3)-C(2)*C(2)   ! determinant
      C(10)= C(3)/C(8)           ! V_11
      C(11)=-C(2)/C(8)           ! V_12
      C(12)= C(1)/C(8)           ! V_22
      C(6)=C(4)*C(10)+C(5)*C(11) ! 1. parameter
      C(7)=C(4)*C(11)+C(5)*C(12) ! 2. parameter
C      WRITE(*,*) 'Solution ',C(6),C(7)
      END  

      SUBROUTINE CHNDPV(CHI2,ND,PVAL)
#include "comcfit.inc"
      DOUBLE PRECISION CHPROB
      REAL CHI2,PVAL
      INTEGER ND
*     ...
      CHI2=CHISQ               ! chi^square
      ND  =NDF                 ! number of degrees of freedom  
      PVAL=CHPROB(CHISQ,NDF)   ! p-value
      END  


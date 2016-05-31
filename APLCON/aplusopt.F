
      SUBROUTINE ADUMMY
*     __________________________________________________________________
      IMPLICIT NONE
      DOUBLE PRECISION STEP,XLOW,XHIG,POW,ARG
      INTEGER IT,LUNP,JPR,NBINOM,I1,I2
#include "comcfit.inc"
#include "nauxfit.inc"
#include "cprofil.inc"   
#include "declarefl.inc"
*     __________________________________________________________________
*
*     parameters for fit method
*     __________________________________________________________________
      ENTRY APRINT(LUNP,JPR)
      LUNSIM=LUNP       ! print unit
      IF(LUNSIM.LE.0)     LUNSIM=6
      IPR=JPR           ! print flag
      RETURN
*     __________________________________________________________________
      ENTRY APDEPS(ARG)                 ! constraint accuracy
      EPSF  =ARG        ! |F| accuracy
      RETURN
*     __________________________________________________________________
      ENTRY APEPSCHI(ARG)               ! chi2 accuracy
      EPSCHI=ARG        
      RETURN      
*     __________________________________________________________________
      ENTRY APDERF(ARG)                 ! factor for step definition
      DERFAC=ARG
      RETURN
*     __________________________________________________________________
      ENTRY APDERU(ARG)                 ! factor for step definition
      DERUFC=ARG
      RETURN
*     __________________________________________________________________
      ENTRY APDLOW(ARG)                 ! factor for step definition
      DERLOW=ARG
      RETURN
*     __________________________________________________________________
      ENTRY APITER(IT)                  ! iteration limit
      ITERMX=MAX(3,IT)  ! max number of iterations
      RETURN
*     __________________________________________________________________
*
*     measured/unmeasured variables
*     __________________________________________________________________
*
      ENTRY APROFL(I1,I2)
      IF(NSECA.GE.MSECA) RETURN
      IF(I1.LT.1.OR.I1.GT.NX) RETURN
      IF(I2.LT.0.OR.I2.GT.NX) RETURN
      NSECA=NSECA+1
      NPSEC(1,NSECA)=I1           ! index for profile analysis
      NPSEC(2,NSECA)=I2
      IF(I2.EQ.0) THEN   ! 1-parameter profil
         NFADD=MAX(1,NFADD)
         IPAK=I1
#include "unpackfl.inc"           
         IF(MOD(NTPRF,2).EQ.0) NTPRF=NTPRF+1
#include "packfl.inc"
      ELSE               ! 2-parameter profile
         NFADD=2
         IPAK=I1
#include "unpackfl.inc"
         IF(NTPRF.LE.1) NTPRF=NTPRF+2
#include "packfl.inc"
         IPAK=I2
#include "unpackfl.inc"
         IF(NTPRF.LE.1) NTPRF=NTPRF+2
#include "packfl.inc"
      END IF
      RETURN
*     __________________________________________________________________      
      ENTRY APSTEP(I,STEP)              ! step size for numdif
      IF(I.LT.1.OR.I.GT.NX) RETURN 
      IPAK=I
#include "unpackfl.inc"     
      AUX(INDST+I)=ABS(STEP)            ! ST(I)= ...
      IF(STEP.NE.0.0D0) THEN
         NTINE=0 ! variable
      ELSE
         NTINE=1 ! fixed by user
      END IF
      GOTO 100
*     __________________________________________________________________
      ENTRY APFIX(I)                    ! fixed parameter
      IF(I.LT.1.OR.I.GT.NX) RETURN 
      IPAK=I
#include "unpackfl.inc"
      NTINE=1    ! fixed by user 
      GOTO 100
*     __________________________________________________________________
      ENTRY APLIMT(I,XLOW,XHIG)        ! range of variable
      IF(I.LT.1.OR.I.GT.NX) RETURN
      IPAK=I
#include "unpackfl.inc"
      AUX(INDLM+2*(I-1)+1)=MIN(XLOW,XHIG) ! lower limit XL(1,I)
      AUX(INDLM+2*(I-1)+2)=MAX(XLOW,XHIG) ! upper limit XL(2,I)
      NTLIM=4 
      GOTO 100
*     __________________________________________________________________
      ENTRY APTRIN(I)                   ! inverse value
      IF(I.LT.1.OR.I.GT.NX) RETURN
      IPAK=I
#include "unpackfl.inc"
      NTVAR=1    ! transformation to inverse
      GOTO 100
*     __________________________________________________________________
      ENTRY APOISS(I)                   ! Poisson distributed variable
      IF(I.LT.1.OR.I.GT.NX) RETURN
      IPAK=I
#include "unpackfl.inc"
      NTVAR=2    ! Poisson distributed variable
      NTLIM=1    
c      WRITE(*,*) 'APOISS I,IPAK,NTVAR,NTLIM ',I,IPAK,NTVAR,NTLIM    
      GOTO 100
*     __________________________________________________________________
      ENTRY ABINOM(I,NBINOM)           ! Binomial distributed variable
      IF(I.LT.1.OR.I.GT.NX) RETURN
      IPAK=I
#include "unpackfl.inc"
      NTVAR=3    ! Binomial distributed variable
      AUX(INDLM+2*I)=NBINOM
      NTLIM=1  
      GOTO 100
*     __________________________________________________________________
      ENTRY APLOGN(I)                  ! Lognormal distributed variable
      IF(I.LT.1.OR.I.GT.NX) RETURN
      IPAK=I
#include "unpackfl.inc"
      NTVAR=4    ! Lognormal distributed variable
      NTLIM=1
      GOTO 100
*     __________________________________________________________________
      ENTRY APSQRT(I)                  ! SQRT transformation
      IF(I.LT.1.OR.I.GT.NX) RETURN
      IPAK=I
#include "unpackfl.inc"
      NTVAR=5    ! SQRT transformation
      NTLIM=1 
      GOTO 100
*     __________________________________________________________________
      ENTRY APOWER(I,POW)              ! x^power transformation
      IF(I.LT.1.OR.I.GT.NX) RETURN
      IPAK=I
#include "unpackfl.inc"
      NTVAR=6    ! x^power transformation
      AUX(INDLM+2*I)=POW
      NTLIM=1 
      GOTO 100
*     __________________________________________________________________
      ENTRY APOSIT(I)                  ! positive 
      IF(I.LT.1.OR.I.GT.NX) RETURN
      IPAK=I
#include "unpackfl.inc"
      NTLIM=1    ! positive
c      WRITE(*,*) 'APOSIT I,IPAK, NTLIM=',I,IPAK,NTLIM
      GOTO 100
*     __________________________________________________________________
 100  CONTINUE
#include "packfl.inc"
C#include "unpackfl.inc"
C      WRITE(*,*) 'NTVAR =',NTVAR,I
      END 

      SUBROUTINE APSTAT(FOPT,NFUN,NITER)          ! return Fopt and Nfun      
*     __________________________________________________________________
*     return information after the fit
*     __________________________________________________________________
      IMPLICIT NONE
#include "comcfit.inc"
      INTEGER NFUN,NITER
      DOUBLE PRECISION FOPT
*     ...
      FOPT =CHISQ
      NFUN =NCALLS
      NITER=ITER 
      END 


      SUBROUTINE SDEFIN(X,V,I,VALUE,ERROR,XPLAIN)
      DOUBLE PRECISION X(*),V(*),RHOCOP,RHOMAX
      CHARACTER*(*) XPLAIN
*     ...
      IF(I.LT.1.OR.I.GT.NVAR) RETURN
      X(I)=VALUE  
      V((I*I+I)/2)=ERROR*ERROR
      IF(XPLAIN.EQ.' ') RETURN
      CALL APNAME(I,XPLAIN)
      RETURN

      ENTRY SCORRC(X,V,I1,I2,RHO)
      IF(I1.LT.1.OR.I2.GT.NVAR) RETURN
      IF(I2.LT.1.OR.I2.GT.NVAR) RETURN
      RHOCOP=RHO
      RHOMAX=RHOCOP
      IF(RHO.GT.0.0D0) THEN
         RHOMAX=SQRT(V(IJSYM(I1,I1))/V(IJSYM(I2,I2)))
         IF(V(IJSYM(I1,I1)).GT.V(IJSYM(I2,I2))) RHOMAX=1.0D0/RHOMAX
         RHOMAX=MIN(RHOMAX,RHOCOP)
         IF(RHOMAX.NE.RHOCOP) THEN
            WRITE(*,101) ' Correlation coefficient between',I1,I2,
     +      ' changed from',RHOCOP,' to',RHOMAX
         END IF
      END IF
      V(IJSYM(I1,I2))=RHOMAX*SQRT(V(IJSYM(I1,I1))*V(IJSYM(I2,I2)))
      WRITE(*,*) V(IJSYM(I1,I1)),V(IJSYM(I1,I2)),V(IJSYM(I2,I2))
      RETURN

      ENTRY SCORRF(X,V,I1,I2,RHO)
      IF(I1.LT.1.OR.I2.GT.NVAR) RETURN
      IF(I2.LT.1.OR.I2.GT.NVAR) RETURN
      V(IJSYM(I1,I2))=RHO*SQRT(V(IJSYM(I1,I1))*V(IJSYM(I2,I2)))
      RETURN

      ENTRY SRESET(X,V,N)
      NVAR=N
      DO J=1,NVAR 
       X(J)=0.0D0 
      END DO
      DO J=1,(NVAR*NVAR+NVAR)/2
       V(J)=0.0D0 
      END DO

 101  FORMAT(A,I3,',',I3,A,F8.5,A,F8.5)
      END





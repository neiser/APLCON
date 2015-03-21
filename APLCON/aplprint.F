
*     print 


      SUBROUTINE APXOPT(X,VX)
#include "comcfit.inc"
#include "nauxfit.inc"
      DOUBLE PRECISION X(*),VX(*)
*     ...
      CALL CIPRV(6,X,VX,NX,AUX(INSTR+1),AUX(INDFL+1),AUX(INDLM+1))
      END 


      SUBROUTINE CIPRV(LUP,X,VX,N) ! print initial status 
*     ==================================================================
*     Print initial values of variables and properties
*     ==================================================================
*
      IMPLICIT NONE
#include "comcfit.inc"
C#include "nauxfit.inc"      
#include "cprofil.inc"
      INTEGER LUP,N,II,JS,I
      DOUBLE PRECISION X(*),VX(*),ERR
      CHARACTER*16 PANAME 
      CHARACTER*44 VATEXT               ! text for variable
*     ...
      WRITE(LUP,*)  ' '
      WRITE(LUP,101)
      WRITE(LUP,103)
      DO I=1,N
       II=(I*I+I)/2                     ! diagonal index
       ERR=0.0
       IF(VX(II).GT.0.0) ERR=SQRT(VX(II)) ! uncertainty
       CALL APNAME(-I,PANAME)           ! name of variable
       CALL APERTY(I,VATEXT,JS) 
       WRITE(LUP,102) PANAME,I,X(I),ERR,VATEXT(1:JS)
      END DO
      WRITE(LUP,*) ' '
      WRITE(LUP,104) NF,NFADD

 101  FORMAT(21X,'i     value         std.dev.   properties')
 102  FORMAT(1X,A16,I5,2G14.6,2X,A)
 103  FORMAT(18X,'----  ------------  ------------  ----------')
 104  FORMAT(I17,' +',I3,' constraint equations')
      END

      SUBROUTINE APERTY(I,TEXT,JS)
      IMPLICIT NONE
      INTEGER JS
      CHARACTER*(*) TEXT
      CHARACTER*44 VATEXT               ! text for variable
#include "comcfit.inc"
#include "nauxfit.inc"
#include "declarefl.inc"
*     ...
       VATEXT=' '
       IPAK=I
#include "unpackfl.inc"
       IF(NTMES.NE.1) THEN
          VATEXT='meas'
          JS=5
       ELSE
          VATEXT='unmeas'
          JS=7
       END IF
       IF(NTINE.EQ.1) THEN
          VATEXT(JS+1:JS+3)='fix'    
          JS=JS+4
       END IF
       IF(NTLIM.EQ.1) THEN
          VATEXT(JS+1:JS+2)='>0'
          JS=JS+3
       ELSE IF(NTLIM.EQ.2) THEN
          VATEXT(JS+1:JS+1)='>'
          JS=JS+2
       ELSE IF(NTLIM.EQ.3) THEN
          VATEXT(JS+1:JS+1)='<'
          JS=JS+2
       ELSE IF(NTLIM.EQ.4) THEN
          VATEXT(JS+1:JS+3)='< >'
          JS=JS+4
       END IF  
       IF(NTVAR.EQ.1) THEN
          VATEXT(JS+1:JS+3)='1/x'
          JS=JS+4
       ELSE IF(NTVAR.EQ.2) THEN   
          VATEXT(JS+1:JS+7)='Poisson'
          JS=JS+8
       ELSE IF(NTVAR.EQ.3) THEN
          VATEXT(JS+1:JS+8)='Binomial'
          JS=JS+9
       ELSE IF(NTVAR.EQ.4) THEN
          VATEXT(JS+1:JS+10)='log-normal'
          JS=JS+11
       ELSE IF(NTVAR.EQ.5) THEN
          VATEXT(JS+1:JS+4)='sqrt'
          JS=JS+5
       ELSE IF(NTVAR.EQ.6) THEN
          VATEXT(JS+1:JS+5)='power'   
          JS=JS+6
       END IF
       IF(NTMES.EQ.2) THEN
          VATEXT(JS+1:JS+8)='M-Cauchy'
          JS=JS+9
       ELSE IF(NTMES.EQ.3) THEN
          VATEXT(JS+1:JS+7)='M-Huber'
          JS=JS+8
       END IF
       IF(NTDER.GE.4) THEN
          VATEXT(JS+1:JS+9)='der=const'
          JS=JS+10  
       END IF 
       IF(NTPRF.EQ.1) THEN
          VATEXT(JS+1:JS+3)='[1]'
          JS=JS+4
       ELSE IF(NTPRF.EQ.2) THEN
          VATEXT(JS+1:JS+3)='[2]' 
          JS=JS+4
       ELSE IF(NTPRF.EQ.3) THEN
          VATEXT(JS+1:JS+5)='[1,2]'
          JS=JS+6
       END IF
      TEXT=VATEXT
      IF(JS.GT.LEN(TEXT)+1) TEXT(LEN(TEXT)-2:LEN(TEXT))='...'
      JS=MIN(JS,LEN(TEXT))    
      END 


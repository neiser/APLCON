

      DOUBLE PRECISION FUNCTION CHPROB(CHISQ,N) ! prob from chisquare
*     chi square probability for N degrees of freedom at CHISQ
*     =integral from 0 ... CHISQ   
      IMPLICIT NONE
      DOUBLE PRECISION CHISQ,DGAMIN
      INTEGER N
*     ...
      IF(CHISQ.LE.0.0) THEN
         CHPROB=1.0D0
      ELSE
         CHPROB=1.0D0-DGAMIN(0.5D0*DFLOAT(N),0.5D0*CHISQ)
      END IF 
      END 

      DOUBLE PRECISION FUNCTION DGAMIN(A,X)
*     incomplete gamma function P(a,x)
*     returns -1.0 for x < 0 or a =< 0
      IMPLICIT NONE
      DOUBLE PRECISION A,X,EPS,FPMIN
      INTEGER ITMAX,I,N
      PARAMETER (ITMAX=300,EPS=3.0E-7,FPMIN=1.0E-30)
      DOUBLE PRECISION AN,AP,B,C,D,DEL,H,SUM,DGAMML,GLN
*     ...
      IF(X.LT.0.0D0.OR.A.LE.0.0D0) THEN
         DGAMIN=-1.0D0  ! error
      ELSE IF(X.EQ.0.0) THEN
         DGAMIN=0.0D0
      ELSE IF(X.LT.A+1.0D0) THEN ! series representation
         GLN=DGAMML(A)         ! ln[Gamma(a)]
         AP=A  
         SUM=1.0D0/A
         DEL=SUM
         DO N=1,ITMAX
          AP=AP+1.
          DEL=DEL*X/AP
          SUM=SUM+DEL
          IF(ABS(DEL).LT.ABS(SUM)*EPS) GOTO 10
         END DO
         STOP 'DGAMIN:  error 1'
 10      DGAMIN=SUM*EXP(-X+A*LOG(X)-GLN)
      ELSE                     ! continued fraction representation
         GLN=DGAMML(A)         ! ln[Gamma(a)]
         B=X+1.0D0-A
         C=1.0D0/FPMIN
         D=1.0D0/B
         H=D
         DO I=1,ITMAX
          AN=-I*(I-A)
          B=B+2.0D0
          D=AN*D+B
          IF(ABS(D).LT.FPMIN) D=FPMIN
          C=B+AN/C
          IF(ABS(C).LT.FPMIN) C=FPMIN
          D=1.0D0/D
          DEL=D*C
          H=H*DEL
          IF(ABS(DEL-1.0D0).LT.EPS) GOTO 20
         END DO
         RETURN
         STOP 'DGAMIN:  error 2'
 20      DGAMIN=1.0D0-EXP(-X+A*LOG(X)-GLN)*H
      END IF
      END

      DOUBLE PRECISION FUNCTION DGAMML(X)       ! ln[Gamma(x)]
      IMPLICIT NONE
      DOUBLE PRECISION X
      INTEGER J
      DOUBLE PRECISION SER,STP,TMP,XX,YY,COF(6)
      SAVE COF,STP
      DATA COF,STP/76.18009172947146D0,-86.50532032941677D0,
     +24.01409824083091D0,-1.231739572450155D0,0.1208650973866179D-2,
     +-.05395239384953D-5,2.5066282746310005D0/
*     ...
      XX=X
      YY=XX
      TMP=(XX+0.5D0)*LOG(XX+5.5D0)-XX-5.5D0
      SER=1.000000000190015D0
      DO J=1,6
        YY=YY+1.0D0
        SER=SER+COF(J)/YY
      END DO
      DGAMML=TMP+LOG(STP*SER/XX)
      END

      DOUBLE PRECISION function DINGAU(p) ! inverse normal distribution
*     routine written by john herrero
      DOUBLE PRECISION p,p_low,p_high
      DOUBLE PRECISION a1,a2,a3,a4,a5,a6
      DOUBLE PRECISION b1,b2,b3,b4,b5
      DOUBLE PRECISION c1,c2,c3,c4,c5,c6
      DOUBLE PRECISION d1,d2,d3,d4
      DOUBLE PRECISION z,q,r
*     ...
      a1=-39.6968302866538
      a2=220.946098424521
      a3=-275.928510446969
      a4=138.357751867269
      a5=-30.6647980661472
      a6=2.50662827745924
      b1=-54.4760987982241
      b2=161.585836858041
      b3=-155.698979859887
      b4=66.8013118877197
      b5=-13.2806815528857
      c1=-0.00778489400243029
      c2=-0.322396458041136
      c3=-2.40075827716184
      c4=-2.54973253934373
      c5=4.37466414146497
      c6=2.93816398269878
      d1=0.00778469570904146
      d2=0.32246712907004
      d3=2.445134137143
      d4=3.75440866190742
      p_low=0.02425
      p_high=1-p_low
      if(p.lt.p_low) goto 201
      if(p.ge.p_low) goto 301
201   q=dsqrt(-2*dlog(p))
      z=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/
     &((((d1*q+d2)*q+d3)*q+d4)*q+1)
      goto 204
301   if((p.ge.p_low).and.(p.le.p_high)) goto 202
      if(p.gt.p_high) goto 302
202   q=p-0.5
      r=q*q
      z=(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/
     &(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1)
      goto 204
302   if((p.gt.p_high).and.(p.lt.1)) goto 203
203   q=dsqrt(-2*dlog(1-p))
      z=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/
     &((((d1*q+d2)*q+d3)*q+d4)*q+1)
204   DINGAU=z
      return
      end

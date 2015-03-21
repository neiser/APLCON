#ifdef weg

      CALL FPS
      RADIUX=0.2E4
      RADIUY=0.1E-3  
      DO I=0,11
       PHI=3.1415926*FLOAT(I)/6.0
       X=RADIUX*COS(PHI)
       Y=RADIUY*SIN(PHI)
       WRITE(*,*) I,X,Y 
       CALL FPNXY('*',1,X,Y)
c       CALL FPXY('O',X,Y)
      END DO 
      CALL FPL

       DO I=1,31
        X=0.1*FLOAT(I-1)
        Y=EXP(-X)
c        CALL FPXY('*',Y,X)
       CALL FPNXY('*',1,Y,X)
      END DO 
      CALL FPL

      END
#endif

      SUBROUTINE FPL
**********************************************************************
*
*     SUBROUTINE FPL
*     --------------
*     PRINTER PLOTS OF DATA Y,X,... ON ONE PAGE (60 * 120  UNITS)
*                                                50 * 120  UNITS)
*
*     ADDITION OF DATA TO PLOT BUFFER (IN ANY ORDER)
*     ----------------------------------------------
*
*     CH = CHARACTER, E.G. '*'
*
*     CALL FPXY(CH,X,Y)               ADD POINT X,Y
*
*     CALL FPNXY(CH,N,X,Y)            ADD POLYGON X(1),Y(1)  X(2),Y(2)
*                                     ... (X(N),Y(N))
*                                     A SEQUENCE OF CALLS WITH N=1 AND
*                                     IDENTICAL CHARACTER CH WILL ALSO
*                                     BE PLOTTED AS A POLYGON. IF CH=' '
*                                     THE X Y VALUES ARE TAKEN INTO
*                                     ACCOUNT IN LIMIT CALCULATION,
*                                     BUT THE POINT IS NOT PLOTTED.
*
*     CALL FLX(X1,X2)                 LIMITS X1 AND X2
*
*     CALL FLY(Y1,Y2)                 LIMITS Y1 AND Y2
*
*
*     PRINTOUT
*     --------
*
*     CALL FPL                        PLOT DATA AND RESET BUFFER
*
*     CALL FPS                        RESET DATA (WITHOUT PLOTTING)
*
*     CALL FPLUN(LUN)                 CHANGE PRINT-UNIT TO LUN
*                                     LUN = 0 MEANS SUPPRESS PRINTOUT
*
*     ALL  ARGUMENTS  ARE  INPUT  ARGUMENTS.   SCALES   ARE   DETERMINED
*     AUTOMATICALLY, TO INCLUDE ALL DATA OR THE LIMITS CAN BE DEFINED
*     BY USER CALLS (FLX, FLY).           FPL USES A  BUFFER  WITH  1602
*     WORDS, WHICH IS SUFFICIENT FOR 800 COORDINATES OF ONE  POLYGON  OR
*     400 COORDINATE PAIRS FOR ONE CALL OF FPNXY.
*
*     EXAMPLE: PLOT OF FUNCTION EXP(-X) FOR X BETWEEN 0 AND 3
*
*         DO 10 I=1,31
*         X=0.1*FLOAT(I-1)
*         Y=EXP(-X)
*      10 CALL FPNYX('*',1,Y,X)
*         CALL FPL
**********************************************************************

*     PARAMETER
      PARAMETER (LINENR = 40)
*     ARGUMENTS
         REAL X(*),Y(*)
         CHARACTER*1 CH
*     STORAGE FOR CHARACTERS
         CHARACTER*70 LINE,SLINE
c         CHARACTER*120 LINE,SLINE
         CHARACTER*1   C
         COMMON/CHACOM/LINE(LINENR),C
         EQUIVALENCE (SLINE,LINE(1))
*     ARRAYS FOR SCALE NUMBERS
         REAL RNUM(8)
         CHARACTER*10 NUMB(8)
*     STORAGE FOR COORDINATES
         COMMON/FICOM/IST(1602)
         REAL         RST(1602)
         EQUIVALENCE (IST(1),RST(1))
         COMMON/CFPL1/XA,XB,YA,YB,XLA,XLB,YLA,YLB
         DATA    ILIM/1602/,I/0/,ILAST/0/,ILX/0/,ILY/0/,ILZ/0/
         DATA IC/0/,ICLAST/0/,LUP/6/,IM/0/

*     ...
      IF(IM.EQ.0.OR.LUP.EQ.0) GOTO 60
      CALL FPLOST(I,ILX,ILY,ILZ)  
      IC=0

*
      XSCAL=MAX(ABS(XLA),ABS(XLB))      
      DO J=1,8
       RNUM(J)=(XLA*FLOAT(8-J)+XLB*FLOAT(J-1))/7.0
        IF(ABS(RNUM(J)).LT.10.0D-5*XSCAL) RNUM(J)=0.0D0
       CALL GFORM(RNUM(J),NUMB(J))
      END DO 
      WRITE(LUP,105) (NUMB(J),J=2,7)
      YSCAL=MAX(ABS(YLA),ABS(YLB))
      DO J=1,LINENR/10+1
      RNUM(J)=(YLB*FLOAT(LINENR/10+1-J)+YLA*FLOAT(J-1))/FLOAT(LINENR/10)
       IF(ABS(RNUM(J)).LT.10.0D-5*YSCAL) RNUM(J)=0.0D0
       CALL GFORM(RNUM(J),NUMB(J))
      END DO
      WRITE(LUP,102)         ! axis 

      DO J=1,LINENR                ! lines
       C='|'
       IF(MOD(J,10).EQ.1) C='T'
       IF(MOD(J,10).EQ.1) THEN
          K=1 
 13       IF(NUMB(1+J/10)(K:K).EQ.' ') THEN
             K=K+1
             GOTO 13
          END IF
          LINE(J)(1:10)=NUMB(1+J/10)(K:10)
       END IF 
       IF(MOD(J,10).EQ.6) C='T'
       WRITE(LUP,104) C,LINE(J),C  ! plot array
      END DO
      WRITE(LUP,102)         ! axis 
      DO J=1,8
       RNUM(J)=(XLA*FLOAT(8-J)+XLB*FLOAT(J-1))/7.0
        IF(ABS(RNUM(J)).LT.10.0D-5*XSCAL) RNUM(J)=0.0D0
       CALL GFORM(RNUM(J),NUMB(J)) 
      END DO
      WRITE(LUP,105) (NUMB(J),J=2,7)

      ENTRY FPS ! reset buffer
   60 ILX=0
      ILY=0
      ILZ=0
      IM=0
      IC=0
      I =0
      LCN=0
      RETURN

      ENTRY FPXY(CH,X,Y)                       ! store single pair
      N=1
      LCN=0
      IF(I+6.GT.ILIM) CALL FPLOST(I,ILX,ILY,ILZ) ! store user Y limits
      GOTO 85

      ENTRY FPNXY(CH,NN,X,Y)                   ! store NN pairs
      IF(NN.LE.0) GOTO 100
      N=NN
      IF(I+4+2*N.GT.ILIM) CALL FPLOST(I,ILX,ILY,ILZ)
      IF(N.NE.1.OR.LCN.NE.1)     GOTO 80
      IF(II.GT.I) GOTO 80
      IF(IST(II-1).NE.ICHAR(CH)) GOTO 80
      IF(I+2.GT.ILIM) GOTO 100
      RST(I+1)=X(1)                            ! append ro previous area
      RST(I+2)=Y(1)
      I=I+2
      IST(II)=IST(II)+1
      XA=MIN(XA,X(1))
      XB=MAX(XB,X(1))
      YA=MIN(YA,Y(1))
      YB=MAX(YB,Y(1))
      RETURN

   80 LCN=NN
   85 IF(IM.EQ.0) THEN
         YA=Y(1)
         XA=X(1)
         XB=X(1)
         YB=Y(1)
         IM=1
      END IF
      IF(I+4.GT.ILIM) GOTO 100
      IF(CH.EQ.' ') GOTO 100
*     STORAGE OF X,Y PAIRS
*     ICHAR / NR OF PAIRS / X / Y / X / Y / ...

      IST(I+1)=ICHAR(CH)
      IST(I+2)=0
      I=I+2
      II=I
      DO 90 J=1,N
      IF(I+2.GT.ILIM) GOTO 100
      RST(I+1)=X(J)
      RST(I+2)=Y(J)
      I=I+2
      IST(II)=IST(II)+1
      XA=MIN(XA,X(J))
      XB=MAX(XB,X(J))
      YA=MIN(YA,Y(J))
   90 YB=MAX(YB,Y(J))
      IF(LCN.EQ.1) RETURN
      GOTO 100

      ENTRY FLX(X1,X2) ! store user limits in X
      XLA=MIN(X1,X2)
      XLB=MAX(X1,X2)
      IF(XLA.EQ.XLB) XLB=XLA+1.0
      ILX=1
      RETURN

      ENTRY FLY(Y1,Y2) ! store user limits in Y
      YLA=MIN(Y1,Y2)
      YLB=MAX(Y1,Y2)
      IF(YLA.EQ.YLB) YLB=YLA+1.0
      ILY=1
      RETURN

  100 LCN=0
      RETURN

      ENTRY FPLUN(LUN) ! print unit
      LUP=LUN
      RETURN
  102 FORMAT('+',7('L----L----'),'+')
  104 FORMAT(A1,A70,A1)
  105 FORMAT(6X,6A10)
      END



      SUBROUTINE FPLOST(ILAST,ILX,ILY,ILZ) ! convert num to char 
 
      PARAMETER (LINENR = 40)
*     STORAGE FOR CHARACTERS
         CHARACTER*70 LINE,SLINE
         CHARACTER*1   C
         COMMON/CHACOM/LINE(LINENR),C
         EQUIVALENCE (SLINE,LINE(1))
*     STORAGE FOR COORDINATES
         COMMON/FICOM/IST(1602)
         REAL         RST(1602)
         EQUIVALENCE (IST(1),RST(1))
      COMMON/CFPL1/XA,XB,YA,YB,XLA,XLB,YLA,YLB
*     ...
      IF(ILX.EQ.0) THEN
         CALL BROUND(XA,XB,70,XLA,XLB)     ! scale in X
         ILX=1
      END IF
      IF(ILY.EQ.0) THEN
         CALL BROUND(YA,YB,40,YLA,YLB)  ! scale in Y
         ILY=1
      END IF

      IF(ILZ.EQ.0) THEN
         ILZ=1
         DO J=1,LINENR
          LINE(J)=' '        ! blank buffer
         END DO 
      END IF

      I=0
   10 IF(I.GE.ILAST) THEN
         ILAST=0 ! clear numerical buffer
         RETURN 
      END IF
      C=CHAR(IST(I+1))       ! loop over numerical data
      M=IST(I+2)
      I=I+2

      IF(M.EQ.1) THEN
         X1=RST(I+1)         ! point
         Y1=RST(I+2) 
         RJ=1.0+ FLOAT(LINENR)*(Y1-YLA)/(YLB-YLA)
         IF(RJ.LT.0.99.OR.RJ.GT.FLOAT(LINENR)+1.01) GOTO 10
         LJ=RJ
         LJ=MIN(LINENR,MAX(1,LJ))
         RI=1.0+70.0*(X1-XLA)/(XLB-XLA)
         IF(RI.LT.0.99.OR.RI.GT.121.01) GOTO 10
         LI=RI
         LI=MIN(70,MAX0(1,LI))
         LINE(LINENR+1-LJ)(LI:LI)=C
         I=I+2
      ELSE
         DO J=1,M-1          ! polygon
          X0=RST(I+1)
          Y0=RST(I+2)
          X2=RST(I+3)
          Y2=RST(I+4)  
          SI=2.0+70.0*(X2-X0)/(XLB-XLA)
          SJ=2.0+FLOAT(LINENR)*(Y2-Y0)/(YLB-YLA)
          KP=SQRT(SI*SI+SJ*SJ)+1.5
          IF(KP.GT.1000) KP=2
          DO K=1,KP
           XX=(X0*FLOAT(KP-K)+X2*FLOAT(K-1))/FLOAT(KP-1)
           YY=(Y0*FLOAT(KP-K)+Y2*FLOAT(K-1))/FLOAT(KP-1)
           LJ=1.0+FLOAT(LINENR)*(YY-YLA)/(YLB-YLA)
           IF(LJ.GE.1.AND.LJ.LE.LINENR) THEN
              LI=1.0+70.0*(XX-XLA)/(XLB-XLA)
              IF(LI.GE.1.AND.LI.LE.70) LINE(LINENR+1-LJ)(LI:LI)=C
           END IF 
          END DO 
         I=I+2
        END DO 
       I=I+2
      END IF
      GOTO 10
      END

      SUBROUTINE BROUND(XA,XB,NB,A,B)
*
*     Determine rounded limits A and B of a region with NB bins, which
*     includes the values XA and XB.
*
*                 -- -- --
*     CALL ROUNDB(XA,XB,NB,A,B)
*                          - -
*
*     Example: CALL BROUND(0.0172,2.5726,60,A,B)
*              resulting A=0.0  B=3.0
*
      IMPLICIT NONE
      INTEGER I,NB,NEX   
      REAL XA,XB,A,B,BINS(10),ZA,ZB,BIN,FEX   
      DATA BINS/0.1,0.15,0.2,0.25,0.4,0.5,0.6,0.8,1.0,1.5/
*     ...
      ZA=XA
      ZB=XB
      BIN=(ZB-ZA)/FLOAT(NB)
      IF(BIN.LE.0.0) BIN=1.0E-4
      NEX=INT(ALOG10(BIN)+100.0)-99
      FEX=10.0**NEX  
      BIN=BIN/FEX
      DO I=1,9
       IF(BINS(I).GE.BIN) GOTO 10
      END DO
      I=9
   10 BIN=BINS(I)*FEX
      A=10.0*BIN*AINT(ZA/(10.0*BIN))
      IF(A.GT.ZA) A=A-10.0*BIN
      B=A+FLOAT(NB)*BIN
      IF(B.LT.ZB.AND.I.NE.10) THEN
         I=I+1
         GOTO 10
      END IF
      END

      SUBROUTINE GFORM(VALUE,TEXT)
      CHARACTER*(*) TEXT,GTEXT(5)*10
*     ...
      DO I=5,2,-1
       IF(I.EQ.1) THEN
          WRITE(GTEXT(1),'(1P,G10.1)') VALUE
          GOTO 50
       ELSE IF(I.EQ.2) THEN
          WRITE(GTEXT(2),'(1P,G10.2)') VALUE
          GOTO 50
       ELSE IF(I.EQ.3) THEN
          WRITE(GTEXT(3),'(1P,G10.3)') VALUE
       ELSE IF(I.EQ.4) THEN
          WRITE(GTEXT(4),'(1P,G10.4)') VALUE
       ELSE IF(I.EQ.5) THEN
          WRITE(GTEXT(5),'(1P,G10.5)') VALUE
       END IF
c       WRITE(*,*) I,GTEXT(I)
       IF(GTEXT(I)(1:1).EQ.'*') GOTO 40
       DO IA=1,10
        IF(GTEXT(I)(IA:IA).NE.' ') GOTO 10
       END DO
       IA=11
 10    IE=0
       DO IB=IA,10
        IF(GTEXT(I)(IB:IB).EQ.' ') GOTO 30
        IF(GTEXT(I)(IB:IB).EQ.'E'.OR.GTEXT(I)(IB:IB).EQ.'e') GOTO 20
       END DO
       IB=11
 20    IE=IB   ! exponent E
       IB=IB+4
 30    IZ=IB-1 ! last non-blank
       IF(IE.EQ.0) THEN ! no exponent
          IF(GTEXT(I)(IZ:IZ).NE.'0') GOTO 50
       ELSE             ! exponent
          IF(GTEXT(I)(IZ-4:IZ-4).NE.'0') GOTO 50
       END IF
 40    CONTINUE
      END DO
 50   I1=(10-IZ+IA+1)/2
      I2=I1+IZ-IA
      TEXT=' '
      TEXT(I1:I2)=GTEXT(I)(IA:IZ)
c      WRITE(*,*) VALUE
      END

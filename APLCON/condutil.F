
      SUBROUTINE SMTOS(V,I,W,J,N)
*     Copy symmetric N-by-N matrix or a N-by-N submatrix of a  symmetric
*     matrix to another symmetric matrix
*
*                     - -     -
*          CALL SMTOS(V,I,W,J,N)
*                         - -
*
*     N rows and columns of the matrix V, starting from diagonal element
*     (i,i), are copied to the matrix W, starting  in  diagonal  element
*     (j,j). thus if a complete symmetric matrix has to be copied, i=j=1
*     has to be used.
*    
      IMPLICIT NONE
      INTEGER I,J,K,L,N,IM,JM
      DOUBLE PRECISION V(*),W(*)
*     ...
      WRITE(*,*) 'SMTOS ',I,J,N
      IM=(I*I+I)/2-1
      JM=(J*J+J)/2-1
      DO K=1,N
       DO L=1,K
        IM=IM+1
        JM=JM+1
        W(JM)=V(IM)
       END DO
       IM=IM+I-1
       JM=JM+J-1
      END DO 
      END

      INTEGER FUNCTION IJSYM(I,J) ! index (I,J)=(J,I) in symmetric matrix
      IMPLICIT NONE
      INTEGER I,J !,NCOUNT
c      DATA NCOUNT/0/
      IF(I.LE.J) THEN
         IJSYM=(J*J-J)/2+I
      ELSE
         IJSYM=(I*I-I)/2+J
      END IF
c      NCOUNT=NCOUNT+1
c      WRITE(*,*) 'IJSYM ',NCOUNT
      END

      SUBROUTINE DBMINV(W,B,NARG,MB, NRANK, AUX,QNEXT) ! matrix inversion
*     Obtain solution of a system of linear equations V *  X  =  B  with
*     symmetric matrix V and inverse (for M =  1)  or  matrix  inversion
*     only (for M = 0).
*
*                   - - - -
*        CALL SMINV(W,B,N,M,NRANK,AUX,QNEXT)
*                   - -     -----
*
*           W = symmetric N-by-N matrix in symmetric storage mode
*               W(1) = W11, W(2) = W12, W(3) = W22, W(4) = W13, . . .
*               replaced by inverse matrix
*           B = N-vector   (for M = 0 use a dummy argument)
*               replaced by solution vector
*           M = see above
*
*
*     Method of solution is by elimination selecting the  pivot  on  the
*     diagonal each stage. The rank of the matrix is returned in  NRANK.
*     For NRANK ne N, all remaining  rows  and  cols  of  the  resulting
*     matrix V and the corresponding elements of  B  are  set  to  zero.
*     SMINV can be used for a dimension up to 100 (see INVCDR).
*
      IMPLICIT NONE
      DOUBLE PRECISION VKK 
C      INTEGER NEXT(*)
C      DOUBLE PRECISION W(*),B(*),VX(*),A(*),AUX(*),QNEXT(*)
      DOUBLE PRECISION W(*),B(*),A(*),AUX(*),QNEXT(*)
      DOUBLE PRECISION SUM,EPS,VJK
      INTEGER NX,NF
      INTEGER I,J,IJ,IA,JFIRST,JLAST,NARG,MB,NRANK,
     +        K,LAST,JK,NMEAS,N,JJ,L,KK,M,JL,LK  

      LOGICAL SOLVE
       
      DATA EPS/1.0D-6/
*     ...
*     normal inversion entry 
      N=NARG                     ! dimension parameter 

      JFIRST=1                  
      DO I=1,N
       QNEXT(I)=I+1               ! set all "next" pointer
       AUX(I)=ABS(W((I*I+I)/2))  ! save abs of diagonal elements
      END DO
      QNEXT(N)=-1                 ! end flag
      NRANK=0
c      GOTO 10
c
c      ENTRY DCMINV(W,B,NARG,MB, NRANK, AUX,QNEXT)
c      N=NARG 
      
*     common code for inversion and (M=1) solution of matrix equation  

 10   SOLVE=.TRUE. 
      IF(MB.EQ.0) SOLVE=.FALSE.  ! solution flag

*     loop begin (loop on all remaining rows/cols)

      DO I=1,N                  ! loop on all remaining elements
       VKK=0.0D0                ! search for pivot element
       K=0                      ! pivot index
       J=JFIRST                 ! first candidate index
       LAST=0
 20    IF(J.GT.0) THEN          ! test for linearity and zero matrix 
          JJ=(J*J+J)/2          ! diagonal index
          IF(ABS(W(JJ)).GT.MAX(ABS(VKK),EPS*AUX(J))) THEN
             VKK=W(JJ)          ! largest pivot candidate so far
             K=J                ! index of largest
             L=LAST
          END IF
          LAST=J
          J=QNEXT(J)             ! index of next candidate
          GOTO 20
       END IF
       IF(K.NE.0) THEN          ! pivot element found - proceed
          NRANK=NRANK+1         ! increase rank counter
          KK=(K*K+K)/2
          IF(L.EQ.0) THEN
             JFIRST=QNEXT(K)    ! new first index
          ELSE
             QNEXT(L)=QNEXT(K)  ! bridge used index
          END IF
          QNEXT(K)=0.0D0        ! reset used index
          NRANK=NRANK+1         ! increase rank

          VKK    =1.0/VKK       ! invert pivot 
          W(KK)  =-VKK          
          IF(SOLVE) B(K)=B(K)*VKK
          JK     =KK-K
          JL     =0
          DO J=1,N              ! elimination
           IF(J.EQ.K) THEN
              JK=KK
              JL=JL+J
           ELSE
              IF(J.LT.K) THEN
                 JK=JK+1
              ELSE
                 JK=JK+J-1
              END IF
              VJK  =W(JK)
              W(JK)=VKK*VJK
              IF(SOLVE) B(J) =B(J)-B(K)*VJK
              LK   =KK-K
              DO L=1,J
               JL=JL+1
               IF(L.EQ.K) THEN
                  LK=KK
               ELSE
                  IF(L.LT.K) THEN
                     LK=LK+1
                  ELSE
                     LK=LK+L-1
                  END IF
                  W(JL)=W(JL)-W(LK)*VJK
               END IF
              END DO
           END IF
          END DO
       ELSE                     ! no pivot candadate found - reset
          DO K=1,N
           IF(QNEXT(K).NE.0.0D0) THEN   ! undefined variable
              IF(SOLVE) B(K)=0.0D0 ! clear undefined vector element
              DO J=1,K
               IF(QNEXT(J).NE.0.0D0) W((K*K-K)/2+J)=0.0D0  ! clear matrix row/col
              END DO
           END IF
          END DO
          GOTO 30
       END IF

      END DO                    ! end of inversion loop

 30   DO I=1,(N*N+N)/2
       W(I)=-W(I)               ! finally reverse sign
      END DO
      RETURN 

*     special entry for partial inversion ******************************

C      ENTRY DUMINV(VX,A, W,B,NX,NF,MB, NRANK, AUX,QNEXT)
      ENTRY DUMINV(A, W,B,NX,NF,MB, NRANK, AUX,QNEXT)

*     ...  
      NRANK=0 
      N=NX+NF                   ! dimension parameter 

*     -VX(NX-sym) is already inserted in W(NX+NF-sym) ------------------ 

C*     insert -VX(NX-sym) and A(NX,NF) into W(NX+NF-sym) ----------------

      IJ=(NX*NX+NX)/2           ! number of elements of V(.) 
c      DO I=1,IJ
c       W(I)=-VX(I)              ! copy -VX(.) into W_11
c      END DO

      DO I=1,N
       QNEXT(I)=0.0D0           ! reset pointer
      END DO

      IA=0
      DO J=1,NF
       DO I=1,NX
        W(IJ+I)=A(IA+I)         ! copy A(.) into W_12
       END DO
       DO I=1,J
        W(IJ+NX+I)=0.0          ! reset last submatrix W_22 of W(.)
       END DO
       IJ=IJ+NX+J 
       IA=IA+NX
      END DO

*     distinguish between measured and unmeasured variables ------------

      JFIRST=0                  ! first index of measured variable 
      NMEAS=0                   ! number of measured variables
      DO I=1,NX
       IF(W((I*I+I)/2).LT.0.0) THEN  ! measured variable
          IF(JFIRST.EQ.0) THEN
             JFIRST=I           ! first index of measured variable
          ELSE
             QNEXT(JLAST)=I     ! insert index at previous index
          END IF
          JLAST=I               ! save index
          NMEAS=NMEAS+1
       END IF
      END DO
c      IF(JLAST.EQ.0) RETURN     ! nothing to do
      IF(JLAST.EQ.0) GOTO 10    ! nothing to do
      QNEXT(JLAST)=-1           ! stop index for last measured variable 


*     apply exchange algorithm to sub-matrices -------------------------

      DO I=NX+1,N               ! loop I over constraint equations

       J=JFIRST                 ! first index of unmeasured variable
       DO M=1,NMEAS             ! already inverted element index J
        SUM=0.0D0
        JK=(J*J-J)/2            ! index of diagonal element before
        DO K=1,NX
         IF(K.LE.J) JK=JK+1     ! index in j column
         IF(QNEXT(K).NE.0.0D0) SUM=SUM+W(JK)*W((I*I-I)/2+K)
         IF(K.GE.J) JK=JK+K     ! index in j row
        END DO
        AUX(J)=SUM              ! = A-row * VX-row/col
        J=QNEXT(J)              ! next index of unmeasured variable
       END DO

       DO K=I,N 
        SUM=0.0D0
        J=JFIRST                ! first index of unmeasured variable 
        DO M=1,NMEAS            ! already inverted element index J
         SUM=SUM+W((K*K-K)/2+J)*AUX(J) ! = A-row * H
         J=QNEXT(J)             ! next index of unmeasured variable
        END DO
        W((K*K-K)/2+I)=W((K*K-K)/2+I)+SUM   ! add to diagonal W_22
       END DO

       J=JFIRST                 ! first index of unmeasured variable
       DO M=1,NMEAS
        W((I*I-I)/2+J)=-AUX(J)  ! add to off-diagonal W_22
        J=QNEXT(J)              ! next index of unmeasured variable
       END DO

      END DO

*     set pointer for unmeasured variables ----------------------------- 

      JFIRST=0
      JLAST=0
      DO I=1,N                  
       IF(QNEXT(I).EQ.0.0D0) THEN ! unmeasured variable
          IF(JFIRST.EQ.0) THEN  
             JFIRST=I           ! first index of unmeasured variable
          ELSE
             QNEXT(JLAST)=I     ! next index of unmeasured variable
          END IF
          JLAST=I
       ELSE
          QNEXT(I)=0.0D0        ! reset index for measured variable
       END IF 
      END DO
c      IF(JLAST.EQ.0) RETURN    ! no unmeasured variable
      IF(JLAST.EQ.0) GOTO 10    ! no unmeasured variable
      QNEXT(JLAST)=-1           ! end flag
c      RETURN                   ! go to invert remaining parts  
      GOTO 10                   ! go to invert remaining parts  
      END 


      SUBROUTINE SMAVAT(V,A,W,N,M)
*     Multiply symmetric N-by-N matrix from the left with general M-by-N
*     matrix and from the right with the transposed of the same  general
*     matrix  to  form  symmetric  M-by-M   matrix   (used   for   error
*     propagation).
*                   - -   - -
*        CALL SMAVA(V,A,W,N,M)
*                       -
*                                  T
*         W   =   A   *   V   *   A
*        M*M     M*N     N*N     N*M
*
*
*        where V = symmetric N-by-N matrix
*              A = general N-by-M matrix
*              W = symmetric M-by-M matrix
*
      IMPLICIT NONE
      INTEGER N,M,I,IL,IJS,LKL,JK,IJ,K,L,J,LK
      DOUBLE PRECISION V(*),A(*),W(*)
      DOUBLE PRECISION CIK
*     ...
      DO I=1,(M*M+M)/2
       W(I)=0.0
      END DO
*
      IL=-N
      IJS=0
      DO I=1,M
       IJS=IJS+I-1
       IL=IL+N
       LKL=0
       DO K=1,N
        CIK=0.0
        LKL=LKL+K-1
        LK=LKL
        DO L=1,K
         LK=LK+1
         CIK=CIK+A(IL+L)*V(LK)
        END DO
        DO L=K+1,N
         LK=LK+L-1
         CIK=CIK+A(IL+L)*V(LK)
        END DO
        JK=K
        IJ=IJS
        DO J=1,I
         IJ=IJ+1
         W(IJ)=W(IJ)+CIK*A(JK)
         JK=JK+N
        END DO
       END DO
      END DO
      END

#ifdef mmm 
      DOUBLE PRECISION FUNCTION SCALXY(X,Y,N) ! scalar vector product
*     Scalar product of two vectors
*                - - -                 T
*        S = VXY(X,Y,N)           S = X  * Y (scalar product)
*
      IMPLICIT NONE
      INTEGER N,I,M,J
      DOUBLE PRECISION X(*),Y(*)
      DOUBLE PRECISION SUM
*     ...
      SUM=0.0D0
      M=N/4
      I=0
      DO J=1,M
       SUM=SUM+X(I+1)*Y(I+1)+X(I+2)*Y(I+2)+X(I+3)*Y(I+3)+X(I+4)*Y(I+4)
       I=I+4
      END DO
      DO J=I+1,N
       SUM=SUM+X(J)*Y(J)
      END DO
      SCALXY=SUM
      END
#endif

      DOUBLE PRECISION FUNCTION SCALXY(X,Y,N) ! scalar vector product
*     Scalar product of two vectors
*                - - -                 T
*        S = VXY(X,Y,N)           S = X  * Y (scalar product)
*
      IMPLICIT NONE
      INTEGER N,J ! I,M
      DOUBLE PRECISION X(*),Y(*)
      DOUBLE PRECISION SUM
*     ...
      SUM=0.0D0
      DO J=1,N
       SUM=SUM+X(J)*Y(J)
      END DO
      SCALXY=SUM
      END


      SUBROUTINE CFGMPR(LUP,A,N,M,TEXT)
*     prints the N-by-M general matrix A
*                        - - - ----
*        CALL CFGMPR(LUP,A,N,M,TEXT)
*
      IMPLICIT NONE
      INTEGER LUP,N,M,L,IJ,I
      DOUBLE PRECISION A(*)
      CHARACTER*(*) TEXT
*
*     general matrix or vector
*
*     ...
      IF(N.EQ.1.OR.M.EQ.1) THEN
*        special case of vector
         WRITE(LUP,101) N,M,TEXT
         WRITE(LUP,102) (A(L),L=1,N*M)
      ELSE 
*        matrix
         WRITE(LUP,103) N,M,TEXT
         IJ=0
         DO I=1,N
          WRITE(LUP,104) I,(A(IJ+L),L=1,M)
          IJ=IJ+M
         END DO
      END IF
      WRITE(LUP,102)
  101 FORMAT(1X,I4,' *',I3,'  Vector ',A/)
  102 FORMAT(10X,10G12.5)
  103 FORMAT(1X,I4,' *',I3,'  Matrix ',A/)
  104 FORMAT(1X,I4,5X,10G12.5/(10X,10G12.5/))
      END



      SUBROUTINE CFCORR(LUP,V,N) ! print matrix of correlation coefficients
      DOUBLE PRECISION V(*),RHO,PRD
      INTEGER MC(15)
      CHARACTER*9 MODTX
*     
      DO I=1,MIN(N,15)
       MC(I)=I  
      END DO
      WRITE(LUP,*) ' '
      WRITE(LUP,*) 'Matrix of correlation coefficients'
      WRITE(LUP,*) ' '
      MODTX='    i    '
      IF(N.GT.15) MODTX='mod(i,15)' 
      WRITE(LUP,103) MODTX,(MC(K),K=1,MIN(N,15))
c      WRITE(LUP,*) '   i'
      DO I=1,N
       L=0
       DO J=1,I
        RHO=0.0D0
        PRD=V(IJSYM(I,I))*V(IJSYM(J,J))
        IF(PRD.GT.0.0D0) RHO=V(IJSYM(I,J))/SQRT(PRD)
        IF(L.EQ.15) THEN
           IF(J.EQ.L+1) THEN
              WRITE(LUP,102) I,MC
           ELSE
              WRITE(LUP,101) (MC(K),K=1,L)
           END IF
           L=0
        END IF  
        L=L+1
        MC(L)=100.0*ABS(RHO)+0.5
        IF(RHO.LT.0.0) MC(L)=-MC(L)
       END DO 
       IF(I.LE.15) THEN 
          WRITE(LUP,102) I,(MC(K),K=1,L)
       ELSE
          WRITE(LUP,101) (MC(K),K=1,L)
       END IF
      END DO
      DO I=1,MIN(N,15)
       MC(I)=I  
      END DO
      WRITE(LUP,103) MODTX,(MC(K),K=1,MIN(N,15))
      WRITE(LUP,*) '(coefficients in %)'
      WRITE(LUP,*) ' ' 
 101  FORMAT(7X,   15I4) 
 102  FORMAT(I5,2X,15I4)
 103  FORMAT(A9,I2,14I4)
      END

#ifdef notyet
      SUBROUTINE CGCORR(LUP,V,N) ! print matrix of correlation coefficients
      DOUBLE PRECISION V(*),RHO,PRD
      INTEGER MC(15)
      CHARACTER*9 MODTX
*     
      DO I=1,MIN(N,15)
       MC(I)=I  
      END DO
      WRITE(LUP,*) ' '
      WRITE(LUP,*) 'Correlation coefficients'
      WRITE(LUP,*) ' '
      MODTX='    i    '
      IF(N.GT.15) MODTX='mod(i,15)' 
      WRITE(LUP,103) MODTX,(MC(K),K=1,MIN(N,15))
c      WRITE(LUP,*) '   i'
      DO I=1,N
       L=0
       DO J=1,I
        RHO=0.0D0
        PRD=V(IJSYM(I,I))*V(IJSYM(J,J))
        IF(PRD.GT.0.0D0) RHO=V(IJSYM(I,J))/SQRT(PRD)
        IF(L.EQ.15) THEN
           IF(J.EQ.L+1) THEN
              WRITE(LUP,102) I,MC
           ELSE
              WRITE(LUP,101) (MC(K),K=1,L)
           END IF
           L=0
        END IF  
        L=L+1
        MC(L)=100.0*ABS(RHO)+0.5
        IF(RHO.LT.0.0) MC(L)=-MC(L)
       END DO 
       IF(I.LE.15) THEN 
          WRITE(LUP,102) I,(MC(K),K=1,L)
       END IF
      END DO
      DO I=1,MIN(N,15)
       MC(I)=I  
      END DO
      WRITE(LUP,103) MODTX,(MC(K),K=1,MIN(N,15))
      WRITE(LUP,*) '(coefficients in %)'
      WRITE(LUP,*) ' ' 
 101  FORMAT(7X,   15I4) 
 102  FORMAT(I5,'|  ',15I4)
 103  FORMAT(A9,I2,14I4)
      END
#endif


      SUBROUTINE CFPRV(LUP,X,V,N) ! print X and uncertainties from V
*     Prints the n-vector X and uncertainties from V
*
      IMPLICIT NONE
      INTEGER LUP,N,II,I,JS
      DOUBLE PRECISION X(*),V(*),P(*),ERR
      LOGICAL NP
      CHARACTER*16 PANAME,TEXT*17
*     ...
      NP=.TRUE.
      WRITE(LUP,*)  ' '
      WRITE(LUP,101)' '
      GOTO 10

      ENTRY CFPRVP(LUP,X,V,P,N) ! print X, pulls and uncertainties from V 
      NP=.FALSE.
      WRITE(LUP,*) ' ' 
      WRITE(LUP,101) 'pull'
      WRITE(LUP,103)
 10   DO I=1,N
       II=(I*I+I)/2
       ERR=0.0
       IF(V(II).GT.0.0) ERR=SQRT(V(II))
       CALL APNAME(-I,PANAME)
       IF(NP) THEN 
          WRITE(LUP,102) PANAME,I,X(I),ERR
       ELSE 
          CALL APERTY(I,TEXT,JS)
          WRITE(LUP,102) PANAME,I,X(I),ERR,P(I),TEXT
       END IF
      END DO
c      WRITE(LUP,*) ' '
c      CALL CFCORR(LUP,V,N)
 101  FORMAT(21X,'i     value        std.dev',7X,A)
 102  FORMAT(1X,A16,I5,2G14.6,F8.2,2X,A17)  
 103  FORMAT(18X,'----  ------------  ------------','    ----')
      END 


      SUBROUTINE SIMSEL(X,VX,NY,LIST,Y,VY)
*     Copy subset of X,VX ---(list)--> Y,VY 
      IMPLICIT NONE
      DOUBLE PRECISION X(*),VX(*),Y(*),VY(*)
      INTEGER NY,LIST(NY)
      INTEGER I,J,IJ,IJSYM  
*     ... 
      IJ=0
      DO I=1,NY
       Y(I)=X(LIST(I)) 
       DO J=1,I
        IJ=IJ+1
        VY(IJ)=VX(IJSYM(LIST(I),LIST(J)))
       END DO
      END DO
      END 


      SUBROUTINE SIMTRN(X,VX,NX) ! extract normalization factor from V
*     extract normalization error
      DOUBLE PRECISION X(*),VX(*),EPS2
*     ...
      EPS2=0.0D0
      DO I=1,NX
       DO J=1,I-1
        IJ=IJSYM(I,J)
        EPS2=EPS2+VX(IJ)/(X(I)*X(J))
       END DO
      END DO
      EPS2=EPS2/DFLOAT((NX*NX-NX)/2)

      X(NX+1)=1.0D0
      DO I=1,NX
       VX(IJSYM(I,NX+1))=0.0D0
      END DO
      DO I=1,NX
       DO J=1,I
        IJ=IJSYM(I,J)
        VX(IJ)=VX(IJ)-EPS2*(X(I)*X(J))
       END DO
      END DO
      VX(IJSYM(NX+1,NX+1))=EPS2
      END

      SUBROUTINE DVELLP(APP,APQ,AQQ,   R1,R2,CPHI,SPHI) 
*     parameters of covariance ellipse
*     input:  APP,APQ,AQQ = matrix elements
*     result: R1,R2     = half axes
*             CPHI,SPHI = cos and sin phi
      IMPLICIT NONE
      DOUBLE PRECISION APP,APQ,AQQ,R1,R2,CPHI,SPHI,THETA,TPHI
*     ...
      IF(APQ.EQ.0.0D0) THEN
         CPHI=1.0D0
         SPHI=0.0D0
      ELSE
         THETA=0.5D0*(AQQ-APP)/APQ
         IF(THETA.EQ.0.0D0) THEN
            TPHI=1.0D0
         ELSE  
            TPHI=1.0D0/(THETA+SIGN(SQRT(1.0D0+THETA*THETA),THETA))
         END IF
         CPHI=1.0D0/SQRT(1.0D0+TPHI*TPHI)
         SPHI=CPHI*TPHI
      END IF   
      R1=SQRT(APP*CPHI*CPHI-2.0D0*APQ*CPHI*SPHI+AQQ*SPHI*SPHI)
      R2=SQRT(APP*SPHI*SPHI+2.0D0*APQ*CPHI*SPHI+AQQ*CPHI*CPHI)
c      WRITE(*,*) 'DVELLP:'
c      WRITE(*,*) SQRT(APP),SQRT(AQQ),APQ/(SQRT(APP)*SQRT(AQQ))
      END


      SUBROUTINE APNAME(I,NAME) ! store/get/rest variable name
*     __________________________________________________________________
*     CALL APNAME(0,'....')        reset
*     CALL APNAME(0,'text')        print text (up to 71 chars) and reset 
*     CALL APNAME(I,'name')  I>0   store name (right adjusted)
*     CALL APNAME(I,'name')  I<0   get name (up to 16 chars, right-adj.)
*     __________________________________________________________________
      IMPLICIT NONE  
      INTEGER NAN,NCH
      PARAMETER (NAN=128,NCH=16)
      INTEGER I,J,K,NUS,NTIT,NIN(NAN)
      CHARACTER*(*) NAME, VNAMES(NAN)*16, TITLE*71
      DATA NUS/0/ 
*     ...
      IF(I.EQ.0) THEN      ! reset
         NUS=0
         TITLE=' ' 
         TITLE=NAME
         IF(NAME.EQ.' ') RETURN
         NTIT=1
         DO J=1,71
          IF(TITLE(J:J).NE.' ') NTIT=J
         END DO
         WRITE(*,*) ' '
         WRITE(*,101) TITLE(1:NTIT)
         WRITE(*,*) ('_',J=1,NTIT)
         WRITE(*,*) ' ' 
      ELSE IF(I.GT.0) THEN ! store    
         DO J=1,NUS
          IF(NIN(J).EQ.I) GOTO 10
         END DO
         IF(NUS.EQ.NAN) RETURN   
         NUS=NUS+1
         J=NUS        
         NIN(J)=I
 10      VNAMES(J)=NAME
         DO K=NCH,1,-1 
          IF(VNAMES(J)(K:K).NE.' ') GOTO 20
         END DO
         RETURN
 20      VNAMES(J)=' ' 
         VNAMES(J)(NCH-K+1:NCH)=NAME
      ELSE                 ! get
         NAME=' '
         DO J=1,NUS
          IF(NIN(J)+I.EQ.0) THEN
             NAME=VNAMES(J)
             RETURN
          END IF
         END DO
      END IF
 101  FORMAT(1X,A)
      END


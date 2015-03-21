
      SUBROUTINE AVVSYS

      PARAMETER (NVARX=200,NEQSX=100)             ! max dimensions
      PARAMETER (NSYMX=(NVARX*NVARX+NVARX)/2)
      DOUBLE PRECISION X(NVARX),V(NSYMX),F(NEQSX)
*     ...

      WRITE(*,*) ' '
      WRITE(*,*) 'Branching rations with systematice errors'
      WRITE(*,*) ' (from Valassi, pages 403-404)'
      WRITE(*,*) 'A systematic error is introduced for the (B) exp'
      WRITE(*,*) 'data, which is either +100%, 0% or -100% correlated'      
      WRITE(*,*) ' '
      WRITE(*,*) 'Combining correlated measurements with ', 
     +           'systematic errors'

      DO ICASE=1,4 ! with different systematic errors

      DO IVERS=1,2 ! without/with lepton universality

      IF(IVERS.EQ.1) THEN
         WRITE(*,*) 'No lepton un iversality'
      ELSE
         WRITE(*,*) 'Lepton universality'
      END IF

      DO I=1,21
       V(I)=0.0D0
      END DO

      X(1)=0.1050  ! electron A
      X(2)=0.1350  ! electron A
      X(3)=0.0950  ! tau A
      X(4)=0.1400  ! tau B
      X(5)=0.0     ! systematic error mean
      X(6)=0.0     ! systematic error mean

      V(IJSYM(1,1))=1.0
      V(IJSYM(2,2))=0.0441
      V(IJSYM(3,3))=9.0
      V(IJSYM(4,4))=0.0441
      V(IJSYM(5,5))=8.96
      V(IJSYM(6,6))=8.96

      DO I=1,21
       V(I)=1.0D-4*V(I) ! scale matrix elements by 1.0D-2 ** 2
      END DO

      NVAR=5
      IF(ICASE.EQ.2) NVAR=6 
      IF(ICASE.EQ.4) NVAR=4
      NEQS=1+IVERS 
      CALL APLCON(NVAR,NEQS)
      LUP=6
      IPR=6
      CALL APRINT(LUP,IPR)
      IF(ICASE.EQ.1) THEN
         WRITE(*,*) 'Positive correlation'
      ELSE IF(ICASE.EQ.2) THEN
         WRITE(*,*) 'No correlation'
      ELSE IF(ICASE.EQ.3) THEN
         WRITE(*,*) 'Negative correlation'
      ELSE IF(ICASE.EQ.4) THEN
         WRITE(*,*) 'no systematic errors'
      END IF
      CALL APNAME(1,'BR to e (A)')
      CALL APNAME(2,'BR to e (B)')
      CALL APNAME(3,'BR to tau (A)')
      CALL APNAME(4,'BR to tau (B)')
      CALL APNAME(5,'sys error tau (B)')
      CALL APNAME(6,'sys error tau (B)')
      CALL CIPRV(LUP,X,V,NVAR) ! print initial status 
      CALL CFCORR(LUP,V,NVAR)  ! print correlation matrix
c      CALL APROFL(4,0)
c      CALL APROFL(1,2) 

*     Fit loop =========================================================

 10   IF(ICASE.EQ.1) THEN
         F(1)=X(1)-(X(2)+X(5))
         F(2)=X(3)-(X(4)+X(5))
      ELSE IF(ICASE.EQ.2) THEN
         F(1)=X(1)-(X(2)+X(5))
         F(2)=X(3)-(X(4)+X(6))
      ELSE IF(ICASE.EQ.3) THEN
         F(1)=X(1)-(X(2)+X(5))
         F(2)=X(3)-(X(4)-X(5))
      ELSE IF(ICASE.EQ.4) THEN
         F(1)=X(1)-X(2)
         F(2)=X(3)-X(4)
      END IF
      IF(IVERS.EQ.2) F(3)=X(1)-X(3)
      CALL APLOOP(X,V,F, IREP)
      IF(IREP.LT.0) GOTO 10
      END DO ! univers
      END DO ! cases

      END 

  
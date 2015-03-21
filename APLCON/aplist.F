
      SUBROUTINE APLIST(NVAR,X,V)    ! print index-value-error 
      DOUBLE PRECISION X(NVAR),V(*)
*     ...
      WRITE(*,*) ' '
      WRITE(*,101)
      DO I=1,NVAR
       ERR=SQRT(V(IJSYM(I,I)))
       WRITE(*,102) I,X(I),ERR 
      END DO
 101  FORMAT('   index    X(index)     error(index)')
 102  FORMAT(I8,F12.3,F15.3)
      END   

      SUBROUTINE APPULL(PULLS)
      DOUBLE PRECISION PULLS(*)
#include "comcfit.inc"
#include "nauxfit.inc" 
*     ...
      DO I=1,NX
       PULLS(I)=AUX(INDPU+I)
      END DO
      END


      
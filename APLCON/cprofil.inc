
*     Common for profile likelihood analysis with APLCON
*
*     NDEXT +1     is start of save area for primary X(.) and VX(.)
*     ILR    =     index of current profile point
*     NLR    =     number of profile points
*     ISECA  =     index of current profile analysis
*     NSECA  =     number of profile analyses
*     NFADD  =     max number of additional constraints
*     IPF    =     index of profile variable
*     NPSEC(2,.)   indices for profile analyses
*     XL,XR  =     limits for profile analysis
*     XLR(.) =     X values for profile
*     YLR(.) =     Y values for profile    
*     FLR(.) =     chi^2 - reference value for profile   
*     CONSTR =     fixed value 
*     CONSTX,Y     fixed values
*     CENTER =     fitted value (minimum)
*     SIGMAX =     parabolic X error
*     SIGMAY =     parabolic Y error
*     CENTEX =     X center 
*     CENTEY =     Y center
*     RED    =     reduction factor, e.g. 1/2
*     REDL   =     reduction factor left
*     REDR   =     reduction factor right
*     CPHI   =     cos phi of eigenvector
*     SPHI   =     sin phi  
*     NSTAR  =     number of directions
*     ISTAR  =     direction index
*     CRT    =     cos of rotation transformation
*     SRT    =     sin of rotation transformation
*     CTH    =     cos rotation
*     STH    =     sin rotation 
*     CHT    =     copy of cos rotation  
*     XCONT(24,3)  x contur data, max 12 points, 3 contours
*     YCONT(24,3)  y contur data, max 12 points, 3 contours
      INTEGER MLR,MSECA,MMLR
      PARAMETER (MLR=36,MSECA=100,MMLR=MLR+MLR+1)
      INTEGER ILR,NLR,ISECA,NSECA,NFADD,NDEXT,IPF,NPSEC
      INTEGER IPFX,IPFY,NSTAR,ISTAR 
      DOUBLE PRECISION XL,XR,CLR,XLR,YLR,FLR,CONSTR,CONSTX,CONSTY,CHIREF 
      DOUBLE PRECISION CSP,CENTER,SIGMAX,SIGMAY,CENTEX,CENTEY,RED,
     + REDL,REDR
      DOUBLE PRECISION CPHI,SPHI,CRT,SRT,CTH,STH,CHT
      REAL XCONT,YCONT,XC,YC       
      
      COMMON/CPROFL/XL,XR,CLR(-MLR:+MLR),
     +              XLR(-MLR:+MLR),YLR(-MLR:+MLR),FLR(-MLR:+MLR),
     +              CONSTR,CONSTX,CONSTY,
     +              CHIREF,CSP(5,MMLR),CENTER,SIGMAX,SIGMAY,
     +              CENTEX,CENTEY,RED,REDL,REDR,
     +              CPHI,SPHI,CRT,SRT,CTH,STH,CHT,
     +              XCONT(24,3),YCONT(24,3),XC(101),YC(101),
     +              NPSEC(2,MSECA),
     +              ILR,NLR,ISECA,NSECA,NFADD,IPF,NDEXT,
     +              IPFX,IPFY,NSTAR,ISTAR 




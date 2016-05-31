*     =========================== MACRO ================================
*     Parameter statement for basic dimensions
      INTEGER NADFS,NX,NF,NUM,IFLG,INIT,LUNSIM,IPR,IPAK,NFPRIM,
     +        ISTAT,INDST,INDLM,  NDENDA,NDENDE,NDACTL,
     +        INDTR,INDFC,INDHH,INDXS,INDDX,INDXP,INDRH,INDWM,
     +        INDIA,NDTOT,INDQN,NAUXC,NCASE,INDCF,INDPU,INDAS,INDVS,
     +        ICNT,NXF,MXF,NDF,IUNPH,NCST,ITER,NCALLS,NDPDIM,ITERMX,
     +        NDTOTL 
*     Definition of common for APLCON/ERRPRP/SIM... subroutines
      DOUBLE PRECISION EPSF,EPSCHI,CHISQ,FTEST,FTESTP,CHSQP,FRMS,FRMSP
      DOUBLE PRECISION DERFAC,DECXP,DERUFC,DERLOW,WEIGHT,PENALT
      REAL TAB
      COMMON/SIMCOM/EPSF,EPSCHI,CHISQ,FTEST,FTESTP,CHSQP,FRMS,FRMSP,
     +      DERFAC,DECXP,DERUFC,DERLOW,WEIGHT,PENALT,
     +      NADFS,NX,NF,NUM,IFLG,INIT,LUNSIM,IPR,NCASE,IPAK,INDCF,
     +      ISTAT,INDST,INDLM,    NDENDA,NDENDE,NDACTL,INDAS,INDVS,
     +      INDTR,INDFC,INDHH,INDXS,INDDX,INDXP,INDRH,INDWM,INDIA,
     +      NDTOT,ICNT,NXF,MXF,NDF,IUNPH,NCST,ITER,NCALLS,NDPDIM,INDQN,
     +      ITERMX,NAUXC,INDPU,NFPRIM,NDTOTL,  
     +      TAB(1000,10) 

*     NADFS  = 0   initial (primary) fit (set by APLCON)
*            = 1   1-parameter analysis  (set by BPLCON)
*            = 2   2-parameter analysis      "..."

*     
*     NX       = number of parameters
*     NF       = number of constraint equations
*     NUM      = flag for numerical differentiation
*                = 0  analytical derivatives
*                = 1  numerical derivatives
*                = 2  numerical derivatives, compare with analytical
*     IFLG     = flag for first case (check derivative)
*     INIT     = flag 
*                = 1  deriving
*                = 0  else
*     LUNSIM   = printout unit (see parameter statement)
*     IDEBUG   = debug flag value
*                = 0   no printout
*                = 1   only warnings are printed (default)
*                = 2,3 more and more printout
*
*     EPSF     = required |F| accuracy
*     ISTAT    = status flag
*                =0  
*                =1  derivative matrix finished
*                =2  X(.) corrections applied
* 
*     XL(2,.)  = lower and upper values of parameters
*     ST(.)    = step sizes for numerical differentiation
*     FC(.)    = central values of parameters
*     H(.)     = copy 
*
*     A        = derivative matrix a/flags during matrix inversion/pulls
*     A(NX,NF)
*
*     ******************************************************************
*
*     ND       = number of degrees of freedom
*              = number of constraint equations minus number of 
*                unmeasured parameters 
*     CHSQ     = chi square      
*
*     DERFAC   = factor for standard deviation in numerical derivatives
*
*     NDENDE   = index of last used word incl. single-precision array
*     NDACTL   = index of actual single-precision array


*     =========================end=of=macro=============================



*     Auxiliary array for all arrays used in APLCON 
      INTEGER    NAUX
      PARAMETER (NAUX=125 000)   ! 1 Mega byte
      DOUBLE PRECISION AUX
      COMMON/NAUXCM/   AUX(NAUX)


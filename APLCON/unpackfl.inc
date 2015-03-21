
*     unpackfl.inc = code for flag unpacking

#ifdef comment
      Flags for variables in TRFLAG(.), range/values in XL(2,.)
   
      Transformation flag NTVAR:

      1     inverse transformation      no limit
      2     Poisson distributed         range > 0
      3     Binomial distributed        range > 0 and < 1; XL(2,.)=N
      4     log-normal transformation   range > 0
      5     sqrt transformation         range > 0
      6     power transformation        range > 0; xl(2,.)=power


      Flag for M-estimate weighting NTMES:

      0     no weighting
      1     unmeasured
      2     Cauchy 
      3     Huber
      4 



      Derivative flag NTDER:

      0    unknown
      1    unrelated (all derivatives = 0)
      2    single variable with constant derivative
      3    linear (derivatives constant)


      Inequality flag NTINE:

      0    variable 
      1    fixed by user     
      2    fixed by inequality constraint
       

      Limit flag NTLIM 

      0    unlimited
      1    positive
      2    lower limit
      3    upper limit
      4    lower & upper limit 

      
      Profile flag NTPRF 
 
      0    no profile
      1    1-dim profile
      2    2-dim profile
      3    1- & 2-dim profile


#endif
      NTRFL=AUX(INDTR+IPAK)        ! get packed flags
      NTVAR=MOD(NTRFL,10)          ! transformation flag
      NTMES=MOD(NTRFL/10,10)       ! M-estimate flag
      NTDER=MOD(NTRFL/100,10)      ! derivative type flag
      NTINE=MOD(NTRFL/1000,10)     ! inequality flag
      NTLIM=MOD(NTRFL/10000,10)    ! limit flag
      NTPRF=MOD(NTRFL/100000,10)   ! profile flag

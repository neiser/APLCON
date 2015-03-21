subroutine C_APLCON_APLCON(NVAR,MCST) bind(c)
  use iso_c_binding
  integer(c_int), value :: NVAR, MCST
  CALL APLCON(NVAR,MCST)  
end subroutine

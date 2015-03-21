cmake_minimum_required (VERSION 2.6)
project(APLCON Fortran)
set(SRCS_LIB
  a12prof.F
  aplist.F
  aploop.F
  aplprint.F
  aplusopt.F
  chprob.F
  condutil.F
  dcsplines.F
  fpl.F
  )
add_library(aplcon ${SRCS_LIB})

set(SRCS_TEST
  test/avlass.F
  test/avvsys.F
  test/avmain.F
  )
add_executable(test_avmain ${SRCS_TEST})
target_link_libraries(test_avmain aplcon)

cmake_minimum_required (VERSION 2.6)
project(APLCON Fortran)

# those flags (from original Makefile) are crucial to make the program run properly
# otherwise you end up in an endless loop
set(CMAKE_Fortran_FLAGS "-Wunused  -fno-automatic -fno-backslash -O")

# build the library (static)
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

# link the test program against it
set(SRCS_TEST
  test/avlass.F
  test/avvsys.F
  test/avmain.F
  )
add_executable(test_avmain ${SRCS_TEST})
target_link_libraries(test_avmain aplcon)

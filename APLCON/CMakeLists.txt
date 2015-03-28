cmake_minimum_required (VERSION 2.6)
project(APLCON Fortran)

if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  message(FATAL_ERROR "Only GNU Fortran compiler supported")
endif()

# set default build type if unspecified so far
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, default to Release")
  set(CMAKE_BUILD_TYPE "Release")
endif()

# add test program if built standalone
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set(STANDALONE YES)
endif()

# those flags (from original Makefile)
# are crucial to make the program run properly
# otherwise you end up in an endless loop...
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-automatic -fno-backslash")

# add test program if built standalone
if(STANDALONE)
  # enable all the many warnings of APLCON
  set(CMAKE_Fortran_FLAGS_DEBUG  "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -Wall -Wextra")
  # but suppress warnings in release mode (just not to confuse people)
  set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -w")
  string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE)
  set(DEFAULT_Fortran_COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_${BUILD_TYPE}})
  message(STATUS "*** Build Type:      " ${CMAKE_BUILD_TYPE})
  message(STATUS "*** Compiler Flags:  " ${DEFAULT_Fortran_COMPILE_FLAGS})
endif()

# build the library
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

add_library(aplcon SHARED ${SRCS_LIB})

# build APLCON's test prgram only in standalone
if(STANDALONE)
  # link the test program against the Fortran library
  set(SRCS_TEST
    test/avlass.F
    test/avvsys.F
    test/avmain.F
    )
  add_executable(test_avmain ${SRCS_TEST})
  target_link_libraries(test_avmain aplcon)
endif()
# really no optimization in debug mode
if(CMAKE_CXX_COMPILER_ID STREQUAL GNU AND
    CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -Wall -Wextra")
else()
  message(FATAL_ERROR "Only GNU CXX/Fortran compiler supported")
endif()

# set default build type if unspecified so far
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, default to Release")
  set(CMAKE_BUILD_TYPE "Release")
endif()

# figure out compile flags
string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE)
set(DEFAULT_CXX_COMPILE_FLAGS ${CMAKE_CXX_FLAGS_${BUILD_TYPE}})
set(DEFAULT_Fortran_COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_${BUILD_TYPE}})

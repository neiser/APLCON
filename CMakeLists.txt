cmake_minimum_required (VERSION 2.6)
project(kinfit CXX Fortran)

add_subdirectory(APLCON)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include_directories(.)

add_library(aplcon_wrapper APLCON_wrapper.f90)
target_link_libraries(aplcon_wrapper aplcon)

add_executable(APLCON_test APLCON_test.cc)
target_link_libraries(APLCON_test aplcon_wrapper)

add_library(kinfit KinFit.cc KinFit.h)
target_link_libraries(kinfit aplcon_wrapper)

add_executable(kinfit_test KinFit_test.cc)
target_link_libraries(kinfit_test kinfit)

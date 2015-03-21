cmake_minimum_required (VERSION 2.6)
project(kinfit CXX Fortran)

add_subdirectory(APLCON)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include_directories(.)
add_executable(APLCON_test APLCON_test.cpp APLCON_wrapper.f90)
target_link_libraries(APLCON_test aplcon)

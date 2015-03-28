cmake_minimum_required (VERSION 2.6)
project(APLCON++ CXX Fortran)

# check for in-source build, forbid it!
if(CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR AND NOT MSVC_IDE)
  message(FATAL_ERROR "\nIn-source build attempt detected!\n"
    "Please create a new directory (e.g. build) and run `cmake ..`.\n"
    "Also don't forget to delete the created CMakeCache.txt and CMakeFiles dir"
    )
endif()

include (cmake/settings.cmake)

message(STATUS "*** Build Type:      " ${CMAKE_BUILD_TYPE})
message(STATUS "*** Compiler Flags:  " ${DEFAULT_COMPILE_FLAGS})
message(STATUS "*** Install libs to: " ${LIBRARY_OUTPUT_PATH})
message(STATUS "*** Install bin to:  " ${EXECUTABLE_OUTPUT_PATH})

add_subdirectory(APLCON)

include_directories(src)

add_library(aplcon++ SHARED
  src/APLCON.cc
  src/wrapper/APLCON.f90
  src/wrapper/APLCON.h
  # add header files to show them in IDEs
  src/APLCON.hpp
  src/detail/APLCON_hpp.hpp
  src/detail/APLCON_cc.hpp
  )
target_link_libraries(aplcon++ aplcon)

# build some examples
add_executable(APLCON_example_01 src/example/01_simple.cc)
target_link_libraries(APLCON_example_01 aplcon++)

add_executable(APLCON_example_02 src/example/02_linker.cc)
target_link_libraries(APLCON_example_02 aplcon++)

add_executable(APLCON_example_03 src/example/03_kinfit.cc)
target_link_libraries(APLCON_example_03 aplcon++)

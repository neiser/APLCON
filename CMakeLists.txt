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

add_library(aplcon++ SHARED src/APLCON.cc src/wrapper/APLCON.f90)
target_link_libraries(aplcon++ aplcon)

add_executable(aplcon_example src/APLCON_example.cc)
target_link_libraries(aplcon_example aplcon++)

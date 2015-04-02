# look for Doxygen
find_package(Doxygen)
if(NOT DOXYGEN_FOUND)
  return()
endif()

set(DOXYGEN_INPUT_DIR "${CMAKE_SOURCE_DIR}")

set(DOXYGEN_OUTPUT_DIR ${CMAKE_SOURCE_DIR}/doc)

# we look for Graphviz, if found, we enable the
# include graph generation
find_program(DOXYGEN_DOT_EXE dot)
if(DOXYGEN_DOT_EXE)
  set(DOXYGEN_DOT_FOUND "YES")
else()
  message(STATUS "Install Graphviz package to enable include graph generation")
  set(DOXYGEN_DOT_FOUND "NO")
endif()

# we use one global Doxyfile
configure_file(${CMAKE_SOURCE_DIR}/cmake/Doxyfile.in
  ${CMAKE_BINARY_DIR}/Doxyfile @ONLY)

# and define how to run Doxygen
# it always outputs a symbolic doc dir to avoid clashing
# with the target "doc"
add_custom_command(OUTPUT "DoxyfileDocDir"
  COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMENT "Running Doxygen...")
set_source_files_properties("DoxyfileDocDir" PROPERTIES SYMBOLIC "on")

# so we can make it saying "make doc"
add_custom_target(doc DEPENDS "DoxyfileDocDir")

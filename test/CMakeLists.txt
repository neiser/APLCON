# the catch++ library is the workhorse for our tests
include_directories(.)
add_custom_target(build_and_test
  COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -j${CTEST_PARALLEL_JOBS}
  COMMENT "Running build_and_test with ${CTEST_PARALLEL_JOBS} workers"
  )
add_library(catch EXCLUDE_FROM_ALL catch.cc catch.hpp)
target_link_libraries(catch)

macro(add_aplcon_test name)
  set(TESTNAME "test_${name}")
  set(TESTFILE "Test${name}.cc")
  add_executable(${TESTNAME} EXCLUDE_FROM_ALL ${TESTFILE})
  target_link_libraries(${TESTNAME} catch aplcon++)
  add_test(NAME ${TESTNAME} COMMAND ${TESTNAME})
  set_tests_properties(${TESTNAME} PROPERTIES TIMEOUT 30)
  add_dependencies(build_and_test ${TESTNAME})
endmacro()

include_directories(${CMAKE_SOURCE_DIR}/src)

add_aplcon_test(VerySimple)
add_aplcon_test(Simple)
add_aplcon_test(Linker)


add_executable(missoula_cpucuda
	missoula_user.cpp
	missoula_user.h
	missoula.cpp
)

target_include_directories(missoula_cpucuda PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
)

target_link_libraries(missoula_cpucuda PRIVATE 
  GEOFLOOD::CPU_CUDA
)

# -- set environment string used for regression tests
set(FCLAW_TEST_ENVIRONMENT "FCLAW_EXAMPLES_BUILD_DIR=${PROJECT_BINARY_DIR};FCLAW_EXAMPLES_SRC_DIR=${PROJECT_SOURCE_DIR}")

add_test(NAME geoflood_missoula_cpu COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(geoflood_missoula_cpu PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")

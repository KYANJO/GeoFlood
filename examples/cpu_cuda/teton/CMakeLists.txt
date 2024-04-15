add_library(teton_cpu_f OBJECT
	reservoir_module.f90
	qinit.f90
)

target_link_libraries(teton_cpu_f PRIVATE 
  GEOFLOOD::CPU_CUDA
)

add_executable(teton_cpu_cuda 
	teton_user.cpp
	teton_user.h
	teton.cpp
)

target_include_directories(teton_cpu_cuda PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
)

target_link_libraries(teton_cpu_cuda PRIVATE 
  teton_cpu_f
  GEOFLOOD::CPU_CUDA
)

# -- set environment string used for regression tests
set(FCLAW_TEST_ENVIRONMENT "FCLAW_EXAMPLES_BUILD_DIR=${PROJECT_BINARY_DIR};FCLAW_EXAMPLES_SRC_DIR=${PROJECT_SOURCE_DIR}")

add_test(NAME geoflood_teton_cpu COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(geoflood_teton_cpu PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")

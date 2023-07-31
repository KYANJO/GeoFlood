add_library(teton_cpu_f OBJECT
	reservoir_module.f90
	qinit.f90
)

add_library(teton_cuda_f OBJECT
  user_cuda/qinit.f90
)

add_executable(teton_cuda 
	teton_user.cpp
	teton_user.h
	teton.cpp
  teton_options.c
	$<TARGET_OBJECTS:teton_cuda_f>
  $<TARGET_OBJECTS:teton_cpu_f>
)

target_include_directories(teton_cuda PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
)

target_link_libraries(teton_cuda PRIVATE 
  GEOFLOOD::CUDACLAW
)

# -- set environment string used for regression tests
set(FCLAW_TEST_ENVIRONMENT "FCLAW_EXAMPLES_BUILD_DIR=${PROJECT_BINARY_DIR};FCLAW_EXAMPLES_SRC_DIR=${PROJECT_SOURCE_DIR}")

add_test(NAME geoflood_teton_cuda COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(geoflood_teton_cuda PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")

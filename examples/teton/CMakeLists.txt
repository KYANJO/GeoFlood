add_library(teton_f OBJECT
	reservoir_module.f90
	qinit.f90
)

target_link_libraries(teton_f PRIVATE 
  GEOFLOOD::GEOCLAW
)

add_executable(teton 
	teton_user.cpp
	teton_user.h
	teton.cpp
)

target_include_directories(teton PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
)

target_link_libraries(teton PRIVATE 
  teton_f
  GEOFLOOD::GEOCLAW
)

# -- set environment string used for regression tests
set(FCLAW_TEST_ENVIRONMENT "FCLAW_EXAMPLES_BUILD_DIR=${PROJECT_BINARY_DIR};FCLAW_EXAMPLES_SRC_DIR=${PROJECT_SOURCE_DIR}")

add_test(NAME geoflood_teton COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(geoflood_teton PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")

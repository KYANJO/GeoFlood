# Add boundary conditions as a library
add_library(dam_break_f OBJECT
  # fortran/fc2d_geoflood_bc2.f90
  ../../../../src/fc2d_geoclaw/hydrograph/hydrograh_module.f90
  ../../../../src/fc2d_geoclaw/geolib_source/geoclaw_module.f90 
  ../../../../src/fc2d_geoclaw/geolib_source/utility_module.f90
  ../../../../src/fc2d_geoclaw/amrlib_source/amr_module.f90
  qinit.f90
)

add_executable(dam_break 
	dam_break_user.cpp
	dam_break_user.h
	dam_break.cpp
	$<TARGET_OBJECTS:dam_break_f>
)

target_include_directories(dam_break PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
)

target_link_libraries(dam_break PRIVATE 
  GEOFLOOD::GEOCLAW
)

# -- set environment string used for regression tests
set(FCLAW_TEST_ENVIRONMENT "FCLAW_EXAMPLES_BUILD_DIR=${PROJECT_BINARY_DIR};FCLAW_EXAMPLES_SRC_DIR=${PROJECT_SOURCE_DIR}")

add_test(NAME geoflood_dam_break_cpu COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(geoflood_dam_break_cpu PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")

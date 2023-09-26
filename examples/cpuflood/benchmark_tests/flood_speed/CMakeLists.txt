# Add boundary conditions as a library
add_library(flood_speed_f OBJECT
  qinit.f90
  fortran/fc2d_geoflood_bc2.f90
  ../../../../src/fc2d_geoclaw/hydrograph/hydrograh_module.f90
  ../../../../src/fc2d_geoclaw/geolib_source/geoclaw_module.f90 
  ../../../../src/fc2d_geoclaw/geolib_source/utility_module.f90
  ../../../../src/fc2d_geoclaw/amrlib_source/amr_module.f90
  ../../../../src/fc2d_geoclaw/fortran_source/extract_dt.f90
)

add_executable(flood_speed 
	flood_speed_user.cpp
	flood_speed_user.h
	flood_speed.cpp
	$<TARGET_OBJECTS:flood_speed_f>
)

target_include_directories(flood_speed PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
)

target_link_libraries(flood_speed PRIVATE 
  GEOFLOOD::GEOCLAW
)

# -- set environment string used for regression tests
set(FCLAW_TEST_ENVIRONMENT "FCLAW_EXAMPLES_BUILD_DIR=${PROJECT_BINARY_DIR};FCLAW_EXAMPLES_SRC_DIR=${PROJECT_SOURCE_DIR}")

add_test(NAME geoflood_flood_speed_cpu COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(geoflood_flood_speed_cpu PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")

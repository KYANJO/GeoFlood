# Add boundary conditions as a library
add_library(valley_flood_f OBJECT
  fortran/fc2d_geoflood_bc2.f90
  ../../../../src/fc2d_geoclaw/hydrograph/hydrograh_module.f90
  ../../../../src/fc2d_geoclaw/geolib_source/geoclaw_module.f90 
  ../../../../src/fc2d_geoclaw/geolib_source/utility_module.f90
  ../../../../src/fc2d_geoclaw/amrlib_source/amr_module.f90
)

add_executable(valley_flood 
	valley_flood_user.cpp
	valley_flood_user.h
	valley_flood.cpp
	$<TARGET_OBJECTS:valley_flood_f>
)

target_include_directories(valley_flood PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
)

target_link_libraries(valley_flood PRIVATE 
  GEOFLOOD::GEOCLAW
)

# -- set environment string used for regression tests
set(FCLAW_TEST_ENVIRONMENT "FCLAW_EXAMPLES_BUILD_DIR=${PROJECT_BINARY_DIR};FCLAW_EXAMPLES_SRC_DIR=${PROJECT_SOURCE_DIR}")

add_test(NAME geoflood_valley_flood_cpu COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(geoflood_valley_flood_cpu PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")


# Add boundary conditions as a library
add_library(disconnected_water_body_f OBJECT
  fortran/fc2d_geoflood_bc2.f90
  ../../../../src/fc2d_geoclaw/hydrograph/hydrograh_module.f90
  ../../../../src/fc2d_geoclaw/geolib_source/geoclaw_module.f90 
  ../../../../src/fc2d_geoclaw/geolib_source/utility_module.f90
  ../../../../src/fc2d_geoclaw/amrlib_source/amr_module.f90
)

add_executable(disconnected_water_body 
	disconnected_water_body_user.cpp
	disconnected_water_body_user.h
	disconnected_water_body.cpp
	$<TARGET_OBJECTS:disconnected_water_body_f>
)

target_include_directories(disconnected_water_body PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_BINARY_DIR}/src/fc2d_geoclaw/include
)

target_link_libraries(disconnected_water_body PRIVATE 
  GEOFLOOD::GEOCLAW
)

# -- set environment string used for regression tests
set(FCLAW_TEST_ENVIRONMENT "FCLAW_EXAMPLES_BUILD_DIR=${PROJECT_BINARY_DIR};FCLAW_EXAMPLES_SRC_DIR=${PROJECT_SOURCE_DIR}")

add_test(NAME geoflood_disconnected_water_body_cpu COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(geoflood_disconnected_water_body_cpu PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")

# Add boundary conditions as a library
add_library(filling_depressions_cpu_f OBJECT
  # fortran/fc2d_geoflood_bc2.f90
  ../../../../src/fc2d_geoclaw/hydrograph/hydrograh_module.f90
  ../../../../src/fc2d_geoclaw/geolib_source/geoclaw_module.f90 
  ../../../../src/fc2d_geoclaw/geolib_source/utility_module.f90
  ../../../../src/fc2d_geoclaw/amrlib_source/amr_module.f90
)

add_executable(filling_depressions_cpucuda
	filling_depressions_user.cpp
	filling_depressions_user.h
	filling_depressions.cpp
	$<TARGET_OBJECTS:filling_depressions_cpu_f>
)

target_include_directories(filling_depressions_cpucuda PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
)

target_link_libraries(filling_depressions_cpucuda PRIVATE 
  GEOFLOOD::CPU_CUDA
)

# -- set environment string used for regression tests
set(FCLAW_TEST_ENVIRONMENT "FCLAW_EXAMPLES_BUILD_DIR=${PROJECT_BINARY_DIR};FCLAW_EXAMPLES_SRC_DIR=${PROJECT_SOURCE_DIR}")

add_test(NAME geoflood_filling_depressions_cpu COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(geoflood_filling_depressions_cpu PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")

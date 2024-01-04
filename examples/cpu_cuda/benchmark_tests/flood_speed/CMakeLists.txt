# Add boundary conditions as a library
add_library(flood_speed_cpu_f OBJECT
  # fortran/fc2d_geoflood_bc2.f90
  fortran/flood_speed_bc.f90
  ../../../../src/fc2d_cudaflood/hydrograph/hydrograh_module.f90
  ../../../../src/fc2d_cudaflood/geolib_source/geoclaw_module.f90 
  ../../../../src/fc2d_cudaflood/geolib_source/refinement_module.f90
  ../../../../src/fc2d_cudaflood/geolib_source/utility_module.f90
  ../../../../src/fc2d_cudaflood/geolib_source/topo_module.f90
  ../../../../src/fc2d_cudaflood/amrlib_source/amr_module.f90
  ../../../../src/fc2d_cudaflood/geolib_source/friction_module.f90
  ../../../../src/fc2d_cudaflood/fortran_source/extract_topo_variables.f90
  qinit.f90
)

add_executable(flood_speed_cpucuda 
$<TARGET_OBJECTS:flood_speed_cpu_f>
	flood_speed_user.cpp
	flood_speed_options.c
	flood_speed.cpp
)

target_include_directories(flood_speed_cpucuda PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
)

target_link_libraries(flood_speed_cpucuda PRIVATE 
  GEOFLOOD::CUDA_FLOOD
  # GEOFLOOD::GEOCLAW
  # GEOFLOOD::CPU_CUDA
)

# -- set environment string used for regression tests
set(FCLAW_TEST_ENVIRONMENT "FCLAW_EXAMPLES_BUILD_DIR=${PROJECT_BINARY_DIR};FCLAW_EXAMPLES_SRC_DIR=${PROJECT_SOURCE_DIR}")

add_test(NAME geoflood_flood_speed_cpu COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(geoflood_flood_speed_cpu PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")

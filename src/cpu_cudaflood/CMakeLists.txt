# -- fortran library
add_library(cpu_cuda_f OBJECT
	fortran_source/cudaclaw_inlinelimiter.f
	fortran_source/cudaclaw_set_capacity.f
	fortran_source/cudaclaw_block.f
	fortran_source/cudaclaw_bc2_default.f
	fortran_source/cudaclaw_flux2.f
	fortran_source/cudaclaw_flux2fw.f
	fortran_source/cudaclaw_step2.f
	fortran_source/cudaclaw_step2_wrap.f
	geoclaw_solver/amrlib_source/amr_module.f90 
	geoclaw_solver/geolib_source/utility_module.f90 
	geoclaw_solver/geolib_source/geoclaw_module.f90 
	geoclaw_solver/amrlib_source/regions_module.f90 
	geoclaw_solver/geolib_source/topo_module.f90 
	geoclaw_solver/geolib_source/qinit_module.f90 
	geoclaw_solver/geolib_source/fixedgrids_module.f90 
	geoclaw_solver/geolib_source/refinement_module.f90 
	geoclaw_solver/geolib_source/friction_module.f90 
	geoclaw_solver/hydrograph/hydrograh_module.f90
	geoclaw_solver/riemann_source/rpn2_geoclaw.f 
	geoclaw_solver/riemann_source/rpt2_geoclaw.f 
	geoclaw_solver/riemann_source/geoclaw_riemann_utils.f 
	geoclaw_solver/amrlib_source/opendatafile.f 
	geoclaw_solver/amrlib_source/check4nans.f90 
	geoclaw_solver/amrlib_source/quick_sort1.f 
	geoclaw_solver/geolib_source/topointegral.f 
	geoclaw_solver/geolib_source/topo_update.f90 
	geoclaw_solver/geolib_source/cellgridintegrate2.f 
	geoclaw_solver/geolib_source/allowflag.f
	geoclaw_solver/geolib_source/bilinearintegral.f 
	geoclaw_solver/fortran_source/philim.f 
	geoclaw_solver/fortran_source/geoclaw_set_modules_fort.f90 
	geoclaw_solver/fortran_source/geoclaw_step2_wrap_fort.f 
	geoclaw_solver/fortran_source/geoclaw_flux2fw_fort.f 
	geoclaw_solver/fortran_source/geoclaw_limiter_fort.f 
	geoclaw_solver/fortran_source/geoclaw_block_fort.f 
	geoclaw_solver/fortran_source/geoclaw_gauges_fort.f90  
	geoclaw_solver/fortran_source/geoclaw_topo_update_fort.f 
	geoclaw_solver/fortran_source/geoflood_bc2_fort.f90
	geoclaw_solver/fortran_source/geoclaw_step2_fort.f90 
	geoclaw_solver/fortran_source/geoclaw_setaux_fort.f90 
	geoclaw_solver/fortran_source/geoclaw_src2_fort.f90 
	geoclaw_solver/fortran_source/extract_dt.f90
	geoclaw_solver/fortran_source/geoclaw_b4step2_fort.f90 
	geoclaw_solver/fortran_source/geoclaw_qinit_fort.f90 
	geoclaw_solver/fclaw2d_source/fc2d_geoclaw_copy_fort.f 
	geoclaw_solver/fclaw2d_source/fc2d_geoclaw_average_fort.f 
	geoclaw_solver/fclaw2d_source/fc2d_geoclaw_interpolate_fort.f 
	geoclaw_solver/fclaw2d_source/fc2d_geoclaw_output_ascii_fort.f 
	geoclaw_solver/fclaw2d_source/fc2d_geoclaw_tag4refinement_fort.f90 
	geoclaw_solver/fclaw2d_source/fc2d_geoclaw_tag4coarsening_fort.f90 
	geoclaw_solver/fclaw2d_source/fc2d_geoclaw_test_regions.f90 
	geoclaw_solver/fclaw2d_source/fc2d_geoclaw_flag2refine_fort.f90 
	geoclaw_solver/fclaw2d_source/fc2d_geoclaw_local_ghost_pack_fort.f 
	geoclaw_solver/fclaw2d_source/fc2d_geoclaw_local_ghost_pack_aux_fort.f 
	geoclaw_solver/fclaw2d_source/fc2d_geoclaw_diagnostics_fort.f 
	geoclaw_solver/fclaw2d_source/fc2d_geoclaw_timeinterp_fort.f
)

set_target_properties(cpu_cuda_f PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

# -- add library 
add_library(cpu_cuda STATIC)
set_target_properties(cpu_cuda PROPERTIES EXPORT_NAME CPU_CUDA)

# device code needs to relcoateable since it is in a shared library
#target_compile_options(cudaclaw PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:
#                       --relocatable-device-code=true
#                       >)
# needed because we have device code in a library
set_target_properties(cpu_cuda PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

target_sources(cpu_cuda PRIVATE
  fc2d_cudaclaw.cpp
	fc2d_cudaclaw_options.c
	cuda_source/cudaclaw_step2.cu
	cuda_source/cudaclaw_initialize.cu
	cuda_source/cudaclaw_limiters.cu
	cuda_source/cudaclaw_allocate.cu
	cuda_source/cudaclaw_flux2.cu
	fc2d_cuda_profiler.cu
	cuda_source/cudaclaw_store_patches.cu
	cuda_source/data_swap.c
	geoclaw_solver/fc2d_geoclaw_gauges_default.c 
	geoclaw_solver/fc2d_geoclaw_output_ascii.c
	geoclaw_solver/fc2d_geoclaw_run.c
  $<TARGET_OBJECTS:cpu_cuda_f>
)

target_link_libraries(cpu_cuda PUBLIC 
	FORESTCLAW::FORESTCLAW 
	FORESTCLAW::CLAWPATCH 
	CUDA::nvToolsExt
)

target_include_directories(cpu_cuda
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)

if(CUDAToolkit_VERSION_MAJOR LESS 11)
	target_include_directories(cpu_cuda PRIVATE ${PROJECT_SOURCE_DIR}/cub)
endif()

# -- install
install(FILES
	fc2d_cudaclaw.h
	fc2d_cudaclaw_cuda.h
	fc2d_cudaclaw_check.h
	fc2d_cuda_profiler.h
	fc2d_cudaclaw_options.h
	cudaclaw_user_fort.h
	fc2d_cudaclaw_fort.h
	cuda_source/data_swap.h
	geoclaw_solver/types.h
	geoclaw_solver/fc2d_geoclaw.h
	geoclaw_solver/fc2d_geoclaw_fort.h
	geoclaw_solver/fc2d_geoclaw_gauges_default.h
  DESTINATION include
)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/include/geoclaw_module.mod
  # ${CMAKE_CURRENT_BINARY_DIR}/include/hydrograph_module.mod
  DESTINATION include
)

install(TARGETS cpu_cuda
  EXPORT ${PROJECT_NAME}-targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
	RUNTIME DESTINATION bin)

# imported target, for use from FetchContent
add_library(GEOFLOOD::CPU_CUDA INTERFACE IMPORTED GLOBAL)
target_link_libraries(GEOFLOOD::CPU_CUDA INTERFACE cpu_cuda)

if(BUILD_TESTING)
  add_executable(fc2d_cudaclaw.TEST
    fc2d_cudaclaw.h.TEST.cpp
    fc2d_cudaclaw_options.h.TEST.cpp
  )
  target_link_libraries(fc2d_cudaclaw.TEST testutils cpu_cuda)
#   register_unit_tests(fc2d_cudaclaw.TEST)
endif(BUILD_TESTING)

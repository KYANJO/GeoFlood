# -- fortran library
add_library(cpu_cuda_f OBJECT
    ../fc2d_geoclaw/amrlib_source/amr_module.f90 
    ../fc2d_geoclaw/geolib_source/utility_module.f90 
    ../fc2d_geoclaw/geolib_source/geoclaw_module.f90 
    ../fc2d_geoclaw/amrlib_source/regions_module.f90 
    ../fc2d_geoclaw/geolib_source/topo_module.f90 
    ../fc2d_geoclaw/geolib_source/qinit_module.f90 
    ../fc2d_geoclaw/geolib_source/fixedgrids_module.f90 
    ../fc2d_geoclaw/geolib_source/refinement_module.f90 
    ../fc2d_geoclaw/geolib_source/friction_module.f90 
    ../fc2d_geoclaw/hydrograph/hydrograh_module.f90
    ../fc2d_geoclaw/riemann_source/rpn2_geoclaw.f 
    ../fc2d_geoclaw/riemann_source/rpt2_geoclaw.f 
    ../fc2d_geoclaw/riemann_source/geoclaw_riemann_utils.f 
    ../fc2d_geoclaw/amrlib_source/opendatafile.f 
    ../fc2d_geoclaw/amrlib_source/check4nans.f90 
    ../fc2d_geoclaw/amrlib_source/quick_sort1.f 
    ../fc2d_geoclaw/geolib_source/topointegral.f 
    ../fc2d_geoclaw/geolib_source/topo_update.f90 
    ../fc2d_geoclaw/geolib_source/cellgridintegrate2.f 
    ../fc2d_geoclaw/geolib_source/allowflag.f
    ../fc2d_geoclaw/geolib_source/bilinearintegral.f 
    ../fc2d_geoclaw/fortran_source/philim.f 
    ../fc2d_geoclaw/fortran_source/geoclaw_set_modules_fort.f90 
    ../fc2d_geoclaw/fortran_source/geoclaw_step2_wrap_fort.f 
    ../fc2d_geoclaw/fortran_source/geoclaw_flux2fw_fort.f 
    ../fc2d_geoclaw/fortran_source/geoclaw_limiter_fort.f 
    ../fc2d_geoclaw/fortran_source/geoclaw_block_fort.f 
    ../fc2d_geoclaw/fortran_source/geoclaw_gauges_fort.f90  
    ../fc2d_geoclaw/fortran_source/geoclaw_topo_update_fort.f 
    ../fc2d_geoclaw/fortran_source/extract_geoclaw_parameters.f90
    ../fc2d_geoclaw/fortran_source/geoflood_bc2_fort.f90
    ../fc2d_geoclaw/fortran_source/geoclaw_step2_fort.f90 
    ../fc2d_geoclaw/fortran_source/geoclaw_setaux_fort.f90 
    ../fc2d_geoclaw/fortran_source/geoclaw_src2_fort.f90 
    ../fc2d_geoclaw/fortran_source/extract_dt.f90
    ../fc2d_geoclaw/fortran_source/topointegral_geo.f
    ../fc2d_geoclaw/fortran_source/geoclaw_b4step2_fort.f90 
    ../fc2d_geoclaw/fortran_source/geoclaw_qinit_fort.f90 
    ../fc2d_geoclaw/fclaw2d_source/fc2d_geoclaw_copy_fort.f 
    ../fc2d_geoclaw/fclaw2d_source/fc2d_geoclaw_average_fort.f 
    ../fc2d_geoclaw/fclaw2d_source/fc2d_geoclaw_interpolate_fort.f 
    ../fc2d_geoclaw/fclaw2d_source/fc2d_geoclaw_output_ascii_fort.f 
    ../fc2d_geoclaw/fclaw2d_source/fc2d_geoclaw_tag4refinement_fort.f90 
    ../fc2d_geoclaw/fclaw2d_source/fc2d_geoclaw_tag4coarsening_fort.f90 
    ../fc2d_geoclaw/fclaw2d_source/fc2d_geoclaw_test_regions.f90 
    # ../fc2d_geoclaw/fclaw2d_source/fc2d_geoclaw_flag2refine_fort.f90 
    ../fc2d_geoclaw/fclaw2d_source/fc2d_geoflood_flag2refine.f90
    ../fc2d_geoclaw/fclaw2d_source/fc2d_geoclaw_local_ghost_pack_fort.f 
    ../fc2d_geoclaw/fclaw2d_source/fc2d_geoclaw_local_ghost_pack_aux_fort.f 
    ../fc2d_geoclaw/fclaw2d_source/fc2d_geoclaw_diagnostics_fort.f 
    ../fc2d_geoclaw/fclaw2d_source/fc2d_geoclaw_timeinterp_fort.f
    cuda_source/data_swap.f90
)

target_link_Libraries(cpu_cuda_f PRIVATE 
    FORESTCLAW::CLAWPATCH
)

set_target_properties(cpu_cuda_f PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

# -- add library
add_library(cpu_cuda_c)
set_target_properties(cpu_cuda_c PROPERTIES EXPORT_NAME CPU_CUDA)

# needed because we have device code in a library
set_target_properties(cpu_cuda_c PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

target_sources(cpu_cuda_c PRIVATE
    $<TARGET_OBJECTS:cpu_cuda_f>
    fc2d_cpucuda.cpp 
    fc2d_cpucuda_options.c 
    geoflood_options_user.c
    ../fc2d_geoclaw/fc2d_geoclaw_gauges_default.c 
    ../fc2d_geoclaw/fc2d_geoclaw_run.c 
    ../fc2d_geoclaw/fc2d_geoclaw_output_ascii.c
    cuda_source/cudaclaw_step2.cu
    cuda_source/cudaclaw_initialize.cu
    cuda_source/cudaclaw_limiters.cu
    cuda_source/cudaclaw_allocate.cu
    cuda_source/cudaclaw_flux2.cu
    fc2d_cuda_profiler.cu
    cuda_source/cudaclaw_store_patches.cu
    riemann_source/variables.cu
    riemann_source/geoflood_src2.cu
    riemann_source/geoflood_rp2.cu
)

target_link_libraries(cpu_cuda_c PUBLIC FORESTCLAW::FORESTCLAW FORESTCLAW::CLAWPATCH CUDA::nvToolsExt)

target_include_directories(cpu_cuda_c PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
)

if(CUDAToolkit_VERSION_MAJOR LESS 11)
	target_include_directories(cpu_cuda_c PRIVATE ${PROJECT_SOURCE_DIR}/cub)
endif()

# -- install
install(FILES
  fc2d_cpucuda.h
  types.h
  fc2d_cpucuda_options.h
  geoflood_options_user.h
  fc2d_cpucuda_fort.h
  ../fc2d_geoclaw/fc2d_geoclaw_gauges_default.h
  ../fc2d_geoclaw/fc2d_geoclaw_output_ascii.h
  fc2d_cudaclaw_cuda.h
  fc2d_cudaclaw_check.h
	fc2d_cuda_profiler.h
  DESTINATION include
)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/include/geoclaw_module.mod
  ${CMAKE_CURRENT_BINARY_DIR}/include/refinement_module.mod
  ${CMAKE_CURRENT_BINARY_DIR}/include/hydrograph_module.mod
  ${CMAKE_CURRENT_BINARY_DIR}/include/topo_module.mod
  DESTINATION include
)

install(TARGETS cpu_cuda_c
  EXPORT ${PROJECT_NAME}Targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin)

# imported target, for use from FetchContent
add_library(GEOFLOOD::CPU_CUDA INTERFACE IMPORTED GLOBAL)
target_link_libraries(GEOFLOOD::CPU_CUDA INTERFACE cpu_cuda_c)
# add_library(GEOFLOOD::CUDAFOOD ALIAS cpu_cuda_c)

if(TARGET GEOFLOOD::CPU_CUDA)
    # If the target exists, print a message
    message(STATUS "Target GEOFLOOD::CPU_CUDA exists.")
else()
    # If the target does not exist, print a different message
    message(WARNING "Target GEOFLOOD::CPU_CUDA does not exist.")
endif()


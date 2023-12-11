# -- fortran library
add_library(cudaflood_f OBJECT
    amrlib_source/amr_module.f90 
    geolib_source/utility_module.f90 
    geolib_source/geoclaw_module.f90 
    amrlib_source/regions_module.f90 
    geolib_source/topo_module.f90 
    geolib_source/qinit_module.f90 
    geolib_source/fixedgrids_module.f90 
    geolib_source/refinement_module.f90 
    geolib_source/friction_module.f90 
    hydrograph/hydrograh_module.f90
    riemann_source/rpn2_geoclaw.f 
    riemann_source/rpt2_geoclaw.f 
    riemann_source/geoclaw_riemann_utils.f 
    amrlib_source/opendatafile.f 
    amrlib_source/check4nans.f90 
    amrlib_source/quick_sort1.f 
    geolib_source/topointegral.f 
    geolib_source/topo_update.f90 
    geolib_source/cellgridintegrate2.f 
    geolib_source/allowflag.f
    geolib_source/bilinearintegral.f 
    fortran_source/philim.f 
    fortran_source/geoclaw_set_modules_fort.f90 
    fortran_source/geoclaw_step2_wrap_fort.f 
    fortran_source/geoclaw_flux2fw_fort.f 
    fortran_source/geoclaw_limiter_fort.f 
    fortran_source/geoclaw_block_fort.f 
    fortran_source/geoclaw_gauges_fort.f90  
    fortran_source/geoclaw_topo_update_fort.f 
    # fortran_source/geoclaw_bc2_fort.f 
    fortran_source/geoflood_bc2_fort.f90
    fortran_source/geoclaw_step2_fort.f90 
    fortran_source/geoclaw_setaux_fort.f90 
    fortran_source/geoclaw_src2_fort.f90 
    fortran_source/extract_dt.f90
    fortran_source/geoclaw_b4step2_fort.f90 
    fortran_source/geoclaw_qinit_fort.f90 
    fclaw2d_source/fc2d_geoclaw_copy_fort.f 
    fclaw2d_source/fc2d_geoclaw_average_fort.f 
    fclaw2d_source/fc2d_geoclaw_interpolate_fort.f 
    fclaw2d_source/fc2d_geoclaw_output_ascii_fort.f 
    fclaw2d_source/fc2d_geoclaw_tag4refinement_fort.f90 
    fclaw2d_source/fc2d_geoclaw_tag4coarsening_fort.f90 
    fclaw2d_source/fc2d_geoclaw_test_regions.f90 
    fclaw2d_source/fc2d_geoclaw_flag2refine_fort.f90 
    fclaw2d_source/fc2d_geoclaw_local_ghost_pack_fort.f 
    fclaw2d_source/fc2d_geoclaw_local_ghost_pack_aux_fort.f 
    fclaw2d_source/fc2d_geoclaw_diagnostics_fort.f 
    fclaw2d_source/fc2d_geoclaw_timeinterp_fort.f
    cuda_source/data_swap.f90
)

target_link_Libraries(cudaflood_f PRIVATE 
    FORESTCLAW::CLAWPATCH
)

set_target_properties(cudaflood_f PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

# -- add library
add_library(cudaflood_c)
set_target_properties(cudaflood_c PROPERTIES EXPORT_NAME CUDA_FLOOD)

# needed because we have device code in a library
set_target_properties(cudaflood_c PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

target_sources(cudaflood_c PRIVATE
    $<TARGET_OBJECTS:cudaflood_f>
    fc2d_geoclaw.cpp 
    fc2d_geoclaw_options.c 
    fc2d_geoclaw_gauges_default.c 
    fc2d_geoclaw_run.c 
    fc2d_geoclaw_output_ascii.c
    cuda_source/cudaclaw_step2.cu
    cuda_source/cudaclaw_initialize.cu
    cuda_source/cudaclaw_limiters.cu
    cuda_source/cudaclaw_allocate.cu
    cuda_source/cudaclaw_flux2.cu
    fc2d_cuda_profiler.cu
    cuda_source/cudaclaw_store_patches.cu
)

target_link_libraries(cudaflood_c PUBLIC FORESTCLAW::FORESTCLAW FORESTCLAW::CLAWPATCH CUDA::nvToolsExt)

target_include_directories(cudaflood_c PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
)

if(CUDAToolkit_VERSION_MAJOR LESS 11)
	target_include_directories(cudaflood_c PRIVATE ${PROJECT_SOURCE_DIR}/cub)
endif()

# -- install
install(FILES
  fc2d_geoclaw.h
  types.h
  fc2d_geoclaw_options.h
  fc2d_geoclaw_fort.h
  fc2d_geoclaw_gauges_default.h
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

install(TARGETS cudaflood_c
  EXPORT ${PROJECT_NAME}Targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin)

# imported target, for use from FetchContent
add_library(GEOFLOOD::CUDA_FLOOD INTERFACE IMPORTED GLOBAL)
target_link_libraries(GEOFLOOD::CUDA_FLOOD INTERFACE cudaflood_c)
# add_library(GEOFLOOD::CUDAFOOD ALIAS cudaflood_c)

if(TARGET GEOFLOOD::CUDA_FLOOD)
    # If the target exists, print a message
    message(STATUS "Target GEOFLOOD::CUDA_FLOOD exists.")
else()
    # If the target does not exist, print a different message
    message(WARNING "Target GEOFLOOD::CUDA_FLOOD does not exist.")
endif()


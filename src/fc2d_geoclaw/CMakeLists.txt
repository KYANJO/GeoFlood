# -- fortran library
add_library(geoflood_f OBJECT
    amrlib_source/amr_module.f90 
    geolib_source/utility_module.f90 
    geolib_source/geoclaw_module.f90 
    amrlib_source/regions_module.f90 
    geolib_source/topo_module.f90 
    geolib_source/qinit_module.f90 
    geolib_source/fixedgrids_module.f90 
    geolib_source/refinement_data_module.f90
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
    # fclaw2d_source/fc2d_geoclaw_flag2refine_fort.f90
    fclaw2d_source/fc2d_geoflood_flag2refine.f90 
    fclaw2d_source/fc2d_geoclaw_local_ghost_pack_fort.f 
    fclaw2d_source/fc2d_geoclaw_local_ghost_pack_aux_fort.f 
    fclaw2d_source/fc2d_geoclaw_diagnostics_fort.f 
    fclaw2d_source/fc2d_geoclaw_timeinterp_fort.f
)

target_link_Libraries(geoflood_f PRIVATE 
    FORESTCLAW::CLAWPATCH
)

set_target_properties(geoflood_f PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

# -- add library
add_library(geoflood)
set_target_properties(geoflood PROPERTIES EXPORT_NAME GEOCLAW)

target_sources(geoflood PRIVATE 
    $<TARGET_OBJECTS:geoflood_f>
    fc2d_geoclaw.cpp 
    fc2d_geoclaw_options.c 
    fc2d_geoclaw_gauges_default.c 
    fc2d_geoclaw_run.c 
    fc2d_geoclaw_output_ascii.c
)

target_link_libraries(geoflood PUBLIC FORESTCLAW::FORESTCLAW FORESTCLAW::CLAWPATCH)

target_include_directories(geoflood PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
)

# -- install
install(FILES
  fc2d_geoclaw.h
  types.h
  fc2d_geoclaw_options.h
  fc2d_geoclaw_fort.h
  fc2d_geoclaw_gauges_default.h
  DESTINATION include
)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/include/geoclaw_module.mod
  ${CMAKE_CURRENT_BINARY_DIR}/include/refinement_data_module.mod
  ${CMAKE_CURRENT_BINARY_DIR}/include/hydrograph_module.mod
  DESTINATION include
)

install(TARGETS geoflood
  EXPORT ${PROJECT_NAME}Targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)


# imported target, for use from FetchContent
add_library(GEOFLOOD::GEOCLAW INTERFACE IMPORTED GLOBAL)
target_link_libraries(GEOFLOOD::GEOCLAW INTERFACE geoflood)
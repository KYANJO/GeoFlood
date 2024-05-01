# set path to modules
set(mod_path "../../../src/fc2d_geoclaw")

add_library(bathymetry_f OBJECT
	${mod_path}/geolib_source/geoclaw_module.f90
	${mod_path}/geolib_source/utility_module.f90
	${mod_path}/geolib_source/friction_module.f90
	${mod_path}/geolib_source/topo_module.f90
	${mod_path}/amrlib_source/amr_module.f90
	src/bilinearintegral.f
	src/topointegral.f
	src/geoclaw_setaux_fort.f90
)

add_executable(bathymetry 
	bathymetry_user.cpp
	bathymetry.cpp
	$<TARGET_OBJECTS:bathymetry_f>
)

target_include_directories(bathymetry PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
)

target_link_libraries(bathymetry PRIVATE 
  GEOFLOOD::GEOCLAW
)

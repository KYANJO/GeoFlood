add_library(bathymetry_f OBJECT
	reservoir_module.f90
	qinit.f90
)

add_executable(bathymetry 
	bathymetry_user.cpp
	bathymetry_user.h
	bathymetry.cpp
)

target_include_directories(bathymetry PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
)


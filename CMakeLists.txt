# --------------------------------------------
# @author:  Brian Kyanjo
# @contact: briankyanjo@u.boisestate.edu
# @date:    2022-10-16
# @version: 1.4
# 
# List of available applications and examples
#
# Note : Any files below which are uncommented
# need to be available at compile time, even 
# the configuration flag is not set.  
# --------------------------------------------


cmake_minimum_required(VERSION 3.13...3.20)

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES 70)
endif()

project(GEOFLOOD 
    DESCRIPTION "Flood inundation model"
    HOMEPAGE_URL https://github.com/KYANJO/GeoFlood
    LANGUAGES C CXX Fortran
    VERSION 1.4)

# --- ctest
include(CTest)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)

set(CMAKE_CXX_STANDARD 17)

# --- user options

include(cmake/options.cmake)
message(STATUS "${PROJECT_NAME} ${PROJECT_VERSION}  install prefix: ${CMAKE_INSTALL_PREFIX}")

# -- check if this being included as a subdirectory
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/ForestClawExternalProject)

# --- get installed ForestClaw package

# Set to True if you want to download ForestClaw and build from scratch
set(external_build False)
if(NOT external_build)
    # 'clawpack' here refers to the "-Dclawpack=on" cmake flag
    # This uses -DFORESTCLAW_ROOT to find ForestClaw installation
    find_package(FORESTCLAW REQUIRED COMPONENTS clawpack mpi)
else()
    # --- fall back to external project build
    include(ForestClawExternalProject)
    ForestClawExternalProject(TAG develop COMPONENTS clawpack mpi)
endif()
 
if(cudaclaw)
  enable_language(CUDA)
  set(CMAKE_CUDA_STANDARD 14)
  find_package(CUDAToolkit REQUIRED)
endif()

# --- configure files
include(cmake/config.cmake)

# --- src directory
if (cpuflood)
  add_subdirectory(src/fc2d_geoclaw)
endif(cpuflood)

# --- cudaclaw
if(cudaflood)
  enable_language(CUDA)
  set(CMAKE_CUDA_STANDARD 14)
  find_package(CUDAToolkit REQUIRED)
endif(cudaflood)
 

  if(cudaflood)
    add_subdirectory(src/fc2d_cudaclaw)
  endif(cudaflood)

# --- examples
add_subdirectory(examples)

# --- test
if(BUILD_TESTING)
  add_subdirectory(test)
endif(BUILD_TESTING)


# -- set environment string used for regression tests
set(FCLAW_TEST_ENVIRONMENT "FCLAW_APPLICATIONS_BUILD_DIR=${PROJECT_BINARY_DIR};FCLAW_APPLICATIONS_SRC_DIR=${PROJECT_SOURCE_DIR}")

if(TARGET MPI::MPI_C)
    list(APPEND FCLAW_TEST_ENVIRONMENT "FCLAW_MPIRUN=${MPIEXEC_EXECUTABLE}")
    list(APPEND FCLAW_TEST_ENVIRONMENT "FCLAW_MPI_TEST_FLAGS=${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS}")
endif()  



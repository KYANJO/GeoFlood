# --------------------------------------------
# List of available applications and examples
#
# Note : Any files below which are uncommented
# need to be available at compile time, even 
# the configuration flag is not set.  
# --------------------------------------------


cmake_minimum_required(VERSION 3.13...3.20)

project(GEOFLOOD LANGUAGES C CXX Fortran)

set(CMAKE_CXX_STANDARD 17)

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
 
add_subdirectory(src)
add_subdirectory(examples)

# -- set environment string used for regression tests
set(FCLAW_TEST_ENVIRONMENT "FCLAW_APPLICATIONS_BUILD_DIR=${PROJECT_BINARY_DIR};FCLAW_APPLICATIONS_SRC_DIR=${PROJECT_SOURCE_DIR}")
if(TARGET MPI::MPI_C)
    list(APPEND FCLAW_TEST_ENVIRONMENT "FCLAW_MPIRUN=${MPIEXEC_EXECUTABLE}")
    list(APPEND FCLAW_TEST_ENVIRONMENT "FCLAW_MPI_TEST_FLAGS=${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS}")
endif()  



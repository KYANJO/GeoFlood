# consumes ForestClaw as ExternalProject, provided imported target forestclaw::forestclaw
include(ExternalProject)

# target_link_libraries(... forestclaw::forestclaw)
# for user programs
add_library(FORESTCLAW::FORESTCLAW INTERFACE IMPORTED)
add_library(FORESTCLAW::CLAWPATCH INTERFACE IMPORTED)

if(NOT FCLAW_ROOT)
  if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(FCLAW_ROOT ${PROJECT_BINARY_DIR} CACHE PATH "default ForestClaw ROOT")
  else()
    set(FCLAW_ROOT ${CMAKE_INSTALL_PREFIX})
  endif()
endif()

find_package(MPI COMPONENTS C REQUIRED)
find_package(ZLIB)
find_package(SC)
find_package(P4EST)

# - forestclaw libraries

set(FORESTCLAW_LIBRARIES ${FCLAW_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}forestclaw${CMAKE_STATIC_LIBRARY_SUFFIX})

set(CLAWPATCH_LIBRARIES ${FCLAW_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}clawpatch${CMAKE_STATIC_LIBRARY_SUFFIX})

set(FORESTCLAW_INCLUDE_DIRS ${FCLAW_ROOT}/include)

# - depedencies

set(SC_LIBRARIES ${FCLAW_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}sc${CMAKE_STATIC_LIBRARY_SUFFIX})

set(P4EST_LIBRARIES ${FCLAW_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}p4est${CMAKE_STATIC_LIBRARY_SUFFIX})

# - determine byproducts

set(forestclaw_byproducts
${FORESTCLAW_LIBRARIES}
${CLAWPATCH_LIBRARIES}
)

if(NOT SC_FOUND)
  list(APPEND forestclaw_byproducts ${SC_LIBRARIES})
endif()

if(NOT P4EST_FOUND)
  list(APPEND forestclaw_byproducts ${P4EST_LIBRARIES})
endif()

# - cmake arguments for external project
set(forestclaw_args
-DCMAKE_INSTALL_PREFIX:PATH=${FCLAW_ROOT}
-DBUILD_SHARED_LIBS:BOOL=false
-DCMAKE_BUILD_TYPE=Release
-DBUILD_TESTING:BOOL=false
-Dclawpack:BOOL=true
-Dmpi:BOOL=true
)

# - external project
ExternalProject_Add(FORESTCLAW
GIT_REPOSITORY ${forestclaw_git}
GIT_TAG ${forestclaw_tag}
CMAKE_ARGS ${forestclaw_args}
BUILD_BYPRODUCTS ${forestclaw_byproducts}
INACTIVITY_TIMEOUT 15
UPDATE_DISCONNECTED true
CONFIGURE_HANDLED_BY_BUILD true)

target_link_libraries(FORESTCLAW::FORESTCLAW INTERFACE ${FORESTCLAW_LIBRARIES})
target_include_directories(FORESTCLAW::FORESTCLAW INTERFACE ${FORESTCLAW_INCLUDE_DIRS})

if(SC_FOUND)
  target_link_libraries(FORESTCLAW::FORESTCLAW INTERFACE SC::SC)
else()
  target_link_libraries(FORESTCLAW::FORESTCLAW INTERFACE ${SC_LIBRARIES})
endif()

if(P4EST_FOUND)
  target_link_libraries(FORESTCLAW::FORESTCLAW INTERFACE P4EST::P4EST)
else()
  target_link_libraries(FORESTCLAW::FORESTCLAW INTERFACE ${P4EST_LIBRARIES})
endif()

if(ZLIB_FOUND)
  target_link_libraries(FORESTCLAW::FORESTCLAW INTERFACE ZLIB::ZLIB)
endif()

target_link_libraries(FORESTCLAW::FORESTCLAW INTERFACE MPI::MPI_C)

target_link_libraries(FORESTCLAW::CLAWPATCH INTERFACE ${CLAWPATCH_LIBRARIES})
target_include_directories(FORESTCLAW::CLAWPATCH INTERFACE ${FORESTCLAW_INCLUDE_DIRS})

add_dependencies(FORESTCLAW::FORESTCLAW FORESTCLAW)
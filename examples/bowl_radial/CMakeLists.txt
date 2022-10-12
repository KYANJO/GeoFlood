set(src ${PROJECT_SOURCE_DIR}/src)

add_executable(bowl_radial 
	bowl_user.cpp
	bowl.cpp
)

# ---------------------------------------------------------
# This is needed to include *.h files needed from bowl_radial
# ---------------------------------------------------------
target_include_directories(bowl_radial PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${src} 
)

target_link_libraries(bowl_radial PUBLIC
  FORESTCLAW::GEOFLOOD
  FORESTCLAW::CLAWPACK4.6 
)

target_include_directories(bowl_radial PUBLIC ${FORESSTCLAW_INCLUDE_DIR})
target_link_libraries(bowl_radial PUBLIC ${FORESSTCLAW_LIBRARIES} ${FORESTCLAW_DEPENDENCIES})

add_test(NAME geoflood_bowl_radial COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(geoflood_bowl_radial PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")

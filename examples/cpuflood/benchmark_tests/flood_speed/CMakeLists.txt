
add_executable(flood_speed 
	flood_speed_user.cpp
	flood_speed_user.h
	flood_speed.cpp
)

target_include_directories(flood_speed PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
)

target_link_libraries(flood_speed PRIVATE 
  GEOFLOOD::GEOCLAW
)

# -- set environment string used for regression tests
set(FCLAW_TEST_ENVIRONMENT "FCLAW_EXAMPLES_BUILD_DIR=${PROJECT_BINARY_DIR};FCLAW_EXAMPLES_SRC_DIR=${PROJECT_SOURCE_DIR}")

add_test(NAME geoflood_flood_speed_cpu COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(geoflood_flood_speed_cpu PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")

add_library(malpasset_cpu_f OBJECT 
  setprob.f 
  rpn2cons_update.f
  user_4.6/qinit.f 
  rp/clawpack46_rpn2.f 
  rp/clawpack46_rpt2.f 
)

add_library(malpasset_cuda_f OBJECT 
  user_cuda/qinit.f 
)

add_executable(malpasset_cuda 
  malpasset.cpp 
  malpasset_options.c 
  malpasset_user.cpp 
	user_cuda/rpn2shallow.cu
  $<TARGET_OBJECTS:malpasset_cuda_f>
  $<TARGET_OBJECTS:malpasset_cpu_f>
)

target_include_directories(malpasset_cuda PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
)

target_link_libraries(malpasset_cuda PRIVATE 
 GEOFLOOD::CUDACLAW
 GEOFLOOD::GEOCLAW
)

add_test(NAME geoflood_malpasset_cpu COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(geoflood_malpasset_cpu PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")

add_test(NAME geoflood_malpasset_cuda COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(geoflood_malpasset_cuda PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")
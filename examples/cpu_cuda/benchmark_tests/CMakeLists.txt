# --------------------------------------------
# List of available bench marks
# --------------------------------------------

if(TARGET GEOFLOOD::GEOCLAW)

    include(benchmarks.cmake)

endif()

if(TARGET GEOFLOOD::CUDA_FLOOD)

    include(benchmarks.cmake)

endif()


if(TARGET GEOFLOOD::CPU_CUDA)

    include(benchmarks.cmake)

endif()
# --------------------------------------------
# List of available examples
# --------------------------------------------

#----------------------------------
# GeoFlood library and examples
# ----------------------------------

if(TARGET GEOFLOOD::GEOCLAW)

    include(cpuflood/cpuflood.cmake)

endif()

if(TARGET GEOFLOOD::CUDACLAW)
    enable_language(CUDA)
    include(cudaflood/cudaflood.cmake)

endif()
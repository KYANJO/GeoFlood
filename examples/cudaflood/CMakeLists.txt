# --------------------------------------------
# List of available cudaflood examples
# --------------------------------------------

if(TARGET GEOFLOOD::CUDACLAW)
    enable_language(CUDA)
    include(cudaflood.cmake)

endif()


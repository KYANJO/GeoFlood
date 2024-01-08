# --------------------------------------------
# List of available examples
# --------------------------------------------

#----------------------------------
# GeoFlood library and examples
# ----------------------------------

if(TARGET GEOFLOOD::GEOCLAW)

    include(cpuflood/cpuflood.cmake)

endif()


if(TARGET GEOFLOOD::CUDA_FLOOD)
    enable_language(CUDA)
    # include(cudaflood/cudaflood.cmake)
    include(cpu_cuda/cpu_cuda.cmake)
endif()

# if(GEOFLOOD_CUDAFLOOD_DEFINED)
#     enable_language(CUDA)
    # include(cudaflood/cudaflood.cmake)
    # include(cpu_cuda/cpu_cuda.cmake)
# endif()


# if(TARGET GEOFLOOD::CUDACLAW)
#     enable_language(CUDA)
#     include(cudaflood/cudaflood.cmake)
# endif()

if(TARGET GEOFLOOD::CPU_CUDA)
    enable_language(CUDA)
    # include(cudaflood/cudaflood.cmake)
    include(cpu_cuda/cpu_cuda.cmake)
endif()

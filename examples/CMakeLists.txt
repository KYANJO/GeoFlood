# --------------------------------------------
# List of available examples
# --------------------------------------------

#----------------------------------
# GeoFlood library and examples
# ----------------------------------
# cpu version
if(TARGET GEOFLOOD::GEOCLAW)

    include(cpuflood/cpuflood.cmake)

endif()

# hybrid cpu/cuda version
if(TARGET GEOFLOOD::CPU_CUDA)
    enable_language(CUDA)
    include(cpu_cuda/cpu_cuda.cmake)
endif()


#  -----------------------------
#  Switches to enable/disable cpu and/or cuda versions
#  -----------------------------

if(cpuflood)
    # option(cpuflood "Build cpu version" ON)
    add_subdirectory(fc2d_geoclaw)
endif(cpuflood)

if(cudaflood)
    # option(cudaflood "Build cuda version" ON)
    add_subdirectory(fc2d_cudaclaw)
endif(cudaflood)

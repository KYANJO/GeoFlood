#  -----------------------------
#  Switches to enable/disable cpu and/or cuda versions
#  -----------------------------

if(cpuflood)
    # option(cpuflood "Build cpu version" ON)
    add_subdirectory(fc2d_geoclaw)
endif(cpuflood)

if(cudaflood)
    # option(cudaflood "Build cuda version" ON)
    add_subdirectory(src/fc2d_cudaflood)
endif(cudaflood)

if(cpu_cuda)
    add_subdirectory(src/cpu_cudaflood)
    # add_subdirectory(src/fc2d_cudaflood)
endif(cpu_cuda)
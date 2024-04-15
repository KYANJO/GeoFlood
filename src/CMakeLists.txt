#  -----------------------------
#  Switches to enable/disable cpu and/or cuda versions
#  -----------------------------

if(cpuflood)
    # option(cpuflood "Build cpu version" ON)
    add_subdirectory(fc2d_geoclaw)
endif(cpuflood)

if(cpu_cuda)
    # option(cudaflood "Build cpu/cuda version" ON)
    add_subdirectory(src/fc2d_cpu_cuda)
endif(cpu_cuda)

# ------ examples/geoFlood.apps ----------------

## ----------------------------------
## GeoFlood examples
## ----------------------------------
if(cpu_cuda)
    add_subdirectory(cpu_cuda/benchmark_tests)
    add_subdirectory(cpu_cuda/malpasset) 
    add_subdirectory(cpu_cuda/teton)
    add_subdirectory(cpu_cuda/missoula)
    # add_subdirectory(cpu_cuda/bump)
endif(cpu_cuda)
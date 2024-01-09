# ------ examples/geoFlood.apps ----------------

## ----------------------------------
## GeoFlood examples
## ----------------------------------
if(cpu_cuda)
    add_subdirectory(cpu_cuda/benchmark_tests)
    # add_subdirectory(cpu_cuda/malpasset) #<--will switch on when implementation
endif(cpu_cuda)
# ------ examples/geoFlood.apps ----------------

## ----------------------------------
## GeoFlood examples
## ----------------------------------
# if(cudaflood)
    # add_subdirectory(cudaflood/teton)
    add_subdirectory(cpu_cuda/bump)
    add_subdirectory(cpu_cuda/benchmark_tests)
    # add_subdirectory(cudaflood/malpasset) #<--will switch on when implementation
# endif(cudaflood)
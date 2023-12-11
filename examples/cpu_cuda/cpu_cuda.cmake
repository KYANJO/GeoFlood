# ------ examples/geoFlood.apps ----------------

## ----------------------------------
## GeoFlood examples
## ----------------------------------
if(cudaflood)
    # add_subdirectory(cudaflood/teton)
    # add_subdirectory(cpu_cuda/bump) #!<-- different initial conditions
    add_subdirectory(cpu_cuda/benchmark_tests)
    # add_subdirectory(cudaflood/malpasset) #<--will switch on when implementation
endif(cudaflood)

if(cpu_cuda)
    # add_subdirectory(cudaflood/teton)
    # add_subdirectory(cpu_cuda/bump) #!<-- different initial conditions
    add_subdirectory(cpu_cuda/benchmark_tests)
    # add_subdirectory(cudaflood/malpasset) #<--will switch on when implementation
endif(cpu_cuda)
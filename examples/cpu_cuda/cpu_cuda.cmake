# ------ examples/geoFlood.apps ----------------

## ----------------------------------
## GeoFlood examples
## ----------------------------------
if(cpu_cuda)
    add_subdirectory(cpu_cuda/benchmark_tests)
    add_subdirectory(cpu_cuda/malpasset) 
    add_subdirectory(cpu_cuda/teton)
    add_subdirectory(cpu_cuda/hurricane_harvey)
endif(cpu_cuda)
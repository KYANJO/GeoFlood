# ------ examples/geoFlood.apps ----------------

## ----------------------------------
## GeoFlood examples
## ----------------------------------
if(cudaflood)
    add_subdirectory(teton)
    add_subdirectory(bump)
    # add_subdirectory(triton)
endif(cudaflood)
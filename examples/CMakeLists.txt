# --------------------------------------------
# List of available examples
# --------------------------------------------

#----------------------------------
# GeoClaw library and examples
# ----------------------------------

if(TARGET GEOFLOOD::GEOCLAW)

    include(geoflood.cmake)

endif()


# -- set environment string used for regression tests
# set(FCLAW_TEST_ENVIRONMENT "FCLAW_EXAMPLES_BUILD_DIR=${PROJECT_BINARY_DIR};FCLAW_EXAMPLES_SRC_DIR=${PROJECT_SOURCE_DIR}")



#!/bin/sh
# absolute path to example we are testing
examples=$FCLAW_EXAMPLES_BUILD_DIR/examples/bowl_radial/bowl_radial

# change to source dir for working directory
cd $FCLAW_EXAMPLES_SRC_DIR/examples/bowl_radial/regression

# run programs, exit script with nonzero on failure (or else script will exit with value of last program run)
$FCLAW_MPIRUN $FCLAW_MPI_TEST_FLAGS $examples -F regression.ini || exit 1


#define _SRC_FCLAW_CONFIG_H 1

/* C compiler */
#define FCLAW_CC "/opt/local/bin/mpicc"

/* C compiler flags */
#define FCLAW_CFLAGS "-O2 -g -Wall -pedantic -std=c99 "

/* C preprocessor */
#define FCLAW_CPP "/opt/local/bin/mpicc -E"

/* C preprocessor flags */
#define FCLAW_CPPFLAGS ""

/* Nvidia C/C++ flags */
#define FCLAW_CUDA_CFLAGS 

/* Nvidia LD flags */
#define FCLAW_CUDA_LDFLAGS 

/* C++ compiler */
#define FCLAW_CXX "/opt/local/bin/mpicxx"

/* C++ compiler flags */
#define FCLAW_CXXFLAGS ""

/* Undefine if: use aligned malloc (optionally use --enable-memalign=<bytes>) */
#define FCLAW_ENABLE_MEMALIGN

/* Define to 1 if we are using MPI */
#define FCLAW_ENABLE_MPI

/* Define to 1 if we are using MPI I/O */
#define FCLAW_ENABLE_MPIIO


/* F77 compiler */

#define FCLAW_F77 ""


/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */

#define FCLAW_F77_FUNC(name,NAME) name ## _


/* As F77_FUNC, but for C identifiers containing underscores. */

#define FCLAW_F77_FUNC_(name,NAME) name ## _


/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */

#define FCLAW_FC_FUNC(name,NAME) name ## _


/* As FC_FUNC, but for C identifiers containing underscores. */
#define FCLAW_FC_FUNC_(name,NAME) name ## _


/* F77 compiler flags */
#define FCLAW_FFLAGS "-O2 -g -Wall -Wno-unused-dummy-argument -cpp "


/* Fortran libraries */
#define FCLAW_FLIBS ""


/* Define to 1 if you have the `feenableexcept' function. */
/* #undef FCLAW_HAVE_FEENABLEEXCEPT */

/* Define to 1 if you have the <fenv.h> header file. */
#define FCLAW_HAVE_FENV_H

/* Define to 1 if you have the <signal.h> header file. */
#define FCLAW_HAVE_SIGNAL_H

/* Define to 1 if you have the <unistd.h> header file. */
#define FCLAW_HAVE_UNISTD_H


/* Linker flags */

#define FCLAW_LDFLAGS ""


/* Libraries */

#define FCLAW_LIBS ";/opt/local/lib/libz.dylib;m"



/* Name of package */

#define FCLAW_PACKAGE "forestclaw"


/* Define to the address where bug reports for this package should be sent. */

#define FCLAW_PACKAGE_BUGREPORT "burstedde@ins.uni-bonn.de"


/* Define to the full name of this package. */

#define FCLAW_PACKAGE_NAME "ForestClaw"


/* Define to the full name and version of this package. */

#define FCLAW_PACKAGE_STRING "ForestClaw 0.1.4880-4dae"


/* Define to the one symbol short name of this package. */

#define FCLAW_PACKAGE_TARNAME "forestclaw"


/* Define to the home page for this package. */

#define FCLAW_PACKAGE_URL ""


/* Define to the version of this package. */
#define FCLAW_PACKAGE_VERSION "0.1.4880"

/* Version number of package */
#define FCLAW_VERSION "0.1.4880"

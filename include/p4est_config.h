#define _SRC_P_EST_CONFIG_H 1

/* C compiler */
#define P4EST_CC "/opt/local/bin/mpicc"

/* C compiler flags */
#define P4EST_CFLAGS " "

/* C preprocessor */
#define P4EST_CPP "/opt/local/bin/mpicc -E"

/* C preprocessor flags */
#define P4EST_CPPFLAGS ""

/* Undefine if: disable the 2D library */
#define P4EST_ENABLE_BUILD_2D

/* Undefine if: disable the 3D library */
#define P4EST_ENABLE_BUILD_3D

/* Undefine if: disable hybrid 2D+1D p6est library */
#define P4EST_ENABLE_BUILD_P6EST

/* Define to 1 if we are using debug build type (assertions and extra checks) */
/* #undef P4EST_ENABLE_DEBUG */

/* Undefine if: use aligned malloc (optionally use --enable-memalign=<bytes>) */
#define P4EST_ENABLE_MEMALIGN

/* Define to 1 if we are using MPI */
#define P4EST_ENABLE_MPI

/* Define to 1 if we can use MPI_COMM_TYPE_SHARED */
#define P4EST_ENABLE_MPICOMMSHARED

/* Define to 1 if we are using MPI I/O */
#define P4EST_ENABLE_MPIIO

/* Define to 1 if we are using MPI_Init_thread */
#define P4EST_ENABLE_MPITHREAD

/* Define to 1 if we can use MPI_Win_allocate_shared */
#define P4EST_ENABLE_MPIWINSHARED


/* Undefine if: write vtk ascii file data */
#define P4EST_ENABLE_VTK_BINARY 1

/* Undefine if: disable zlib compression for vtk binary data */
#define P4EST_ENABLE_VTK_COMPRESSION 1

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define P4EST_F77_FUNC(name,NAME) name ## _

/* As F77_FUNC, but for C identifiers containing underscores. */
#define P4EST_F77_FUNC_(name,NAME) name ## _

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define P4EST_FC_FUNC(name,NAME) name ## _

/* As FC_FUNC, but for C identifiers containing underscores. */
#define P4EST_FC_FUNC_(name,NAME) name ## _


/* Define to 1 if you have the <arpa/inet.h> header file. */
#define P4EST_HAVE_ARPA_INET_H

/* Define to 1 if you have the <dlfcn.h> header file. */
#define P4EST_HAVE_DLFCN_H

/* Define to 1 if you have the `fsync' function. */
#define P4EST_HAVE_FSYNC

/* Define to 1 if you have the <inttypes.h> header file. */
#define P4EST_HAVE_INTTYPES_H

/* Have we found function pthread_create. */
#define HAVE_LPTHREAD

/* Have we found function lua_createtable. */
/* #undef HAVE_LUA */

/* Define to 1 if you have the <memory.h> header file. */
#define P4EST_HAVE_MEMORY_H

/* Define to 1 if you have the <netinet/in.h> header file. */
#define P4EST_HAVE_NETINET_IN_H

/* Define if you have the <Winsock2.h> header file. */
/* #undef P4EST_HAVE_WINSOCK2_H */

/* Define to 1 if you have the `posix_memalign' function. */
#define P4EST_HAVE_POSIX_MEMALIGN

/* Define to 1 if you have the <stdint.h> header file. */
#define P4EST_HAVE_STDINT_H

/* Define to 1 if you have the <stdlib.h> header file. */
#define P4EST_HAVE_STDLIB_H

/* Define to 1 if you have the <strings.h> header file. */
#define P4EST_HAVE_STRINGS_H

/* Define to 1 if you have the <string.h> header file. */
#define P4EST_HAVE_STRING_H

/* Define to 1 if you have the <sys/stat.h> header file. */
#define P4EST_HAVE_SYS_STAT_H

/* Define to 1 if you have the <sys/types.h> header file. */
#define P4EST_HAVE_SYS_TYPES_H

/* Define to 1 if you have the <unistd.h> header file. */
#define P4EST_HAVE_UNISTD_H

/* Have we found function adler32_combine. */
#define P4EST_HAVE_ZLIB

/* Linker flags */
#define P4EST_LDFLAGS "-Wl,-syslibroot,/Library/Developer/CommandLineTools/SDKs/MacOSX12.sdk -Wl,-flat_namespace -Wl,-commons,use_dylibs -L/opt/local/lib"

/* Libraries */
#define P4EST_LIBS ";/opt/local/lib/libz.dylib;m"


/* Name of package */
#define P4EST_PACKAGE "p4est"

/* Define to the address where bug reports for this package should be sent. */
#define P4EST_PACKAGE_BUGREPORT "info@p4est.org"

/* Define to the full name of this package. */
#define P4EST_PACKAGE_NAME "p4est"

/* Define to the full name and version of this package. */
#define P4EST_PACKAGE_STRING "p4est 2.0.94-00da"

/* Define to the one symbol short name of this package. */
#define P4EST_PACKAGE_TARNAME "p4est"

/* Define to the home page for this package. */
#define P4EST_PACKAGE_URL ""

/* Define to the version of this package. */
#define P4EST_PACKAGE_VERSION "2.0.94-00da"

/* Version number of package */
#define P4EST_VERSION "2.0.94-00da"

/* Package major version */
#define P4EST_VERSION_MAJOR 2

/* Package minor version */
#define P4EST_VERSION_MINOR 0

/* Package point version */
#define P4EST_VERSION_POINT 94-00da

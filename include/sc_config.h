#ifndef _SRC_SC_CONFIG_H
#define _SRC_SC_CONFIG_H 1

/* C compiler */
#ifndef SC_CC
#define SC_CC "/opt/local/bin/mpicc"
#endif

/* C compiler flags */
#ifndef SC_CFLAGS
#define SC_CFLAGS " "
#endif

/* C preprocessor */
#ifndef SC_CPP
#define SC_CPP "/opt/local/bin/mpicc -E"
#endif

/* C preprocessor flags */
#ifndef SC_CPPFLAGS
#define SC_CPPFLAGS ""
#endif

#define SC_ENABLE_PTHREAD

/* Define to 1 if we are using debug build type (assertions and extra checks) */
/* #undef SC_ENABLE_DEBUG */

/* Undefine if: use aligned malloc (optionally use --enable-memalign=<bytes>) */
#define SC_ENABLE_MEMALIGN

/* Define to 1 if we are using MPI */
#define SC_ENABLE_MPI

/* Define to 1 if we can use MPI_COMM_TYPE_SHARED */
#define SC_ENABLE_MPICOMMSHARED

/* Define to 1 if we are using MPI I/O */
#define SC_ENABLE_MPIIO

/* Define to 1 if we are using MPI_Init_thread */
#define SC_ENABLE_MPITHREAD

/* Define to 1 if we can use MPI_Win_allocate_shared */
#define SC_ENABLE_MPIWINSHARED

/* Undefine if: disable non-thread-safe internal debug counters */
#define SC_ENABLE_USE_COUNTERS 1


/* Undefine if: replace array/dmatrix resize with malloc/copy/free */
#define SC_ENABLE_USE_REALLOC

/* Define to 1 if you have the `aligned_alloc' function. */
#define SC_HAVE_ALIGNED_ALLOC

/* Define to 1 if you have the `_aligned_,alloc' function. */
/* #undef SC_HAVE_ALIGNED_MALLOC */

/* Define to 1 if you have the `backtrace' function. */
#define SC_HAVE_BACKTRACE

/* Define to 1 if you have the `backtrace_symbols' function. */
#define SC_HAVE_BACKTRACE_SYMBOLS

/* Define to 1 if you have the <dlfcn.h> header file. */
/* #undef SC_HAVE_DLFCN_H */

/* Define to 1 if you have the <execinfo.h> header file. */
#define SC_HAVE_EXECINFO_H

/* Define to 1 if you have the `fsync' function. */
#define SC_HAVE_FSYNC

/* Define to 1 if you have the <inttypes.h> header file. */
// cmakedefine SC_HAVE_INTTYPES_H


/* Define to 1 if you have the <memory.h> header file. */
/* #undef SC_HAVE_MEMORY_H */



/* Define to 1 if you have the `posix_memalign' function. */
#define SC_HAVE_POSIX_MEMALIGN

/* Define to 1 if you have the `qsort_r' function. */
#define SC_HAVE_QSORT_R

/* Define to 1 if you have the <signal.h> header file. */
#define SC_HAVE_SIGNAL_H

#define SC_HAVE_STDINT_H

#define SC_HAVE_STDLIB_H

#define SC_HAVE_STRING_H

/* Define to 1 if you have the `strtoll' function. */
#define SC_HAVE_STRTOLL

/* Define to 1 if you have the <sys/time.h> header file. */
#define SC_HAVE_SYS_TIME_H

/* Define to 1 if you have the <sys/types.h> header file. */
#define SC_HAVE_SYS_TYPES_H

/* Define to 1 if you have the <time.h> header file. */
#define SC_HAVE_TIME_H

/* #undef SC_HAVE_WINSOCK2_H */


/* desired alignment of allocations in bytes */
#define SC_MEMALIGN_BYTES (SC_SIZEOF_VOID_P)

/* Define to 1 if you have the <unistd.h> header file. */
#define SC_HAVE_UNISTD_H

/* Have we found function adler32_combine. */
#define SC_HAVE_ZLIB

/* Linker flags */
#ifndef SC_LDFLAGS
#define SC_LDFLAGS "-Wl,-syslibroot,/Library/Developer/CommandLineTools/SDKs/MacOSX12.sdk -Wl,-flat_namespace -Wl,-commons,use_dylibs -L/opt/local/lib"
#endif

/* Libraries */
#ifndef SC_LIBS
#define SC_LIBS "/opt/local/lib/libz.dylib;m"
#endif


/* Name of package */
#ifndef SC_PACKAGE
#define SC_PACKAGE "libsc"
#endif

/* Define to the address where bug reports for this package should be sent. */
#ifndef SC_PACKAGE_BUGREPORT
#define SC_PACKAGE_BUGREPORT "info@SC.org"
#endif

/* Define to the full name of this package. */
#ifndef SC_PACKAGE_NAME
#define SC_PACKAGE_NAME "libsc"
#endif

/* Define to the full name and version of this package. */
#ifndef SC_PACKAGE_STRING
#define SC_PACKAGE_STRING "libsc 2.8.1.57-a813"
#endif

/* Define to the one symbol short name of this package. */
#ifndef SC_PACKAGE_TARNAME
#define SC_PACKAGE_TARNAME "libsc"
#endif

/* Define to the home page for this package. */
#ifndef SC_PACKAGE_URL
#define SC_PACKAGE_URL ""
#endif

/* Define to the version of this package. */
#ifndef SC_PACKAGE_VERSION
#define SC_PACKAGE_VERSION "2.8.1.57-a813"
#endif

/* Byte sizes of standard types */
#ifndef SC_SIZEOF_INT
#define SC_SIZEOF_INT 4
#endif
#ifndef SC_SIZEOF_UNSIGNED_INT
#define SC_SIZEOF_UNSIGNED_INT 4
#endif
#ifndef SC_SIZEOF_LONG
#define SC_SIZEOF_LONG 8
#endif
#ifndef SC_SIZEOF_LONG_LONG
#define SC_SIZEOF_LONG_LONG 8
#endif
#ifndef SC_SIZEOF_UNSIGNED_LONG
#define SC_SIZEOF_UNSIGNED_LONG 8
#endif
#ifndef SC_SIZEOF_UNSIGNED_LONG_LONG
#define SC_SIZEOF_UNSIGNED_LONG_LONG 8
#endif

/* The size of `void *', as computed by sizeof. */
#ifndef SC_SIZEOF_VOID_P
#define SC_SIZEOF_VOID_P 8
#endif


/* Version number of package */
#ifndef SC_VERSION
#define SC_VERSION "2.8.1.57-a813"
#endif

/* Package major version */
#ifndef SC_VERSION_MAJOR
#define SC_VERSION_MAJOR 2
#endif

/* Package minor version */
#ifndef SC_VERSION_MINOR
#define SC_VERSION_MINOR 8
#endif

/* Package point version */
#ifndef SC_VERSION_POINT
#define SC_VERSION_POINT 1.57-a813
#endif

/* once: _SRC_SC_CONFIG_H */
#endif

## Makefile.am in examples/triton

bin_PROGRAMS += examples/triton/triton

examples_triton_triton_SOURCES = \
	examples/triton/triton_user.cpp \
	examples/triton/triton_user.h \
	examples/triton/triton.cpp


examples_triton_triton_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_GEOCLAW_CPPFLAGS)

examples_triton_triton_LDADD = \
        $(FCLAW_GEOCLAW_LDADD) $(LDADD)

## @top_builddir@/src/solvers/fc2d_geoclaw/libgeoclaw.la

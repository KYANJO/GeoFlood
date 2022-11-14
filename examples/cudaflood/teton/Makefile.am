## Makefile.am in examples/teton

bin_PROGRAMS += examples/teton/teton

examples_teton_teton_SOURCES = \
	examples/teton/teton_user.cpp \
	examples/teton/teton_user.h \
	examples/teton/reservoir_module.f90 \
	examples/teton/qinit.f90 \
	examples/teton/teton.cpp


examples_teton_teton_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_GEOCLAW_CPPFLAGS)

examples_teton_teton_LDADD = \
        $(FCLAW_GEOCLAW_LDADD) $(LDADD)

## @top_builddir@/src/solvers/fc2d_geoclaw/libgeoclaw.la

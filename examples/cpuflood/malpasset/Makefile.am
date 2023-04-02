## Makefile.am in examples/malpasset

bin_PROGRAMS += examples/malpasset/malpasset

examples_malpasset_malpasset_SOURCES = \
	examples/malpasset/malpasset_user.cpp \
	examples/malpasset/malpasset_user.h \
	examples/malpasset/malpasset.cpp


examples_malpasset_malpasset_CPPFLAGS = \
        $(AM_CPPFLAGS) $(FCLAW_GEOCLAW_CPPFLAGS)

examples_malpasset_malpasset_LDADD = \
        $(FCLAW_GEOCLAW_LDADD) $(LDADD)

## @top_builddir@/src/solvers/fc2d_geoclaw/libgeoclaw.la

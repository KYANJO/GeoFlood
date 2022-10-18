## Makefile.am in examples/malpasset

bin_PROGRAMS += examples/malpasset/malpasset

TESTS += examples/malpasset/regressions.sh

examples_malpasset_malpasset_SOURCES = \
	examples/malpasset/bowl_user.h \
	examples/malpasset/bowl.cpp

EXTRA_DIST += \
	examples/malpasset/fclaw_options.ini \
	examples/malpasset/setprun.py \
	examples/malpasset/setplot.py \
	examples/malpasset/make_topo.py \
	examples/malpasset/afterframe.m \
	examples/malpasset/setplot2.m \

examples_malpasset_malpasset_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_GEOCLAW_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

examples_malpasset_malpasset_LDADD = \
        $(LDADD) \
        $(FCLAW_GEOCLAW_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)
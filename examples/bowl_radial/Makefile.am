## Makefile.am in examples/bowl_radial

bin_PROGRAMS += examples/bowl_radial/bowl_radial

TESTS += examples/bowl_radial/regressions.sh

examples_bowl_radial_bowl_radial_SOURCES = \
	examples/bowl_radial/bowl_user.cpp \
	examples/bowl_radial/bowl_user.h \
	examples/bowl_radial/bowl.cpp

EXTRA_DIST += \
	examples/bowl_radial/fclaw_options.ini \
	examples/bowl_radial/make_data.py \
	examples/bowl_radial/make_plots.py \
	examples/bowl_radial/make_topo.py \
	examples/bowl_radial/afterframe.m \
	examples/bowl_radial/setplot2.m \
	examples/bowl_radial/bathy.m \
	examples/bowl_radial/bathy1d.m \
	examples/bowl_radial/map1d.m \
	examples/bowl_radial/compare_gauges.m

examples_bowl_radial_bowl_radial_CPPFLAGS = \
        $(AM_CPPFLAGS) \
        $(FCLAW_GEOCLAW_CPPFLAGS) \
        $(FCLAW_CLAWPATCH_CPPFLAGS)

examples_bowl_radial_bowl_radial_LDADD = \
        $(LDADD) \
        $(FCLAW_GEOCLAW_LDADD) \
        $(FCLAW_CLAWPATCH_LDADD) \
        $(FCLAW_LDADD)


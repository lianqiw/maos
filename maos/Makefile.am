SUBDIRS = . 
#AM_CFLAGS mean CFLAGS from automake, it is used in addition to user CFLAGS if per target CFLAGS is not specified.
AM_CFLAGS=$(OPTS)
AM_LDFLAGS=$(LDEXE)

#the following libs can not be in AM_LDFLAGS otherwise in mac the binary will not be updated when the libs change (in MAC)#Notice the order. previous dependent on later

LDADD=../lib/libaos.la ../lib/sys/libsys.la  ../external/cholmod/libchol.la ../external/arpack/libarpack.la $(LDOPTS)

bin_PROGRAMS=../bin/maos #run/test_recon run/test_powfs_ints
___bin_maos_SOURCES=maos.h maos.c parms.h setup_parms.c types.h utils.h utils.c	\
	 setup_powfs.c setup_powfs.h setup_aper.c setup_aper.h setup_recon.c	\
	 setup_recon.h sim.c sim.h recon.c recon.h fdpcg.c fdpcg.h wfsgrad.c	\
	 perfevl.c wfsints.c genseotf.c genseotf.h mtch.c mtch.h cachedm.c	\
	 filter.c ahst.c ahst.h sim_utils.c sim_utils.h setup_surf.c		\
	 setup_surf.h recon_utils.c recon_utils.h cn2est.c cn2est.h


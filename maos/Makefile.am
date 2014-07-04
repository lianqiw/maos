SUBDIRS = . 
#AM_CFLAGS mean CFLAGS from automake, it is used in addition to user CFLAGS if per target CFLAGS is not specified.
AM_CFLAGS=$(OPTS)

#the libs can not be in AM_LDFLAGS otherwise in mac the binary will not be updated when the libs change (in MAC)#Notice the order. previous dependent on later
nobase_include_HEADERS=common.h maos.h parms.h types.h utils.h setup_powfs.h \
	setup_aper.h setup_recon.h recon.h  sim.h fdpcg.h genseotf.h	\
	mtch.h ahst.h sim_utils.h setup_surf.h save.h mvm_client.h moao.h \
	recon_utils.h
lib_LTLIBRARIES=libmaos.la
libmaos_la_SOURCES=maos.c setup_parms.c utils.c setup_powfs.c setup_aper.c	\
setup_recon.c setup_recon_dm.c setup_recon_fit.c setup_recon_lsr.c		\
setup_recon_mvm.c sim.c recon.c fdpcg.c wfsgrad.c perfevl.c wfsints.c		\
genseotf.c mtch.c cachedm.c filter.c ahst.c sim_utils.c setup_surf.c		\
recon_utils.c moao.c save.c mvm_client.c
#The following must be specified like this. do not use LDAOS set in configure.ac that causes linking to be not performed when lib is updated.
libmaos_la_LIBADD=
if USE_CUDA
libmaos_la_LIBADD+=../cuda/libaogpu.la
endif
libmaos_la_LIBADD+=../lib/libaos.la 

bin_PROGRAMS=../bin/maos 
___bin_maos_SOURCES=main.c
___bin_maos_LDADD=libmaos.la $(LDOPTS)
if USE_CUDA
___bin_maos_LDADD+=$(LDCUDA)
endif

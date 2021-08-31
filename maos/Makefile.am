SUBDIRS=../lib
if USE_CUDA
SUBDIRS+=../cuda
endif

AM_CFLAGS=$(OPTS) -I $(top_srcdir) #Compiler flags for each target that doesn't specify per target CFLAGS. Supplements CFLAGS
AM_LDFLAGS=$(LDOPTS)
#the libs can not be in AM_LDFLAGS otherwise in mac the binary will not be updated when the libs change (in MAC)#Notice the order. previous dependent on later
nobase_include_HEADERS=common.h maos.h parms.h types.h utils.h setup_powfs.h \
	setup_aper.h setup_recon.h recon.h  sim.h fdpcg.h powfs_utils.h	\
	ahst.h sim_utils.h setup_surf.h save.h mvm_client.h moao.h \
	recon_utils.h pywfs.h
SRCS=maos.c setup_parms.c utils.c setup_powfs.c setup_aper.c	\
setup_recon.c setup_recon_prep.c setup_recon_fit.c setup_recon_lsr.c		\
setup_recon_mvm.c sim.c recon.c fdpcg.c wfsgrad.c perfevl.c wfsints.c		\
powfs_utils.c cachedm.c filter.c ahst.c sim_utils.c setup_surf.c		\
recon_utils.c moao.c save.c mvm_client.c pywfs.c

lib_LTLIBRARIES=libmaos.la
libmaos_la_SOURCES=$(SRCS)
libmaos_la_LIBADD=
if USE_CUDA
libmaos_la_LIBADD+=../cuda/libaogpu.la
#if USE_STATIC
#libmaos_la_LIBADD+=../cuda/dlink.o
#endif
endif
libmaos_la_LIBADD+=../lib/libaos.la 
bin_PROGRAMS=../bin/maos 
___bin_maos_SOURCES=main.c
___bin_maos_LDADD=libmaos.la  
if USE_CUDA
___bin_maos_LDADD+=$(libmaos_la_LIBADD) $(LIBCUDA)
#if USE_STATIC
#___bin_maos_LDADD+=-lculibos -ldl -lrt
#endif
endif
___bin_maos_LDADD+=$(LDEXE_MAOS) $(LDEXE_FULL) 

AM_CFLAGS=-I $(top_srcdir) #Compiler flags for each target that doesn't specify per target CFLAGS. Supplements CFLAGS
AM_LDFLAGS=$(LDCUDA)
#the libs can not be in AM_LDFLAGS otherwise in mac the binary will not be updated when the libs change (in MAC)#Notice the order. previous dependent on later
nobase_include_HEADERS=common.h maos.h parms.h types.h utils.h powfs.h \
	aper.h recon.h  sim.h fdpcg.h powfs_utils.h	ahst.h sim_utils.h surf.h save.h mvm_client.h moao.h \
	recon_utils.h 
SRCS=setup.c parms.c utils.c powfs.c aper.c	sim.c fdpcg.c wfsgrad.c perfevl.c wfsints.c	\
    recon_setup.c recon_prep.c recon_fit.c recon_lsr.c recon_mvm.c recon.c recon_utils.c\
    powfs_utils.c cachedm.c filter.c ahst.c sim_utils.c surf.c		\
    moao.c save.c mvm_client.c

lib_LTLIBRARIES=libmaos.la
libmaos_la_SOURCES=$(SRCS)
libmaos_la_LIBADD=
if USE_CUDA
libmaos_la_LIBADD+=../cuda/recon/libcurecon.la
#if USE_STATIC
#libmaos_la_LIBADD+=../cuda/dlink.o
#endif
endif
libmaos_la_LIBADD+=../lib/libaos.la

bin_PROGRAMS=../bin/maos
___bin_maos_SOURCES=maos.c
___bin_maos_LDADD=libmaos.la
#if USE_CUDA
#___bin_maos_LDADD+=$(LIBCUDA)
#if USE_STATIC
#___bin_maos_LDADD+=-lculibos -ldl -lrt
#endif
#endif
___bin_maos_LDADD+=$(LDEXE_MAOS) $(LDEXE_FULL)


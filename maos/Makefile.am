SUBDIRS = . 
#AM_CFLAGS mean CFLAGS from automake, it is used in addition to user CFLAGS if per target CFLAGS is not specified.
AM_CFLAGS=$(OPTS)
AM_LDFLAGS=$(LDEXE) 

#the following libs can not be in AM_LDFLAGS otherwise in mac the binary will not be updated when the libs change (in MAC)#Notice the order. previous dependent on later
bin_PROGRAMS=../bin/maos 
___bin_maos_SOURCES=maos.h maos.c parms.h setup_parms.c types.h utils.h utils.c	\
	 setup_powfs.c setup_powfs.h setup_aper.c setup_aper.h setup_recon.c	\
	 setup_recon.h sim.c sim.h recon.c recon.h fdpcg.c fdpcg.h wfsgrad.c	\
	 perfevl.c wfsints.c genseotf.c genseotf.h mtch.c mtch.h cachedm.c	\
	 filter.c ahst.c ahst.h sim_utils.c sim_utils.h setup_surf.c		\
	 setup_surf.h recon_utils.c recon_utils.h moao.c moao.h cn2est.c	\
	 cn2est.h save.c save.h
#The following must be specified like this. do not use LDAOS set in configure.ac that causes linking to be not performed when lib is updated.
if USE_LONG
___bin_maos_LDADD=../lib/libaos_long.la $(LDOPTS)
else
___bin_maos_LDADD=../lib/libaos_int.la $(LDOPTS)
endif
if USE_CUDA
___bin_maos_LDADD+=../cuda/libgpu.a $(CUDA_OPT_L)
endif
#if USE_STATIC
#config.o:../config.tar.gz
#	objcopy --input binary --output $(objcopy_output) --binary-architecture $(objcopy_arch) ../config.tar.gz config.o
#___bin_maos_LDADD+=config.o
#endif
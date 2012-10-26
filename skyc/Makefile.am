SUBDIRS = . 
#AM_CFLAGS is used in addition to user CFLAGS if per target CFLAGS is not specified.
AM_CFLAGS=$(OPTS)
AM_LDFLAGS=$(LDEXE)

bin_PROGRAMS=../bin/skyc 
___bin_skyc_SOURCES=parms.h types.h setup_parms.c  skyc.c  skyc.h  utils.c  utils.h skysim.c skysim.h setup_powfs.c setup_powfs.h setup_aster.c setup_aster.h photon.c photon.h skysim_utils.c skysim_utils.h mtch.c mtch.h  nafocus.c nafocus.h genstars.c genstars.h genpistat.c genpistat.h setup_star.c setup_star.h
___bin_skyc_LDADD=
if USE_LONG
___bin_skyc_LDADD+=../lib/libaos_long.la
else
___bin_skyc_LDADD+=../lib/libaos_int.la 
endif
___bin_skyc_LDADD+=../sys/libaosys.la  
___bin_skyc_LDADD+=$(LDOPTS)

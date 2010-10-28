SUBDIRS = . 
#AM_CFLAGS is used in addition to user CFLAGS if per target CFLAGS is not specified.
AM_CFLAGS=$(OPTS)
AM_LDFLAGS=$(LDEXE)

#AM_LDFLAGS=$(LDOPTS) 
#Notice the order. previous dependent on later
LDADD=../lib/libaos.la ../lib/sys/libsys.la  ../external/cholmod/libchol.la ../external/arpack/libarpack.la $(LDOPTS)

bin_PROGRAMS=../bin/skyc ../bin/servo_optim
___bin_skyc_SOURCES=parms.h types.h setup_parms.c  skyc.c  skyc.h  utils.c  utils.h skysim.c skysim.h setup_powfs.c setup_powfs.h setup_aster.c setup_aster.h photon.c photon.h skysim_utils.c skysim_utils.h mtch.c mtch.h servo.c servo.h nafocus.c nafocus.h genstars.c genstars.h genpistat.c genpistat.h setup_star.c setup_star.h

___bin_servo_optim_SOURCES=servo.c
___bin_servo_optim_CFLAGS=$(AM_CFLAGS) -DTEST

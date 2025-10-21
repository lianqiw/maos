AM_CFLAGS=-I$(top_srcdir) $(CHOL_CFLAGS)
#if set LDCHOL here instead of in configure.ac compiles "No rule to make target ..."

nobase_include_HEADERS=accphi.h aos.h cn2est.h cure.h draw.h fractal.h genotf.h gensei.h hyst.h kalman.h laplacian.h locfft.h \
		       libmisc.h mkdtf.h mkg.h mkh.h mkw.h mkz.h mtch.h chol.h muv.h pcg.h proj.h psd.h readcfg.h servo.h slaving.h stfun.h \
			   turbulence.h zernike.h petal.h phase.h pywfs.h comb.h

EXTRA_DIST=fractal_do.c accphi_grid_do.c 

AOS_SRC=readcfg.c turbulence.c mkg.c mkh.c mkz.c mkw.c accphi.c proj.c		\
laplacian.c chol.c pcg.c muv.c draw.c genotf.c gensei.c fractal.c stfun.c servo.c hyst.c	\
slaving.c mkdtf.c psd.c kalman.c cn2est.c libmisc.c locfft.c zernike.c cure.c mtch.c petal.c phase.c  pywfs.c \
comb.c

lib_LTLIBRARIES=libaos.la
libaos_la_SOURCES=$(AOS_SRC)
libaos_la_LIBADD=../math/libaomath.la $(CHOL_LIBS)

bin_PROGRAMS=../bin/aolib.so
___bin_aolib_so_SOURCES=aolib.c
___bin_aolib_so_LDADD=libaos.la $(FFTW_LIBS) $(LAPACK) $(LDEXE)
___bin_aolib_so_LDFLAGS=-static -Xcompiler -shared $(LDSO)
#Linker flags depends on system type (native linker is used)
if LINUX
___bin_aolib_so_LDFLAGS+=-lrt -Wl,--whole-archive,.libs/libaos.a,../math/.libs/libaomath.a,../sys/.libs/libaosys.a,--no-whole-archive
else
if DARWIN
___bin_aolib_so_LDFLAGS+=-Wl,-force_load,.libs/libaos.a,../math/.libs/libaomath.a,../sys/.libs/libaosys.a
endif
endif

#BUILT_SOURCES=aolib.h
#aolib.h:$(nobase_include_HEADERS)
#	gcc -E $(SRC_DIR)/lib/aos.h > $(SRC_DIR)/lib/aolib.h


AM_CFLAGS=$(OPTS) -I$(top_srcdir) -I$(top_srcdir)/external/cholmod/Include
if !WIN
AM_CFLAGS+=-fPIC
endif
#if set LDCHOL here instead of in configure.ac compiles "No rule to make target ..."

nobase_include_HEADERS=accphi.h locbin.h mkg.h readcfg.h \
  mkh.h aos.h mkw.h turbulence.h loc.h mkz.h proj.h laplacian.h pcg.h muv.h draw.h \
 genotf.h fractal.h stfun.h servo.h slaving.h mkdtf.h psd.h kalman.h type.h

EXTRA_DIST=fractal_do.c prop_grid_stat.c prop_grid_pts.c prop_grid_map.c 

AOS_SRC=readcfg.c loc.c locbin.c turbulence.c mkg.c mkh.c mkz.c mkw.c accphi.c proj.c laplacian.c \
pcg.c muv.c draw.c genotf.c fractal.c stfun.c servo.c slaving.c mkdtf.c psd.c kalman.c


if USE_LONG
lib_LTLIBRARIES=libaos_long.la
libaos_long_la_SOURCES=$(AOS_SRC)
libaos_long_la_LIBADD=$(LDLIB)
libaos_long_la_LDFLAGS = -avoid-version 
else
lib_LTLIBRARIES=libaos_int.la
libaos_int_la_SOURCES=$(AOS_SRC)
libaos_int_la_LIBADD=$(LDLIB)
libaos_int_la_LDFLAGS = -avoid-version
endif

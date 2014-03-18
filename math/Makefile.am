AM_CFLAGS=$(OPTS) -I$(top_srcdir) -I$(top_srcdir)/external/cholmod/Include
if !WIN
AM_CFLAGS+=-fPIC
endif
#if set LDCHOL here instead of in configure.ac compiles "No rule to make target ..."

nobase_include_HEADERS= type.h  mathmisc.h random.h mathdef.h cell.h  suitesparse.h blas.h \
 matbin.h spbin.h cellarr.h chol.h  mat.h sp.h fft.h imat.h 

EXTRA_DIST=cell.c sp.c matbin.c matmmap.c spbin.c blas.c mat.c cmat.c \
 suitesparse.c  defs.h fft.c bessik.c amoeba.c

AOMATH_SRC= dmath.c cmath.c smath.c zmath.c mathmisc.c random.c cellarr.c imat.c chol.c


if USE_LONG
lib_LTLIBRARIES=libaomath_long.la
libaomath_long_la_SOURCES=$(AOMATH_SRC)
libaomath_long_la_LIBADD=$(LDLIB)
libaomath_long_la_LDFLAGS = -avoid-version 
else
lib_LTLIBRARIES=libaomath_int.la
libaomath_int_la_SOURCES=$(AOMATH_SRC)
libaomath_int_la_LIBADD=$(LDLIB)
libaomath_int_la_LDFLAGS = -avoid-version
endif

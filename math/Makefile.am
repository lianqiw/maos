AM_CFLAGS=$(OPTS) -I$(top_srcdir)
#if set LDCHOL here instead of in configure.ac compiles "No rule to make target ..."

nobase_include_HEADERS= type.h  mathmisc.h random.h mathdef.h cell.h  suitesparse.h blas.h \
 matbin.h spbin.h cellarr.h chol.h  mat.h sp.h fft.h matmath.h loc.h locbin.h

EXTRA_DIST=sp.c matbin.c matmmap.c spbin.c blas.c mat.c matmath.c cmat.c \
 suitesparse.c  defs.h fft.c bessik.c amoeba.c

AOMATH_SRC= dmath.c cmath.c smath.c zmath.c lmath.c mathmisc.c random.c cellarr.c chol.c loc.c locbin.c cell.c

lib_LTLIBRARIES=libaomath.la
libaomath_la_SOURCES=$(AOMATH_SRC)
libaomath_la_LIBADD=../sys/libaosys.la
libaomath_la_LDFLAGS = $(LDLIB)

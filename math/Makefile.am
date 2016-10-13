AM_CFLAGS=$(OPTS) -I$(top_srcdir)
AM_LDFLAGS=$(LDOPTS)

nobase_include_HEADERS= type.h numtype.h mathmisc.h random.h mathdef.h cell.h  suitesparse.h blas.h \
 matbin.h spbin.h zfarr.h chol.h  mat.h sp.h fft.h matmath.h loc.h locbin.h defs.h

EXTRA_DIST=suitesparse.c matmmap.c
EXTRA_DIST+=mat.c matmath.c blas.c matbin.c fft.c  matcomp.c amoeba.c bessik.c sp.c spmm.c spbin.c 

AOMATH_SRC= dmath.c cmath.c smath.c zmath.c lmath.c mathmisc.c random.c zfarr.c chol.c loc.c locbin.c cell.c

lib_LTLIBRARIES=libaomath.la
libaomath_la_SOURCES=$(AOMATH_SRC)
libaomath_la_LIBADD=../sys/libaosys.la $(LDCHOL)
# $(LDFFT) $(LAPACK) don't put fft or lapack here for flexibility when compiling for mex routines

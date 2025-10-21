AM_CFLAGS=-I$(top_srcdir)

nobase_include_HEADERS=blas.h defs.h fft.h loc.h matbin.h mathdef.h matmath.h random.h spbin.h type.h zmath.h\
 cell.h cmath.h dmath.h lmath.h map.h mat.h mathmisc.h numtype.h smath.h sp.h zfarr.h 

EXTRA_DIST=suitesparse.c matmmap.c 
EXTRA_DIST+=vec.c mat.c matmath.c blas.c matbin.c fft.c matcomp.c amoeba.c bessik.c sp.c spmm.c spbin.c 

AOMATH_SRC= dmath.c cmath.c smath.c zmath.c lmath.c mathmisc.c random.c zfarr.c loc.c map.c cell.c

lib_LTLIBRARIES=libaomath.la
libaomath_la_SOURCES=$(AOMATH_SRC)
libaomath_la_LIBADD=../sys/libaosys.la 
# $(FFTW_LIB) $(LAPACK) don't put fft or lapack here for flexibility when compiling for mex routines

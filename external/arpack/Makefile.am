#Specifying per object CFLAGS can avoid name conflict.
AM_FFLAGS=-DNDEBUG=1 -O3 -DUSE_MEM=0 -DHAVE_SSE2=1
#AM_CFLAGS=-DNDEBUG=1 -O3 -DUSE_MEM=0 -DHAVE_SSE2=1 -fPIC 
#AM_LDFLAGS=-static-libgfortran #override the value provided by configure
#LIBS= #while LIBS contains fftw3 by default
#if DO_SO
lib_LTLIBRARIES=libarpack.la
libarpack_la_SOURCES=$(ARPACK_SRC)
libarpack_la_FFLAGS=$(AM_FFLAGS)
#libarpack_la_CFLAGS=$(AM_CFLAGS)
libarpack_la_LDFLAGS=-static-libgfortran
#else
#lib_LIBRARIES=libarpack.a
#libarpack_a_SOURCES=$(ARPACK_SRC)
#libarpack_a_FFLAGS=$(AM_FFLAGS)
#libarpack_a_CFLAGS=$(AM_CFLAGS)
#endif
EXTRA_DIST=SRC/debug.h\
SRC/stat.h\
SRC/version.h
ARPACK_SRC=./SRC/dnaupd.f\
./SRC/dneigh.f\
./SRC/dlaqrb.f\
./SRC/dneupd.f\
./SRC/dsaupd.f\
./SRC/dseupd.f\
./SRC/dsaitr.f\
./SRC/dgetv0.f\
./SRC/dngets.f\
./SRC/dstqrb.f\
./SRC/dnaup2.f\
./SRC/dsesrt.f\
./SRC/dsortc.f\
./SRC/dsapps.f\
./SRC/dnconv.f\
./SRC/dnaitr.f\
./SRC/dsortr.f\
./SRC/dsgets.f\
./SRC/dstats.f\
./SRC/dstatn.f\
./SRC/dsaup2.f\
./SRC/dsconv.f\
./SRC/dnapps.f\
./SRC/dnaupe.f\
./SRC/dseigt.f\
./UTIL/dmout.f\
./UTIL/iswap.f\
./UTIL/dvout.f\
./UTIL/iset.f\
./UTIL/icopy.f\
./UTIL/icnteq.f\
./UTIL/ivout.f
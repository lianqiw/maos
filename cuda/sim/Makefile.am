AM_CFLAGS=-I $(top_srcdir)
AM_LDFLAGS=$(CUDA_LDFLAGS)
#In the order of include:
HEADER_LIB=cudata.h accphi.h gpu_sim.h  perf.h  wfs.h

EXTRA_DIST=$(HEADER_LIB) 
SRCS=atm.cu accphi.cu cudata.cu  perfevl.cu  perfinit.cu  pywfs.cu  wfsgrad.cu  wfsinit.cu  wfsints.cu

if USE_CUDA
lib_LTLIBRARIES=libcusim.la
libcusim_la_SOURCES=$(SRCS) $(EXTRA_DIST)
libcusim_la_LIBADD=../math/libcumath.la # $(LIBCUDA)
endif

include $(top_srcdir)/cuda/cuda.am

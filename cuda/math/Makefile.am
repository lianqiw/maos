AM_CFLAGS=-I $(top_srcdir)
AM_LDFLAGS=$(CUDA_LDFLAGS)
#In the order of include:
HEADER_LIB=common.h cublas.h curmat.h gpu_math.h kernel.h types.h utils.h prop_map.h cumath.h

EXTRA_DIST=$(HEADER_LIB) 
SRCS=cublas.cu curmat.cu kernel.cu types.cu utils.cu prop_map.cu gpu_math.cu

if USE_CUDA
lib_LTLIBRARIES=libcumath.la
libcumath_la_SOURCES=$(SRCS) $(EXTRA_DIST)
libcumath_la_LIBADD=$(CUDA_LIBS)
endif

include $(top_srcdir)/cuda/cuda.am

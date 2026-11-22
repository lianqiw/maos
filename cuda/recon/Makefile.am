AM_CFLAGS=-I $(top_srcdir)
AM_LDFLAGS=$(CUDA_LDFLAGS)
#In the order of include:
HEADER_LIB=fdpcg.h  fit.h  gpu_recon.h  pcg.h  recon.h  recon_geom.h  solve.h  tomo.h moao.h
#cuwfs.cu cuwfs.h
EXTRA_DIST=$(HEADER_LIB)
SRCS=fdpcg.cu  fit.cu  mvm_direct.cu  mvm_trans.cu  pcg.cu  recon.cu  recon_geom.cu  solve.cu  tomo.cu moao.cu

if USE_CUDA
lib_LTLIBRARIES=libcurecon.la
libcurecon_la_SOURCES=$(SRCS) $(EXTRA_DIST)
libcurecon_la_LIBADD=../sim/libcusim.la
endif

include $(top_srcdir)/cuda/cuda.am

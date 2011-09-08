#use -g -G for debugging.
AM_CFLAGS=-I/opt/cuda/include --compiler-options -Wall -g -O -I/opt/cuda/sdk/C/common/inc/ -I/opt/cuda/sdk/shared/inc/ -arch=sm_20 
if USE_LONG
AM_CFLAGS+=-DDLONG
endif
AM_CFLAGS+=-DHAVE_CONFIG_H -I..
AM_LDFLAGS=-L/opt/cuda/lib64/ -lcudart -lcurand 

EXTRA_DIST=gpu.h wfs.h utils.h accphi.h recon.h 
#LTLIBRARIES does not work.
lib_LIBRARIES=libgpu.a
libgpu_a_SOURCES=accphi.cu utils.cu curmat.cu cucmat.cu curmat.h cucmat.h wfsgrad.cu wfsints.cu wfsinit.cu perfevl.cu recon.cu pcg.cu fit.cu tomo.cu $(EXTRA_DIST)
libgpu_a_LIBADD=

.SUFFIXES=.cu
.cu.o:
	nvcc -o $@ -c $< $(AM_CFLAGS)

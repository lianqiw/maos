#use -g -G for debugging.
AM_CFLAGS=-I/opt/cuda/include --compiler-options -Wall -g -O3 -I/opt/cuda/sdk/C/common/inc/ -I/opt/cuda/sdk/shared/inc/ -arch=sm_20 --ptxas-options=-v 
if USE_LONG
AM_CFLAGS+=-DDLONG
endif
AM_CFLAGS+=-DHAVE_CONFIG_H -I..
AM_LDFLAGS=-L/opt/cuda/lib64/ -lcudart -lcurand 

EXTRA_DIST=gpu.h
#LTLIBRARIES does not work.
lib_LIBRARIES=libgpu.a
libgpu_a_SOURCES=accphi.cu utils.cu wfsgrad.cu perfevl.cu
libgpu_a_LIBADD=

.SUFFIXES=.cu.o
.cu.o:
	nvcc -o $@ -c $< $(AM_CFLAGS)

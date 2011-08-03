CC=cc
#use -g -G for debugging.
AM_CFLAGS=-I/opt/cuda/include --compiler-options -Wall -g -O3 -I/opt/cuda/sdk/C/common/inc/ -I/opt/cuda/sdk/shared/inc/ -arch=sm_20 --ptxas-options=-v 
AM_LDFLAGS=-L/opt/cuda/lib64/ -lcudart
#AM_LDFLAGS=-lcudart
#CUDASRC=test.cu

EXTRA_DIST=accphi.cu gpu.h

lib_LIBRARIES=libgpu.a
libgpu_a_SOURCES=gpu.h
libgpu_a_LIBADD=accphi.o

#bin_PROGRAMS=vectorAdd prop
#vectorAdd_SOURCES=
#vectorAdd_LDADD=vectorAdd.o 
#prop_SOURCES=
#prop_LDADD=prop.o

.cu.o:gpu.h
	nvcc -o $@ -c $< $(AM_CFLAGS) 

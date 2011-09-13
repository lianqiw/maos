#use -g -G for debugging.
AM_CFLAGS=-I/opt/cuda/include -arch=sm_20 --compiler-options -Wall -g -O 
if USE_LONG
AM_CFLAGS+=-DDLONG
endif
AM_CFLAGS+=-DHAVE_CONFIG_H -I..
AM_LDFLAGS=-L/opt/cuda/lib64/ -lcudart -lcurand 

EXTRA_DIST=accphi.h  cucmat.h  curmat.h  gpu.h  pcg.h  recon.h  types.h  utils.h  wfs.h
#LTLIBRARIES does not work.
lib_LIBRARIES=libgpu.a
libgpu_a_SOURCES=accphi.cu utils.cu curmat.cu cucmat.cu  wfsgrad.cu wfsints.cu wfsinit.cu perfevl.cu recon.cu pcg.cu fit.cu tomo.cu $(EXTRA_DIST)
libgpu_a_LIBADD=

.SUFFIXES=.cu
.cu.o:
	nvcc $(AM_CFLAGS) -c -o $@ $< 

#In Makefile for .c, there are include $(DEPDIR)/file.Po automatically generated. I don't know how to do for nvcc. So the dependency tracing is not correct.
#include ./$(DEPDIR)/accphi.Po
#The following rules takes care of the header dependency correctly.

#.cu.o:
#	nvcc $(AM_CFLAGS) -Xcompiler " -MT $@ -MD -MP -MF $(DEPDIR)/$*.Tpo" -c -o $@ $<
#	$(AM_V_at)mv -f $(DEPDIR)/$*.Tpo $(DEPDIR)/$*.Po




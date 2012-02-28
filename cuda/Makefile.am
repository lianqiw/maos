#use -g -G for debugging.
WARN="-Wall -Wno-unused-value -Wno-unused-function"
AM_CFLAGS=-arch=sm_$(CUDAVER) --compiler-options $(WARN) -g -O $(CUDAOPT)
if USE_LONG
AM_CFLAGS+=-DDLONG
endif
AM_CFLAGS+=-DHAVE_CONFIG_H -I..
AM_LDFLAGS=$(LDCUDA) -static

EXTRA_DIST=accphi.h  cucmat.h  curmat.h  gpu.h  pcg.h  recon.h  types.h  utils.h  wfs.h kernel.h cumat.h 
SRCS=accphi.cu utils.cu curmat.cu cucmat.cu  wfsgrad.cu wfsints.cu wfsinit.cu perfevl.cu recon.cu pcg.cu tomo.cu fit.cu moao.cu kernel.cu 
#LTLIBRARIES does not work.
lib_LIBRARIES=libgpu.a
libgpu_a_SOURCES=$(SRCS) $(EXTRA_DIST)
libgpu_a_LIBADD=

#.SUFFIXES=.cu

#The following does the depency generation/tracking and verbosity right.
CC=$(NVCC)

.cu.o:
	$(AM_V_CC)$(COMPILE) -c -o $@ $<
	@mkdir -p $(DEPDIR)
	@$(NVCC) $(AM_CFLAGS) $(CFLAGS) -M $< > $(DEPDIR)/$*.Po

-include $(SRCS:%.cu=$(DEPDIR)/%.Po)

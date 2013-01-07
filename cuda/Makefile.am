#use -g -G for debugging.
WARN=--compiler-options -Wall --compiler-options -Wno-unused-value --compiler-options -Wno-unused-function
VER=-gencode arch=compute_20,code=sm_20 -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35
#VER=-gencode arch=compute_35,code=sm_35
VER2=-arch=sm_$(CUDAVER)
AM_CFLAGS=$(VER) $(WARN) -g -O $(CUDAOPT)
AM_CFLAGS+=-Xptxas="-v" #for verbose 
#AM_CFLAGS+=-Xptxas -dlcm=cg
AM_CFLAGS+=-lineinfo
if USE_LONG
AM_CFLAGS+=-DDLONG
endif
AM_CFLAGS+=-DHAVE_CONFIG_H -I..
#AM_LDFLAGS=$(LDCUDA) -static

EXTRA_DIST=accphi.h  common.h  cucmat.h  cumat.h  curmat.h  gpu.h  kernel.h  pcg.h  recon.h  types.h  utils.h  wfs.h
SRCS=accphi.cu utils.cu curmat.cu cucmat.cu  wfsgrad.cu wfsints.cu wfsinit.cu perfevl.cu recon.cu pcg.cu tomo.cu fdpcg.cu fit.cu moao.cu kernel.cu cuchol.cu cumuv.cu prop_grid.cu mvm.cu test.cu mvm_trans.cu mvm_direct.cu mvm_iwfs.cu mvmfull_iwfs.cu mvm_test.cu

#LTLIBRARIES does not work.
lib_LIBRARIES=libgpu.a
libgpu_a_SOURCES=$(SRCS) $(EXTRA_DIST)
libgpu_a_LIBADD=

#The following does the depency generation/tracking and verbosity right.
CC=$(NVCC)

.cu.o:
	$(AM_V_CC)$(COMPILE) -c -o $@ $<
	@mkdir -p $(DEPDIR)
	@$(NVCC) $(AM_CFLAGS) $(CFLAGS) -M $< > $(DEPDIR)/$*.Po

-include $(SRCS:%.cu=$(DEPDIR)/%.Po)

#use -g -G for debugging.
WARN=--compiler-options -Wall --compiler-options -Wno-unused-value --compiler-options -Wno-unused-function
#VER=-gencode arch=compute_$(CUDAVER),code=sm_$(CUDAVER)
#VER+=-gencode arch=compute_20,code=sm_20 -gencode arch=compute_30,code=sm_30
#VER+=-gencode arch=compute_35,code=sm_35.
AM_CFLAGS=$(WARN) -g -O2 $(CUDAOPT)
#AM_CFLAGS+=-Xptxas="-v" #for verbose 
#AM_CFLAGS+=-Xptxas -dlcm=cg
AM_CFLAGS+=-lineinfo
if USE_LONG
AM_CFLAGS+=-DDLONG
endif
AM_CFLAGS+=-DHAVE_CONFIG_H -I..
#AM_LDFLAGS=$(LDCUDA) -static

EXTRA_DIST=accphi.h  common.h  cucmat.h  cumat.h  curmat.h  gpu.h  kernel.h  pcg.h  recon.h  types.h  utils.h  wfs.h recon_geom.h fit.h tomo.h solve.h fdpcg.h prop_wrap.h cudata.h
SRCS=solve.cu accphi.cu utils.cu curmat.cu cucmat.cu  wfsgrad.cu wfsints.cu wfsinit.cu perfevl.cu recon.cu pcg.cu recon_geom.cu tomo.cu fdpcg.cu fit.cu moao.cu kernel.cu prop_grid.cu prop_wrap.cu mvm.cu test.cu mvm_trans.cu mvm_direct.cu mvm_iwfs.cu mvmfull_iwfs.cu mvmfull_pipe.cu cudata.cu

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

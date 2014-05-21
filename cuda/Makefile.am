WARN=-Xcompiler -Wall -Xcompiler -Wno-unused-value -Xcompiler 
AM_CFLAGS=$(WARN) $(CUDAOPT) 
if DEBUG
AM_CFLAGS+=-O0 -g
else
AM_CFLAGS+=-O3
endif
AM_CFLAGS+=-Xptxas="-v" #for verbose 
#AM_CFLAGS+=-Xptxas -dlcm=cg
AM_CFLAGS+=-lineinfo
#AM_CFLAGS+=-Xcompiler -arch -Xcompiler x86_64
if USE_LONG
AM_CFLAGS+=-DDLONG
endif
AM_CFLAGS+=-DHAVE_CONFIG_H -I..
#AM_LDFLAGS=$(LDCUDA) -static
#In the order of include:
HEADER_LIB=common.h types.h cucmat.h  cumat.h  curmat.h gpu.h kernel.h utils.h prop_wrap.h
HEADER_RECON=solve.h recon_base.h pcg.h tomo.h fit.h recon.h fdpcg.h moao.h
#cuwfs.cu cuwfs.h
EXTRA_DIST=$(HEADER_LIB) $(HEADER_RECON) accphi.h wfs.h  cudata.h perf.h
SRC_LIB=utils.cu curmat.cu cucmat.cu kernel.cu prop_grid.cu prop_wrap.cu accphi.cu cudata.cu 
SRC_RECON=solve.cu recon_base.cu pcg.cu tomo.cu fit.cu recon.cu fdpcg.cu moao.cu mvm_trans.cu mvm_direct.cu 
SRC_MVM= mvm_iwfs.cu mvmfull_iwfs.cu mvmfull_real.cu mvmfull_pipe.cu mvm_daemon.cu
SRC_WFS= wfsinit.cu wfsints.cu wfsgrad.cu
SRCS=$(SRC_LIB) $(SRC_RECON) $(SRC_WFS) $(SRC_MVM) perfinit.cu perfevl.cu test.cu
#LTLIBRARIES does not work.
lib_LIBRARIES=libgpu.a
libgpu_a_SOURCES=$(SRCS) $(EXTRA_DIST)
libgpu_a_LIBADD=

#The following does the depency generation/tracking and verbosity right.

.cu.o:
	$(AM_V_CC)$(NVCC) $(AM_CFLAGS) $(CFLAGS) -c -o $@ $<
	@mkdir -p $(DEPDIR)
	@$(NVCC) $(AM_CFLAGS) $(CFLAGS) -M $< > $(DEPDIR)/$*.Po

-include $(SRCS:%.cu=$(DEPDIR)/%.Po)

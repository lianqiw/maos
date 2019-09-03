AM_CFLAGS=$(OPTS) -I $(top_srcdir)
AM_LDFLAGS=$(LDOPTS) 
#In the order of include:
HEADER_LIB=common.h types.h cucmat.h  Array.h  curmat.h gpu.h kernel.h utils.h prop_map.h
HEADER_RECON=solve.h recon_geom.h pcg.h tomo.h fit.h recon.h fdpcg.h moao.h
#cuwfs.cu cuwfs.h
EXTRA_DIST=$(HEADER_LIB) $(HEADER_RECON) accphi.h wfs.h  cudata.h perf.h
SRC_LIB=utils.cu curmat.cu cucmat.cu kernel.cu prop_map.cu accphi.cu cudata.cu atm.cu
SRC_RECON=solve.cu recon_geom.cu tomo.cu fit.cu moao.cu fdpcg.cu pcg.cu recon.cu mvm_trans.cu mvm_direct.cu 
SRC_MVM= mvm_iwfs.cu mvmfull_iwfs.cu mvmfull_real.cu mvmfull_pipe.cu mvm_daemon.cu
SRC_WFS= wfsinit.cu wfsints.cu wfsgrad.cu pywfs.cu
SRCS=$(SRC_LIB) $(SRC_RECON) $(SRC_WFS) $(SRC_MVM) perfinit.cu perfevl.cu 

#must include EXTRA_DIST in SOURCES, otherwise some automake tools are missing

#lib_LIBRARIES=libaogpu.a
#libaogpu_a_SOURCES=$(SRCS) $(EXTRA_DIST)
#libaogpu_a_LDFLAGS=$(LDOPTS)

#else
lib_LTLIBRARIES=libaogpu.la
libaogpu_la_SOURCES=$(SRCS) $(EXTRA_DIST)
libaogpu_la_LIBADD=../lib/libaos.la $(LIBCUDA)
libaogpu_la_LDFLAGS=$(LDOPTS)
#endif
#if USE_STATIC
#all-local:dlink.o
#dlink.o:libaogpu.la
#	$(NVCC) $(CUDAOPT) -dlink .libs/libaogpu.a -o dlink.o
#endif
LDADD=libaogpu.la $(LIBCUDA) $(LDEXE_FULL)
#bin_PROGRAMS=../test/test_cuda ../bin/gpu_mvm
#___test_test_cuda_SOURCES=test.cu 
#___bin_gpu_mvm_SOURCES=gpu_mvm.c

#The following does the depency generation/tracking and verbosity right.
.cu.o:
	$(AM_V_CC)$(NVCC) $(CUDAOPT) -Xcompiler "$(CCBINOPT)" -c -o $@ $<
	@mkdir -p $(DEPDIR)
	@$(NVCC) $(CUDAOPT) -M $< -MT $@ > $(DEPDIR)/$*.Po

.cu.lo:
	$(AM_V_CC)$(LIBTOOL) $(AM_V_lt) --tag=CC $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=compile $(NVCC) $(CUDAOPT) --compiler-options "$(CCBINOPT) " -c $< -o $@ 
	@mkdir -p $(DEPDIR)
	@$(NVCC) $(CUDAOPT) -M $< -MT $@ > $(DEPDIR)/$*.Po

-include $(SRCS:%.cu=$(DEPDIR)/%.Po)


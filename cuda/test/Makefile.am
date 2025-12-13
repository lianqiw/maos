AM_CFLAGS=-I $(top_srcdir)
AM_LDFLAGS=$(CUDA_LDLAGS)
LDADD=../recon/libcurecon.la ../../lib/libaos.la $(EXEFULL_LIBS)
#In the order of include:
HEADER_LIB=gpu_test.h
#cuwfs.cu cuwfs.h
EXTRA_DIST=$(HEADER_LIB)

if USE_CUDA
bin_PROGRAMS=gpu_mvm_daemon gpu_mvm_test gpu_test gpu_ext
gpu_mvm_daemon_SOURCES=gpu_mvm_daemon.c mvm_daemon.cu 
gpu_mvm_test_SOURCES=gpu_mvm_test.c mvm_iwfs.cu mvm_test.cu mvmfull_iwfs.cu mvmfull_pipe.cu mvmfull_real.cu
gpu_test_SOURCES=gpu_test.cu
TESTS=gpu_ext
endif
#The following does the depency generation/tracking and verbosity right.
.cu.o:
	$(AM_V_CC)$(NVCC) $(CUDA_CFLAGS) -c -o $@ $<
	@mkdir -p $(DEPDIR)
	@$(NVCC) $(CUDA_CFLAGS) -M $< -MT $@ > $(DEPDIR)/$*.Po

.cu.lo:
	$(AM_V_CC)$(LIBTOOL) $(AM_V_lt) --tag=CC $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=compile $(NVCC) $(CUDA_CFLAGS) -c $< -o $@ 
	@mkdir -p $(DEPDIR)
	@$(NVCC) $(CUDA_CFLAGS) -M $< -MT $@ > $(DEPDIR)/$*.Po

-include $(SRCS:%.cu=$(DEPDIR)/%.Po)

clean-local:
	-rm -rf $(DEPDIR)/*.Po

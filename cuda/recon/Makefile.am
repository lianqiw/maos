AM_CFLAGS=$(OPTS) -I $(top_srcdir)
AM_LDFLAGS=$(LDOPTS) $(LDCUDA)
#In the order of include:
HEADER_LIB=fdpcg.h  fit.h  gpu_recon.h  pcg.h  recon.h  recon_geom.h  solve.h  tomo.h moao.h
#cuwfs.cu cuwfs.h
EXTRA_DIST=$(HEADER_LIB)
SRCS=fdpcg.cu  fit.cu  mvm_direct.cu  mvm_trans.cu  pcg.cu  recon.cu  recon_geom.cu  solve.cu  tomo.cu moao.cu

if USE_CUDA
lib_LTLIBRARIES=libcurecon.la
libcurecon_la_SOURCES=$(SRCS) $(EXTRA_DIST)
libcurecon_la_LIBADD=../sim/libcusim.la
libcurecon_la_LDFLAGS=$(LDOPTS)
endif

#The following does the depency generation/tracking and verbosity right.
.cu.o:
	$(AM_V_CC)$(NVCC) $(CUDAOPT) $(CCBINOPT) -c -o $@ $<
	@mkdir -p $(DEPDIR)
	@$(NVCC) $(CUDAOPT) -M $< -MT $@ > $(DEPDIR)/$*.Po

.cu.lo:
	$(AM_V_CC)$(LIBTOOL) $(AM_V_lt) --tag=CC $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=compile $(NVCC) $(CUDAOPT) $(CCBINOPT)  -c $< -o $@ 
	@mkdir -p $(DEPDIR)
	@$(NVCC) $(CUDAOPT) -M $< -MT $@ > $(DEPDIR)/$*.Po

-include $(SRCS:%.cu=$(DEPDIR)/%.Po)

clean-local:
	-rm -rf $(DEPDIR)/*.Po

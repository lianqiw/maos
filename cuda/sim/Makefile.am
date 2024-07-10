AM_CFLAGS=-I $(top_srcdir)
AM_LDFLAGS=$(LDCUDA)
#In the order of include:
HEADER_LIB=cudata.h accphi.h gpu_sim.h  perf.h  wfs.h

EXTRA_DIST=$(HEADER_LIB) 
SRCS=atm.cu accphi.cu cudata.cu  perfevl.cu  perfinit.cu  pywfs.cu  wfsgrad.cu  wfsinit.cu  wfsints.cu

if USE_CUDA
lib_LTLIBRARIES=libcusim.la
libcusim_la_SOURCES=$(SRCS) $(EXTRA_DIST)
libcusim_la_LIBADD=../math/libcumath.la # $(LIBCUDA)
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

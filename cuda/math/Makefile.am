AM_CFLAGS=-I $(top_srcdir)
AM_LDFLAGS=$(LDCUDA)
#In the order of include:
HEADER_LIB=common.h cublas.h curmat.h gpu_math.h kernel.h types.h utils.h prop_map.h cumath.h

EXTRA_DIST=$(HEADER_LIB) 
SRCS=cublas.cu curmat.cu kernel.cu types.cu utils.cu prop_map.cu gpu_math.cu

if USE_CUDA
lib_LTLIBRARIES=libcumath.la
libcumath_la_SOURCES=$(SRCS) $(EXTRA_DIST)
libcumath_la_LIBADD=$(LIBCUDA)
endif

#The following does the depency generation/tracking and verbosity right.
.cu.o:
	$(AM_V_CC)$(NVCC) $(CUDAOPT) $(CCBINOPT) -c -o $@ $<
	@mkdir -p $(DEPDIR)
	@$(NVCC) $(CUDAOPT) -M $< -MT $@ > $(DEPDIR)/$*.Po

.cu.lo:
	$(AM_V_CC)$(LIBTOOL) $(AM_V_lt) --tag=CC $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=compile $(NVCC) $(CUDAOPT) $(CCBINOPT) -c $< -o $@ 
	@mkdir -p $(DEPDIR)
	@$(NVCC) $(CUDAOPT) -M $< -MT $@ > $(DEPDIR)/$*.Po

-include $(SRCS:%.cu=$(DEPDIR)/%.Po)

clean-local:
	-rm -rf $(DEPDIR)/*.Po

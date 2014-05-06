LDADD=../lib/libaos.la ../math/libaomath.la ../sys/libaosys.la $(LDOPTS)
AM_CFLAGS=$(OPTS)
AM_LDFLAGS=$(LDEXE) 
bin_PROGRAMS=test_mat test_plot test_draw test_rand test_sparse test_cmat2		\
		test_fft test_proj test_loc test_chol test_accphi convert_psfout	\
		test_fd test_atm test_fractal test_svdlas2 test_fits 
if USE_CUDA
bin_PROGRAMS+=test_gpu
test_gpu_LDADD=../cuda/libgpu.a $(LDADD)  $(LDCUDA)
endif

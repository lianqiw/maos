AM_CFLAGS=$(OPTS) -I $(top_srcdir) -g -O0
AM_LDFLAGS=$(LDOPTS) 
LDADD=../lib/libaos.la $(LDEXE_FULL) 

bin_PROGRAMS=test_mat test_plot test_draw test_rand test_sparse test_cmat2 \
	test_fft test_proj test_loc test_chol test_accphi	\
	test_fd test_atm test_fractal test_fits
#if USE_CUDA
#bin_PROGRAMS+=test_gpu
#test_gpu_LDADD=../cuda/libaogpu.la $(LDADD)  $(LIBCUDA)
#endif

if USE_CMOCKA
bin_PROGRAMS+=cmocka_mat cmocka_matmath cmocka_locfft cmocka_blas \
cmocka_readstr cmocka_threadpool cmocka_psd
LDADD+=-lcmocka
endif

TESTS=$(bin_PROGRAMS)

clean-local:
	@rm -rf *.bin *.fits
if GCOV	
	@rm -fr coverage
	@find . -name "*.gcda" -exec rm {} \;
	@lcov --directory . --zerocounters
endif

if GCOV
check-local:
	@mkdir -p coverage
	@lcov --compat-libtool --directory .. --capture --output-file coverage/app.info  
	@genhtml -o coverage/ coverage/app.info 
endif
AM_LDFLAGS=$(LDOPTS) 
LDADD=../lib/libaos.a 
AM_CFLAGS=$(OPTS)
bin_PROGRAMS=test_mat test_plot test_draw test_rand test_sparse		\
		test_cmat2 test_fft test_proj test_loc	\
		test_chol test_eigs test_accphi convert_psfout test_fd test_atm

test_chol_LDFLAGS=$(AM_LDFLAGS) 
test_chol_LDADD=$(LDADD) ../external/cholmod/libchol.a 
test_eigs_LDFLAGS=$(AM_LDFLAGS) $(LAPACK)
test_eigs_LDADD=$(LDADD) ../external/arpack/libarpack.a


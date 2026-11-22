AM_CFLAGS=-I $(top_srcdir) -g -O0
LDADD=../lib/libaos.la $(EXEFULL_LIBS)
EXTRA_DIST=unittest.h unittest_driver.sh unittest_report.sh unittest_coverage.sh \
	scheduler/run.sh scheduler/common.py scheduler/test_commands.py scheduler/test_auth.py scheduler/test_websocket.py

#All unit tests use cmocka. They must not depend on temporary input or output files.
if HAS_CMOCKA
bin_PROGRAMS=cmocka_sys cmocka_threadpool cmocka_readstr cmocka_fits \
	cmocka_mat cmocka_matmath cmocka_matops cmocka_cmat cmocka_fft cmocka_blas \
	cmocka_rand cmocka_sparse cmocka_chol cmocka_loc cmocka_locfft cmocka_psd \
	cmocka_accphi cmocka_proj cmocka_fdpcg cmocka_atm cmocka_fractal cmocka_rayleigh \
	cmocka_sock cmocka_sysutil cmocka_map cmocka_cell cmocka_mathmisc \
	cmocka_wfs cmocka_recon cmocka_modes cmocka_pywfs \
	cmocka_servo cmocka_cfg cmocka_libmisc cmocka_scheduler cmocka_shwfs cmocka_wfsutil
LDADD+=$(CMOCKA_LIBS)
AM_CFLAGS+=$(CMOCKA_CFLAGS)
endif
TESTS=$(bin_PROGRAMS)
#The driver also updates the unit test status report docs/unittest_status.md in the
#build folder after each test. make doc includes it in the documentation.
LOG_DRIVER=$(SHELL) $(srcdir)/unittest_driver.sh "$(top_srcdir)" "$(top_builddir)"

#Tests of the scheduler protocol (docs/94_scheduler.md). They start a scheduler
#and use the network, so they are not part of make check. Run make check-scheduler.
EXTRA_PROGRAMS=scheduler_probe
scheduler_probe_SOURCES=scheduler/probe.c
scheduler_probe_LDADD=../lib/libaos.la $(EXEFULL_LIBS)
if HAS_LWS
SCHEDULER_LWS=1
endif
check-scheduler: scheduler_probe$(EXEEXT)
	$(SHELL) $(srcdir)/scheduler/run.sh "$(abs_top_srcdir)" "$(abs_top_builddir)/bin" "$(abs_builddir)/scheduler_probe$(EXEEXT)" $(SCHEDULER_LWS)
.PHONY: check-scheduler

CLEANFILES=maos_check_GPU_5.all maos_check_GPU_5.res maos_check_GPU_5.cur \
	maos_check_CPU_5.all maos_check_CPU_5.res maos_check_CPU_5.cur
clean-local:
	@rm -rf *.bin *.fits unittest_report.lock
if GCOV
	@rm -fr coverage
	@find . -name "*.gcda" -exec rm {} \;
	@if command -v lcov >/dev/null ;then lcov --directory . --zerocounters ;fi
endif

if GCOV
check-local:
	@if command -v lcov >/dev/null && command -v genhtml >/dev/null ;then \
		mkdir -p coverage && \
		lcov --compat-libtool --directory .. --capture --output-file coverage/app.info && \
		genhtml -o coverage/ coverage/app.info ;\
	fi
endif

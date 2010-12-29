AM_CFLAGS=$(OPTS) -I$(top_srcdir) -I$(top_srcdir)/lib/sys -I$(top_srcdir)/external/cholmod/Include

nobase_include_HEADERS=accphi.h common.h type.h locbin.h misc.h mathmisc.h	\
random.h sp.h dsp.h csp.h bin.h cell.h mkg.h path.h readcfg.h suitesparse.h	\
blas.h matbin.h cellarr.h mkh.h aos.h chol.h mkw.h turbulence.h cmat_extra.h	\
mat.h dmat.h cmat.h loc.h mkz.h proj.h hashlittle.h laplacian.h pcg.h muv.h	\
draw.h fft.h genotf.h fractal.h stfun.h

nobase_include_HEADERS+=sys/mem.h sys/process.h sys/scheduler_server.h		\
sys/sockio.h sys/thread.h sys/io.h sys/daemonize.h sys/scheduler_client.h	\
sys/shm.h sys/sys.h sys/thread_pool.h

EXTRA_DIST=cell.c sp.c matbin.c mat.c suitesparse.c fractal_do.c

AOS_SRC=readcfg.c dmat.c cmat.c cmat_extra.c hashlittle.c dsp.c csp.c bin.c	\
misc.c mathmisc.c random.c loc.c locbin.c cellarr.c turbulence.c path.c mkg.c	\
mkh.c mkz.c mkw.c accphi.c proj.c chol.c laplacian.c pcg.c muv.c draw.c fft.c	\
genotf.c fractal.c stfun.c

AOS_SRC+=sys/io.c sys/daemonize.c sys/process_bsd.c sys/process_cygwin.c	\
sys/scheduler_server.c sys/sockio.c sys/mem.c sys/process_apple.c sys/process.c	\
sys/process_linux.c sys/scheduler_client.c sys/shm.c sys/thread.c sys/thread_pool.c

lib_LTLIBRARIES=libaos.la
libaos_la_SOURCES=$(AOS_SRC)
libaos_la_LIBADD=$(LDCHOL) $(LDFFT)
libaos_la_LDFLAGS = -avoid-version
#The following enables the generation of documentation for cmat, ccell, csp functions.
dmat.c:tmp.c
tmp.c:$(srcdir)/mat.c $(srcdir)/cell.c $(srcdir)/matbin.c $(srcdir)/sp.c
	@echo "#define USE_COMPLEX" > $(srcdir)/tmp.c
	@cat $(srcdir)/mat.c $(srcdir)/cell.c $(srcdir)/matbin.c $(srcdir)/sp.c >> $(srcdir)/tmp.c

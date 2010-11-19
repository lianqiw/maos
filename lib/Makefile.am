SUBDIRS=sys/ .
AM_CFLAGS=$(OPTS) -I$(top_srcdir)
if COMP_CHOL
AM_CFLAGS+=-I$(top_srcdir)/external/cholmod/Include
endif
#include_HEADERS
nobase_include_HEADERS=accphi.h common.h type.h locbin.h misc.h mathmisc.h	\
random.h sp.h dsp.h csp.h bin.h cell.h mkg.h path.h readcfg.h suitesparse.h	\
blas.h matbin.h cellarr.h mkh.h aos.h chol.h mkw.h turbulence.h cmat_extra.h	\
mat.h dmat.h cmat.h loc.h mkz.h proj.h hashlittle.h laplacian.h pcg.h muv.h	\
draw.h fft.h genotf.h

EXTRA_DIST=cell.c sp.c matbin.c mat.c suitesparse.c

AOS_SRC=readcfg.c dmat.c cmat.c cmat_extra.c hashlittle.c dsp.c csp.c bin.c	\
misc.c mathmisc.c random.c loc.c locbin.c cellarr.c turbulence.c path.c mkg.c	\
mkh.c mkz.c mkw.c accphi.c proj.c chol.c laplacian.c pcg.c muv.c draw.c	\
fft.c genotf.c


lib_LTLIBRARIES=libaos.la
libaos_la_SOURCES=$(AOS_SRC)
libaos_la_LIBADD=sys/libsys.la $(LDCHOL) $(LDFFT)

#The following enables the generation of documentation for cmat, ccell, csp functions.
dmat.c:tmp.c
tmp.c:$(srcdir)/mat.c $(srcdir)/cell.c $(srcdir)/matbin.c $(srcdir)/sp.c
	@echo "#define USE_COMPLEX" > $(srcdir)/tmp.c
	@cat $(srcdir)/mat.c $(srcdir)/cell.c $(srcdir)/matbin.c $(srcdir)/sp.c >> $(srcdir)/tmp.c
PROGS=
if DRAW
PROGS+=../bin/drawbin
___bin_drawbin_SOURCES=drawbin.c
drawbin_CFLAGS=$(OPTS)
___bin_drawbin_LDADD=libaos.la $(LDOPTS)
endif
bin_PROGRAMS=$(PROGS)
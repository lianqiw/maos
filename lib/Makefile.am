AM_CFLAGS=$(OPTS) -I$(top_srcdir) -I$(top_srcdir)/external/cholmod/Include
if !WIN
AM_CFLAGS+=-fPIC
endif
#if set LDCHOL here instead of in configure.ac compiles "No rule to make target ..."

nobase_include_HEADERS=accphi.h type.h locbin.h mathmisc.h random.h sp.h dsp.h	\
 ssp.h csp.h bin.h cell.h mkg.h readcfg.h readstr.h suitesparse.h blas.h	\
 matbin.h spbin.h cellarr.h mkh.h aos.h chol.h mkw.h turbulence.h mat.h dmat.h	\
 smat.h cmat.h zmat.h loc.h mkz.h proj.h laplacian.h pcg.h muv.h draw.h fft.h	\
 genotf.h fractal.h stfun.h servo.h slaving.h imat.h mkdtf.h psd.h

EXTRA_DIST=cell.c sp.c matbin.c matmmap.c spbin.c mat.c suitesparse.c		\
fractal_do.c prop_grid_stat.c prop_grid_pts.c prop_grid_map.c defs.h blas.c

AOS_SRC=readcfg.c readstr.c dmat.c smat.c zmat.c cmat.c dsp.c	\
ssp.c csp.c bin.c mathmisc.c random.c loc.c locbin.c cellarr.c		\
turbulence.c mkg.c mkh.c mkz.c mkw.c accphi.c proj.c chol.c laplacian.c	\
pcg.c muv.c draw.c fft.c genotf.c fractal.c stfun.c servo.c slaving.c imat.c	\
bessik.c amoeba.c mkdtf.c psd.c


if USE_LONG
lib_LTLIBRARIES=libaos_long.la
libaos_long_la_SOURCES=$(AOS_SRC)
libaos_long_la_LIBADD=$(LDLIB)
libaos_long_la_LDFLAGS = -avoid-version -static
else
lib_LTLIBRARIES=libaos_int.la
libaos_int_la_SOURCES=$(AOS_SRC)
libaos_int_la_LIBADD=$(LDLIB)
libaos_int_la_LDFLAGS = -avoid-version -static
endif

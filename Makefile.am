ACLOCAL_AMFLAGS=-I m4 --install

SUBDIRS=bin test
if USE_MEX
SUBDIRS+=mex
endif
# include doxygen definition.
include $(top_srcdir)/aminclude.am

#dist-hook:
#	rm -rf $(distdir)/configure.ac 


#all-local:$(srcdir)/math/dtmp.c
#For documentation only. | enforces ordering in parallel build
.NOTPARALLEL:doc
doc:| tmp clear-doc doxygen-run clear-tmp
#The following tmp files are used to by pass two doxygen limitations
#1. files included in c files are only expanded for dependencies, not original definitions
#2. macros in C files are only expanded once.
tmp:
	@cd $(srcdir)/math && cat defs.h mat.c matmath.c blas.c matbin.c matmmap.c fft.c amoeba.c bessik.c sp.c spmm.c spbin.c > dtmp.c \
	&& echo '#define COMP_COMPLEX\n' > ctmp.c && cat dtmp.c | sed 's/X(/CX(/g' >>ctmp.c \
	&& echo '#define COMP_SINGLE\n' > stmp.c && cat dtmp.c | sed 's/X(/SX(/g'>> stmp.c \
	&& echo '#define COMP_COMPLEX\n#define COMP_SINGLE\n' > ztmp.c && cat dtmp.c | sed 's/X(/ZX(/g' >> ztmp.c \
	&& echo '#define COMP_LONG\n' >ltmp.c && cat defs.h mat.c matbin.c matmmap.c | sed 's/X(/LX(/g' >> ltmp.c

clear-tmp:
	@cd $(srcdir)/math && rm dtmp.c ctmp.c stmp.c ztmp.c ltmp.c

clear-doc:
	@rm -rf doxygen-doc

EXTRA_DIST = doxygen.cfg docs scripts config 
BUILT_SOURCES=version.h #non existent file will always run
if LIBTOOL_FIX_OPENMP
BUILT_SOURCES+=fixlibtool
endif
fixlibtool:libtool
	@sed -i '' 's/|-fopenmp|/|/g' libtool
version.h:FORCE #always runs
	@str="static const char *GIT_VERSION=\"`cd $(srcdir) && git show -s --format='%h %ci' || echo unknown`\";" && \
	if test ! -f version.h -o -n "`echo $$str | diff version.h -`" ;then echo $$str >version.h;fi
FORCE:
.ONESHELL:

CLEANFILES=version.h

MAKEOVERRIDES=--no-print-directory

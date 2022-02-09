ACLOCAL_AMFLAGS=-I m4 --install
AM_MAKEFLAGS=--no-print-directory

SUBDIRS=sys math lib cuda maos skyc tools test mex

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
BUILT_SOURCES=version.h bin/Makefile
if LIBTOOL_FIX_OPENMP
BUILT_SOURCES+=fixlibtool
endif
fixlibtool:libtool
	@sed 's/|-fopenmp|/|/g' libtool > libtool.new
	@mv libtool.new libtool
version.h:FORCE #non existent target will always run
	@str="static const char *GIT_VERSION=\"`cd $(srcdir) && git show -s --format='%h %ci' || echo $(VERSION)`\";" && \
	if test ! -f version.h -o -n "`echo $$str | diff version.h -`" ;then echo $$str >version.h;fi
FORCE:
bin/Makefile:Makefile
	@mkdir -p bin/
	@echo -e "all:\n\t@+make --no-print-directory -C .." > bin/Makefile
.ONESHELL:

CLEANFILES=version.h bin/Makefile

clean-local:
	-rm -rf external

#- before command means ignore error
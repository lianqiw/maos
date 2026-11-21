ACLOCAL_AMFLAGS=-I m4 --install
AM_MAKEFLAGS=--no-print-directory
AM_DISTCHECK_CONFIGURE_FLAGS=--enable-release 
SUBDIRS=sys math lib cuda maos skyc tools test mex

# include doxygen definition.
include $(top_srcdir)/aminclude.am

#dist-hook:
#	rm -rf $(distdir)/configure.ac 


#all-local:$(srcdir)/math/dtmp.c
#For documentation only. | enforces ordering in parallel build
.NOTPARALLEL:doc
doc:| tmp clear-doc unittest-doc doxygen-run clear-tmp
#Run the unit tests and generate the unit test status and coverage report
#docs/unittest_status.md in the build folder. It replaces the placeholder
#docs/96_unittest.md of the source folder in the documentation.
unittest-doc:
if GCOV
	@for dir in sys math lib test; do find $$dir -name "*.gcda" -exec rm -f {} \; ;done
endif
	-@$(MAKE) $(AM_MAKEFLAGS) check
	@$(SHELL) $(srcdir)/test/unittest_report.sh $(abs_top_srcdir) $(abs_top_builddir) --coverage
DX_ENV+=UNITTEST_DOC='$(abs_top_builddir)/docs/unittest_status.md' UNITTEST_EXCLUDE='$(abs_top_srcdir)/docs/96_unittest.md'
#The following tmp files are used to by pass two doxygen limitations
#1. files included in c files are only expanded for dependencies, not original definitions
#2. macros in C files are only expanded once.
tmp:
	@cd $(srcdir)/math && cat defs.h mat.c matmath.c blas.c matbin.c matmmap.c fft.c amoeba.c bessik.c sp.c spmm.c spbin.c > dtmp.c \
	&& printf '#define COMP_COMPLEX\n' > ctmp.c && cat dtmp.c | sed 's/X(/CX(/g' >>ctmp.c \
	&& printf '#define COMP_SINGLE\n' > stmp.c && cat dtmp.c | sed 's/X(/SX(/g'>> stmp.c \
	&& printf '#define COMP_COMPLEX\n#define COMP_SINGLE\n' > ztmp.c && cat dtmp.c | sed 's/X(/ZX(/g' >> ztmp.c \
	&& printf '#define COMP_LONG\n' >ltmp.c && cat defs.h mat.c matbin.c matmmap.c | sed 's/X(/LX(/g' >> ltmp.c

clear-tmp:
	@cd $(srcdir)/math && rm dtmp.c ctmp.c stmp.c ztmp.c ltmp.c

clear-doc:
	@rm -rf doxygen-doc docs/unittest_status.md

EXTRA_DIST = doxygen.cfg docs scripts config 
BUILT_SOURCES=bin/Makefile
if LIBTOOL_FIX_OPENMP
BUILT_SOURCES+=fixlibtool
endif
fixlibtool:libtool
	@sed 's/|-fopenmp|/|/g' libtool > libtool.new
	@mv libtool.new libtool

bin/Makefile:Makefile
	@mkdir -p bin/
	@printf "all:\n\t@+make --no-print-directory -C .." > bin/Makefile
.ONESHELL:

CLEANFILES=bin/Makefile

clean-local:
	-rm -rf external

#- before command means ignore error

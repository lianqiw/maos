ACLOCAL_AMFLAGS=-I m4 --install
SUBDIRS=
SUBDIRS+=sys math lib
if USE_CUDA
SUBDIRS+=cuda
endif
SUBDIRS+=maos skyc
if USE_MEX
SUBDIRS+=mex
endif
SUBDIRS+=tools test 
include $(top_srcdir)/aminclude.am

#dist-hook:
#	rm -rf $(distdir)/configure.ac 

all-local:$(srcdir)/math/tmp.c
$(srcdir)/math/tmp.c:$(srcdir)/math/mat.c $(srcdir)/math/cell.c $(srcdir)/math/matbin.c \
	$(srcdir)/math/spbin.c $(srcdir)/math/sp.c
	@echo "#define USE_COMPLEX" > $(srcdir)/math/tmp.c
	@cat $(srcdir)/math/mat.c $(srcdir)/math/cell.c $(srcdir)/math/matbin.c $(srcdir)/math/spbin.c $(srcdir)/math/sp.c >> $(srcdir)/math/tmp.c

EXTRA_DIST = doxygen.cfg docs version.h scripts config

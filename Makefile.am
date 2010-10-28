ACLOCAL_AMFLAGS=-I m4 --install
SUBDIRS=external/ lib/ maos/ skyc/ config/

include $(top_srcdir)/aminclude.am

EXTRA_DIST = doxygen.cfg doxygen_layout.xml 

#enable GTK-Doc for the distcheckmake target
#DISTCHECK_CONFIGURE_FLAGS = --enable-gtk-doc

#Note: in EXTRA_DIST, don't specify a folder. If you do, everything in that
#folder, including compiled files (in the compiled folder), svn directory, will
#be included.

EXTRA_DIST += mex/Makefile mex/write.c mex/read.c mex/io.c mex/io.h

EXTRA_DIST += test/Makefile.am test/Makefile.in test/convert_psfout.c	\
test/test_cmat2.c test/test_fft.c test/test_proj.c test/test_tomo.c	\
test/test_accphi.c test/test_draw.c test/test_loc.c test/test_rand.c	\
test/test_atm.c test/test_eigs.c test/test_mat.c test/test_sanea.c	\
test/test_chol.c test/test_fd.c test/test_plot.c test/test_sparse.c

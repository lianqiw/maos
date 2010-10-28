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

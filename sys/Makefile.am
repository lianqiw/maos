AM_CFLAGS=-I$(top_srcdir) 

nobase_include_HEADERS=mem.h process.h scheduler.h sock.h sockio.h thread.h	\
 daemonize.h scheduler_client.h sys.h thread_pool.h misc.h hashlittle.h path.h	\
 common.h bin.h readstr.h
EXTRA_DIST=process_apple.c process_linux.c process_bsd.c process_cygwin.c
SYS_SRC=daemonize.c misc.c process.c thread.c hashlittle.c mem.c path.c	\
 scheduler_client.c sockio.c thread_pool.c sock.c bin.c readstr.c version.c


lib_LTLIBRARIES=libaosys.la
libaosys_la_SOURCES=$(SYS_SRC)
libaosys_la_LIBADD=-lz -lm -ldl #don't remove
#we put -lz -lm -ldl here so that libtool will always put these flags after the library when linking exes.

BUILT_SOURCES=version.c
version.c:FORCE #non existent target will always run
	@str="const char *GIT_VERSION=\"`cd $(srcdir) && git show -s --format='%h %ci' || (echo $(VERSION) && date +'%F %r %z') `\";" && \
	if test ! -f version.c -o -n "`echo $$str | diff version.c -`" ;then echo $$str >version.c;fi
FORCE:	
CLEANFILES=version.c

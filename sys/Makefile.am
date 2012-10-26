AM_CFLAGS=$(OPTS) -I$(top_srcdir) 
if !WIN
AM_CFLAGS+=-fPIC
endif

nobase_include_HEADERS=mem.h process.h scheduler_server.h sock.h sockio.h	\
 thread.h  daemonize.h scheduler_client.h sys.h thread_pool.h misc.h	\
 hashlittle.h path.h common.h

SYS_SRC=daemonize.c  misc.c process_apple.c process.c process_linux.c	\
 thread.c hashlittle.c mem.c path.c process_bsd.c process_cygwin.c		\
 scheduler_client.c sockio.c thread_pool.c scheduler_server.c sock.c


lib_LTLIBRARIES=libaosys.la
libaosys_la_SOURCES=$(SYS_SRC)
libaosys_la_LIBADD=$(LDLIB)
libaosys_la_LDFLAGS = -avoid-version -static

if USE_SCHEDULER
#we compile a scheduler, which takes less memory than maos.
bin_PROGRAMS=../bin/scheduler 
___bin_scheduler_SOURCES=scheduler_server.c
___bin_scheduler_CFLAGS=$(AM_CFLAGS) -DSCHEDULER
___bin_scheduler_LDADD=libaosys.la $(LDOPTS)
___bin_scheduler_LDFLAGS=$(LDEXE)
endif

AM_CFLAGS=$(OPTS) -I $(top_srcdir)
LDADD= libsys.la -lm $(LDSYS) 
#It seems if bin_LDADD is specified. LDADD will not be used.

nobase_include_HEADERS=$(ICONS) scheduler_server.h scheduler_client.h		\
sockio.h proc.h process.h misc.h common.h mem.h thread.h	\
io.h sys.h shm.h private.h
SYS_SRC=proc.c proc_linux.c proc_apple.c proc_bsd.c proc_win.c sockio.c \
scheduler_server.c scheduler_client.c process.c misc.c mem.c thread.c io.c shm.c

lib_LTLIBRARIES=libsys.la
libsys_la_SOURCES=$(SYS_SRC)

if DO_THRPOOL #thread pool is only possible on new systems.
nobase_include_HEADERS+=thr_pool.h
SYS_SRC+=thr_pool.c
endif


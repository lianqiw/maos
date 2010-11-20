AM_CFLAGS=$(OPTS) -I $(top_srcdir)
LDADD= libsys.la -lm $(LDSYS) 
#It seems if bin_LDADD is specified. LDADD will not be used.
ICONS=icon-optics.h icon-cancel.h icon-running.h icon-draw.h icon-finished.h	\
icon-failed.h
nobase_include_HEADERS=$(ICONS) scheduler_server.h scheduler_client.h drawdaemon.h	\
sockio.h monitor.h proc.h process.h misc.h common.h mem.h thread.h io.h sys.h	\
shm.h private.h
SYS_SRC=proc.c proc_linux.c proc_apple.c proc_bsd.c proc_win.c sockio.c scheduler_server.c	\
scheduler_client.c process.c misc.c mem.c thread.c io.c shm.c 

lib_LTLIBRARIES=libsys.la
libsys_la_SOURCES=$(SYS_SRC)

if DO_THRPOOL #thread pool is only possible on new systems.
nobase_include_HEADERS+=thr_pool.h
SYS_SRC+=thr_pool.c
endif
PROGS=
if LINUX
PROGS+=../../bin/load ../../bin/record_cpu
______bin_load_SOURCES=load.c
______bin_load_CFLAGS=$(OPTScommon) $(AM_CFLAGS)
______bin_record_cpu_SOURCES=record_cpu.c 
______bin_record_CPU_CFLAGS=$(OPTScommon) $(AM_CFLAGS)
endif
if DRAW
PROGS+=drawdaemon
drawdaemon_SOURCES=drawdaemon.c 
drawdaemon_CFLAGS=$(DRAW_CFLAGS) $(AM_CFLAGS) $(OPTScommon)
drawdaemon_LDFLAGS=
drawdaemon_LDADD=$(LDADD) $(DRAW_LIBS) -lgthread-2.0
endif
if MON
if MONN
MONITOR_FLAGS=$(MONN_CFLAGS) $(AM_CFLAGS) $(OPTScommon) -DWITH_NOTIFY=1
else
MONITOR_FLAGS=$(MON_CFLAGS) $(AM_CFLAGS) $(OPTScommon) -DWITH_NOFITY=0
endif
PROGS+=../../bin/monitor ../../bin/monitor_list ../../bin/monitor_table
______bin_monitor_SOURCES=monitor.c monitor_gui.c
______bin_monitor_CFLAGS=$(MONITOR_FLAGS)
______bin_monitor_LDFLAGS=
______bin_monitor_LDADD=$(LDADD) $(MON_LIBS) 
______bin_monitor_list_SOURCES=monitor_gui_list.c monitor.c
______bin_monitor_list_CFLAGS=$(MONITOR_FLAGS)
______bin_monitor_list_LDFLAGS=
______bin_monitor_list_LDADD=$(LDADD) $(MON_LIBS)
______bin_monitor_table_SOURCES=monitor_gui_table.c monitor.c
______bin_monitor_table_CFLAGS=$(MONITOR_FLAGS)
______bin_monitor_table_LDFLAGS=
______bin_monitor_table_LDADD=$(LDADD) $(MON_LIBS)
endif

bin_PROGRAMS=$(PROGS)

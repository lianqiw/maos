AM_CFLAGS=$(OPTS) -I $(top_srcdir)
AM_LDFLAGS=$(LDOPTS) $(LDEXE)
ICONS=icon-clear.png icon-error.png icon-monitor.png icon-save.png  icon-waiting.png \
icon-connect.png icon-finished.png icon-play.png  icon-skip.png

noinst_HEADERS=monitor.h drawdaemon.h scheduler_ws.h mygtk.h monitor.html gtk3-css.h $(ICONS) 

PROGS=

if USE_SCHEDULER
PROGS+=../bin/scheduler
___bin_scheduler_SOURCES=scheduler.c
___bin_scheduler_LDADD=../sys/libaosys.la #$(LDOPTS)
if HAS_LWS
___bin_scheduler_SOURCES+=scheduler_ws.c
___bin_scheduler_LDADD+=-lwebsockets
endif
endif

if LINUX
PROGS+=../bin/load ../bin/record_cpu
___bin_load_SOURCES=load.c
___bin_record_cpu_SOURCES=record_cpu.c
endif

if GUI
if CAIRO
PROGS+=../bin/drawdaemon 
___bin_drawdaemon_SOURCES=drawdaemon.c drawdaemon_gui.c drawdaemon_io.c drawdaemon_draw.c resources.c
___bin_drawdaemon_CFLAGS=$(CAIRO_CFLAGS) $(GTK_CFLAGS) $(MAC_CFLAGS) $(AM_CFLAGS)
___bin_drawdaemon_LDADD=../sys/libaosys.la $(CAIRO_LIBS) $(GTK_LIBS) $(MAC_LIBS) -lgthread-2.0
endif
if USE_SCHEDULER
PROGS+=../bin/monitor 
___bin_monitor_SOURCES=monitor_gui.c monitor.c monitor_thread.c resources.c
___bin_monitor_CFLAGS=$(GTK_CFLAGS) $(MAC_CFLAGS) $(AM_CFLAGS)
___bin_monitor_LDADD=../sys/libaosys.la  $(GTK_LIBS) $(MAC_LIBS) #$(LDOPTS)
if NOTIFY
___bin_monitor_CFLAGS+=$(NOTIFY_CFLAGS) -DWITH_NOTIFY=1
___bin_monitor_LDADD+=$(NOTIFY_LIBS)
endif
endif
endif

LDADD=../lib/libaos.la $(LDEXE_FULL)
PROGS+=../bin/bin2fits ../bin/psfiris ../bin/psfgc ../bin/ethtest ../bin/mvm_cpu ../bin/drawbin ../bin/drawres 
___bin_bin2fits_SOURCES=bin2fits.c
___bin_psfiris_SOURCES=psfiris.c
___bin_psfgc_SOURCES=psfgc.c
___bin_ethtest_SOURCES=ethtest.c
___bin_mvm_cpu_SOURCES=mvm_cpu.c
___bin_drawbin_SOURCES=drawbin.c
___bin_drawres_SOURCES=drawres.c
bin_PROGRAMS=$(PROGS)

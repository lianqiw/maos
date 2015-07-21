AM_CFLAGS=$(OPTS) -I $(top_srcdir)
AM_LDFLAGS=$(LDOPTS) $(LDEXE)
LDADD=../lib/libaos.la

ICONS=icon-cancel.h icon-failed.h icon-monitor.h icon-running.h icon-draw.h	\
icon-finished.h icon-optics.h icon-waiting.h mouse_black.h mouse_hand.h mouse_white.h mygtk.h\
button-err.png button-ok.png button-play.png button-wait.png icon-waiting.png monitor.png mouse_black.png mouse_hand.png mouse_white.png monitor.html
noinst_HEADERS=monitor.h drawdaemon.h scheduler_ws.h $(ICONS) 

PROGS=

if USE_SCHEDULER
PROGS+=../bin/scheduler
___bin_scheduler_SOURCES=scheduler.c
___bin_scheduler_LDADD=../sys/libaosys.la $(LDOPTS)
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

if DRAW
PROGS+=../bin/drawdaemon ../bin/drawbin ../bin/drawres ../bin/calcenc
___bin_drawdaemon_SOURCES=drawdaemon.c drawdaemon_gui.c drawdaemon_io.c drawdaemon_draw.c
___bin_drawdaemon_CFLAGS=$(DRAW_CFLAGS) $(AM_CFLAGS)
___bin_drawdaemon_LDADD=$(LDADD) $(DRAW_LIBS) -lgthread-2.0
___bin_drawbin_SOURCES=drawbin.c
___bin_drawres_SOURCES=drawres.c
___bin_calcenc_SOURCES=calcenc.c
endif

if USE_SCHEDULER
if MON
if MONN #FULL with notify
MONITOR_FLAGS=$(MONN_CFLAGS) $(AM_CFLAGS) -DWITH_NOTIFY
MONITOR_LIBS=../sys/libaosys.la $(LDOPTS) $(MONN_LIBS) 
else #Monitor wihtout libnotify.
MONITOR_FLAGS=$(MON_CFLAGS) $(AM_CFLAGS)
MONITOR_LIBS=../sys/libaosys.la $(LDOPTS) $(MON_LIBS) 
endif
PROGS+=../bin/monitor ../bin/monitor_list ../bin/monitor_table
___bin_monitor_SOURCES=monitor.c monitor_gui.c monitor_thread.c
___bin_monitor_CFLAGS=$(MONITOR_FLAGS)
___bin_monitor_LDADD=$(MONITOR_LIBS)
___bin_monitor_list_SOURCES=monitor_gui_list.c monitor.c monitor_thread.c
___bin_monitor_list_CFLAGS=$(MONITOR_FLAGS)
___bin_monitor_list_LDADD=$(MONITOR_LIBS)
___bin_monitor_table_SOURCES=monitor_gui_table.c monitor.c monitor_thread.c
___bin_monitor_table_CFLAGS=$(MONITOR_FLAGS)
___bin_monitor_table_LDADD=$(MONITOR_LIBS)
endif
endif

PROGS+=../bin/servo_optim ../bin/slaving ../bin/cov2psf ../bin/bin2fits ../bin/psfiris ../bin/psfgc ../bin/mk2dcov ../bin/ethtest ../bin/mvm_cpu
___bin_servo_optim_SOURCES=servo_optim.c
___bin_slaving_SOURCES=slaving.c
___bin_cov2psf_SOURCES=cov2psf.c
___bin_bin2fits_SOURCES=bin2fits.c
___bin_psfiris_SOURCES=psfiris.c
___bin_psfgc_SOURCES=psfgc.c
___bin_mk2dcov_SOURCES=mk2dcov.c
___bin_ethtest_SOURCES=ethtest.c
___bin_mvm_cpu_SOURCES=mvm_cpu.c
bin_PROGRAMS=$(PROGS)

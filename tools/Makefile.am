AM_CFLAGS=$(OPTS) -I $(top_srcdir) -I $(top_srcdir)/lib/ -I $(top_srcdir)/lib/sys
AM_LDFLAGS=$(LDEXE)
if USE_LONG
LDADD=../lib/libaos_long.la $(LDOPTS)
else
LDADD=../lib/libaos_int.la $(LDOPTS)
endif
ICONS=icon-cancel.h icon-failed.h icon-monitor.h icon-running.h icon-draw.h	\
icon-finished.h icon-optics.h mouse_black.h  mouse_hand.h  mouse_white.h

noinst_HEADERS=monitor.h drawdaemon.h $(ICONS)
PROGS=
if LINUX
PROGS+=../bin/load ../bin/record_cpu
___bin_load_SOURCES=load.c
___bin_load_CFLAGS=$(OPTScommon) $(AM_CFLAGS)
___bin_record_cpu_SOURCES=record_cpu.c
___bin_record_CPU_CFLAGS=$(OPTScommon) $(AM_CFLAGS)
endif
if DRAW
PROGS+=../bin/drawdaemon ../bin/drawbin ../bin/drawres ../bin/calcenc
___bin_drawdaemon_SOURCES=drawdaemon.c drawdaemon_gui.c drawdaemon_io.c drawdaemon_draw.c
___bin_drawdaemon_CFLAGS=$(DRAW_CFLAGS) $(AM_CFLAGS) $(OPTScommon)
___bin_drawdaemon_LDADD=$(LDADD) $(DRAW_LIBS) -lgthread-2.0
___bin_drawbin_SOURCES=drawbin.c
___bin_drawres_SOURCES=drawres.c
___bin_calcenc_SOURCES=calcenc.c
endif
if MON
if MONN  #FULL with notify
MONITOR_FLAGS=$(MONN_CFLAGS) $(AM_CFLAGS) $(OPTScommon) -DWITH_NOTIFY
MONITOR_LIBS=$(LDADD) $(MONN_LIBS)
else     #Monitor wihtout libnotify.
MONITOR_FLAGS=$(MON_CFLAGS) $(AM_CFLAGS) $(OPTScommon)
MONITOR_LIBS=$(LDADD) $(MON_LIBS)
endif
PROGS+=../bin/monitor ../bin/monitor_list ../bin/monitor_table
___bin_monitor_SOURCES=monitor.c monitor_gui.c
___bin_monitor_CFLAGS=$(MONITOR_FLAGS)
___bin_monitor_LDADD=$(MONITOR_LIBS)
___bin_monitor_list_SOURCES=monitor_gui_list.c monitor.c #gtkcellrendererprogressnew.c gtkcellrendererprogressnew.h
___bin_monitor_list_CFLAGS=$(MONITOR_FLAGS)
___bin_monitor_list_LDADD=$(MONITOR_LIBS)
___bin_monitor_table_SOURCES=monitor_gui_table.c monitor.c
___bin_monitor_table_CFLAGS=$(MONITOR_FLAGS)
___bin_monitor_table_LDADD=$(MONITOR_LIBS)
endif
PROGS+=../bin/servo_optim ../bin/slaving ../bin/cov2psf ../bin/bin2fits
bin_PROGRAMS=$(PROGS)
___bin_servo_optim_SOURCES=servo_optim.c
___bin_slaving_SOURCES=slaving.c
___bin_cov2psf_SOURCES=cov2psf.c
___bin_bin2fits_SOURCES=bin2fits.c
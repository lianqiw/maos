AM_CFLAGS=$(OPTS) -I $(top_srcdir)
AM_LDFLAGS=$(LDEXE)
if USE_LONG
LDADD=../lib/libaos_long.la ../math/libaomath_long.la ../sys/libaosys.la $(LDOPTS)
else
LDADD=../lib/libaos_int.la ../math/libaomath_int.la ../sys/libaosys.la $(LDOPTS)
endif

ICONS=icon-cancel.h icon-failed.h icon-monitor.h icon-running.h icon-draw.h	\
icon-finished.h icon-optics.h icon-waiting.h mouse_black.h  mouse_hand.h  mouse_white.h mygtk.h

noinst_HEADERS=monitor.h drawdaemon.h $(ICONS)

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
___bin_load_CFLAGS=$(OPTScommon) $(AM_CFLAGS)
___bin_record_cpu_SOURCES=record_cpu.c
___bin_record_CPU_CFLAGS=$(OPTScommon) $(AM_CFLAGS)
endif

if DRAW
PROGS+=../bin/drawdaemon ../bin/drawbin ../bin/drawres ../bin/calcenc
___bin_drawdaemon_SOURCES=drawdaemon.c drawdaemon_gui.c drawdaemon_io.c drawdaemon_draw.c
___bin_drawdaemon_CFLAGS=$(DRAW_CFLAGS) $(AM_CFLAGS) $(OPTScommon)
___bin_drawdaemon_LDADD=$(LDADD) $(DRAW_LIBS)  -lgthread-2.0
___bin_drawbin_SOURCES=drawbin.c
___bin_drawres_SOURCES=drawres.c
___bin_calcenc_SOURCES=calcenc.c
endif

if USE_SCHEDULER
if MON
if MONN  #FULL with notify
MONITOR_FLAGS=$(MONN_CFLAGS)  $(AM_CFLAGS) $(OPTScommon) -DWITH_NOTIFY
MONITOR_LIBS=../sys/libaosys.la $(LDOPTS) $(MONN_LIBS) 
else     #Monitor wihtout libnotify.
MONITOR_FLAGS=$(MON_CFLAGS)  $(AM_CFLAGS) $(OPTScommon)
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

if USE_CUDA
PROGS+=../bin/gpu_mvm
___bin_gpu_mvm_SOURCES=gpu_mvm.c
___bin_gpu_mvm_LDADD=../cuda/libgpu.a $(LDADD) $(LDCUDA)
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

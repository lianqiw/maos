AM_CFLAGS=-I $(top_srcdir)

ICONS=icon-avg.png icon-cancel.png icon-clear.png icon-connect.png icon-draw-small.png icon-draw.png icon-error.png \
icon-finished.png icon-log.png icon-monitor.png icon-play.png icon-save.png icon-skip.png icon-waiting.png favicon.ico

noinst_HEADERS=monitor.h drawdaemon.h scheduler.h mygtk.h monitor.html gtk3-css.h
EXTRA_DIST=resources.xml monitor.html monitor.js monitor.css monitor_utils.js monitor_draw_utils.js monitor_draw.js $(ICONS)
PROGS=

if USE_SCHEDULER
PROGS+=../bin/scheduler
___bin_scheduler_SOURCES=scheduler.c scheduler_http.c scheduler_proxy.c
___bin_scheduler_LDADD=../sys/libaosys.la $(OPENSSL_LIBS) $(ZSTD_LIBS) $(EXE_LIBS)
___bin_scheduler_LDFLAGS=$(OPENSSL_LDFLAGS) $(AM_LDFLAGS)
___bin_scheduler_CFLAGS=$(OPENSSL_INCLUDES) $(ZSTD_CFLAGS)
if HAS_LWS
___bin_scheduler_SOURCES+=scheduler_lws.c
___bin_scheduler_LDADD+=$(LWS_LIBS)
___bin_scheduler_CFLAGS+=$(LWS_CFLAGS)
endif
endif

if LINUX
PROGS+=../bin/load ../bin/record_cpu
___bin_load_SOURCES=load.c
___bin_record_cpu_SOURCES=record_cpu.c
endif

if GUI
if USE_SCHEDULER
PROGS+=../bin/monitor
___bin_monitor_SOURCES=monitor_gui.c monitor.c monitor_thread.c resources.c mygtk.c
___bin_monitor_CFLAGS=$(GTK_CFLAGS) $(MAC_CFLAGS) $(AM_CFLAGS)
___bin_monitor_LDADD=../sys/libaosys.la $(GTK_LIBS) $(MAC_LIBS) $(EXE_LIBS)
if NOTIFY
___bin_monitor_CFLAGS+=$(NOTIFY_CFLAGS) -DWITH_NOTIFY=1
___bin_monitor_LDADD+=$(NOTIFY_LIBS)
endif
endif
if CAIRO
PROGS+=../bin/drawdaemon
___bin_drawdaemon_SOURCES=drawdaemon.c drawdaemon_gui.c drawdaemon_io.c drawdaemon_draw.c resources.c mygtk.c
___bin_drawdaemon_CFLAGS=$(CAIRO_CFLAGS) $(GTK_CFLAGS) $(MAC_CFLAGS) $(AM_CFLAGS)
___bin_drawdaemon_LDADD=../sys/libaosys.la $(CAIRO_LIBS) $(GTK_LIBS) $(MAC_LIBS) -lgthread-2.0 $(EXE_LIBS)
endif
endif

LDADD=../lib/libaos.la $(EXEFULL_LIBS)
PROGS+=../bin/bin2fits ../bin/psfiris ../bin/psfgc ../bin/ethtest ../bin/mvm_cpu ../bin/drawbin ../bin/drawres ../bin/readbin
___bin_bin2fits_SOURCES=bin2fits.c
___bin_psfiris_SOURCES=psfiris.c
___bin_psfgc_SOURCES=psfgc.c
___bin_ethtest_SOURCES=ethtest.c
___bin_mvm_cpu_SOURCES=mvm_cpu.c
___bin_drawbin_SOURCES=drawbin.c
___bin_drawres_SOURCES=drawres.c
___bin_readbin_SOURCES=readbin.c
bin_PROGRAMS=$(PROGS)

resources.c:$(ICONS) resources.xml
	@cd $(srcdir) && glib-compile-resources resources.xml --generate-source --target=resources.c
if DARWIN
bindir_app = $(abs_top_builddir)/bin
APPS = drawdaemon
if USE_SCHEDULER
APPS += monitor
endif
# per-app config (naming convention: <App>_*)
monitor_EXE = ../bin/monitor
monitor_ICNS = icon-monitor.icns
monitor_BUNDLE_ID = io.github.lianqiw.monitor
monitor_EXEC_NAME = monitor
monitor_ICON_NAME = icon-monitor
drawdaemon_EXE = ../bin/drawdaemon
drawdaemon_ICNS = icon-draw.icns
drawdaemon_BUNDLE_ID = io.github.lianqiw.drawdaemon
drawdaemon_EXEC_NAME = drawdaemon
drawdaemon_ICON_NAME = icon-drawdaemon

# build all bundles
all-local: $(APPS:%=$(bindir_app)/%.app)

define MAKE_APP_BUNDLE
$(bindir_app)/$1.app: $$($1_EXE) $$($1_ICNS)
	$$(MKDIR_P) $$@/Contents/MacOS
	$$(MKDIR_P) $$@/Contents/Resources

	printf '%s\n' \
'<?xml version="1.0" encoding="UTF-8"?>' \
'<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN"' \
' "http://www.apple.com/DTDs/PropertyList-1.0.dtd">' \
'<plist version="1.0">' \
'<dict>' \
'    <key>CFBundleName</key>' \
'    <string>$1</string>' \
'    <key>CFBundleExecutable</key>' \
'    <string>$$($1_EXEC_NAME)</string>' \
'    <key>CFBundleIdentifier</key>' \
'    <string>$$($1_BUNDLE_ID)</string>' \
'    <key>CFBundlePackageType</key>' \
'    <string>APPL</string>' \
'    <key>CFBundleVersion</key>' \
'    <string>1.0</string>' \
'    <key>CFBundleShortVersionString</key>' \
'    <string>1.0</string>' \
'    <key>CFBundleIconFile</key>' \
'    <string>$$($1_ICON_NAME)</string>' \
'</dict>' \
'</plist>' \
> $$@/Contents/Info.plist
	cp $(srcdir)/$$($1_ICNS) $$@/Contents/Resources/$$($1_ICON_NAME).icns
	cp $$($1_EXE) $$@/Contents/MacOS/$$($1_EXEC_NAME)
	touch $$@
endef
$(foreach app,$(APPS),$(eval $(call MAKE_APP_BUNDLE,$(app))))

clean-local:
	rm -rf $(APPS:%=$(bindir_app)/%.app)
endif #DARWIN
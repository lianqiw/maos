AC_DEFUN([MY_CHECK_GTK],[
	AC_ARG_WITH(gtk,  AS_HELP_STRING([--with-gtk[=DIR]],[Specify non-standard GTK+2 directory (obsolete).]), [GTKDIR=${withval}], [GTKDIR=])
	AC_ARG_ENABLE(gtk-3, AS_HELP_STRING([--disable-gtk-3],[Disable GTK+-3]), [use_gtk_3="$enableval"], [use_gtk_3=])
	AC_ARG_ENABLE(gtk-4, AS_HELP_STRING([--disable-gtk-4],[Disable GTK4]), [use_gtk_4="$enableval"], [use_gtk_4=])

	gtk_ver=0
	has_cairo=0
	has_notify=0
	has_mac=0
	m4_ifdef([PKG_PROG_PKG_CONFIG], [
		#Check availability of gtk and libnotify. If available, we will build drawdaemon and the monitor.
		#cairo 1.2 comes with gtk+-2.10
		#make sure there are spaces before and after >=
		if test -n "$GTKDIR" -a -d "$GTKDIR" ;then #User specified GTK directory
			gtk_ver="2"
			if test -d "$GTKDIR/pkgconfig" ;then
				GTKDIR="$GTKDIR/../"
			fi
			if test -d "$GTKDIR/lib/pkgconfig" ;then
				export PKG_CONFIG_PATH="$PKG_CONFIG_PATH:$GTKDIR/lib/pkgconfig"
			fi
			TT=-I$GTKDIR/include
			GTK_CFLAGS="$TT $TT/cairo $TT/atk-1.0 $TT/gtk-${gtk_ver}.0 $TT/glib-2.0 $TT/pango-1.0 $TT/gdk-pixbuf-2.0 -I$GTKDIR/lib/glib-2.0/include -I${GTKDIR}/lib/gtk-2.0/include"
			GTK_LIBS="-L$GTKDIR/lib/ -Wl,-rpath,$GTKDIR/lib/"
			if test -d "$GTKDIR/bin" ;then
				GTK_LIBS="$GTK_LIBS -L$GTKDIR/bin -Wl,-rpath,GTKDIR/bin"
			fi
		else
			AS_IF([test "$gtk_ver" = "0" -a "$use_gtk_4" = "yes" ], [PKG_CHECK_MODULES([GTK], gtk4, [gtk_ver="4"], [gtk_ver="0"])])
			AS_IF([test "$gtk_ver" = "0" -a "$use_gtk_3" != "no" ], [PKG_CHECK_MODULES([GTK], gtk+-3.0, [gtk_ver="3"], [gtk_ver="0"])]) #prefer gtk+3
			AS_IF([test "$gtk_ver" = "0" -a "$use_gtk_4" != "no" ], [PKG_CHECK_MODULES([GTK], gtk4, [gtk_ver="4"], [gtk_ver="0"])])
			AS_IF([test "$gtk_ver" = "0" ],                         [PKG_CHECK_MODULES([GTK], gtk+-2.0 >= 2.2,    [gtk_ver="2"], [gtk_ver="0"])])
		fi
		AS_IF([test "$gtk_ver" != "0"], [
			PKG_CHECK_MODULES(DRAW,   cairo >= 1.4,    [has_cairo="1"], [has_cairo="0"])
			PKG_CHECK_MODULES(NOTIFY, libnotify >= 0.1,[has_notify="1"], [has_notify="0"])
			PKG_CHECK_MODULES(MAC,    gtk-mac-integration-gtk${gtk_ver}, [has_mac="1"], [has_mac="0"])
			])

		if test "$gtk_ver" = "2" ;then
			if ! echo "$GTK_LIBS" | grep -q gthread ;then
				#OLD GLIB has GTHREAD separated
				PKG_CHECK_MODULES(GTHREAD, gthread-2.0 >= 2.0, [has_gthread="yes"], [has_gthread="no"])
				if test -n "$GTHREAD_LIBS" ;then
				GTK_LIBS+=" $GTHREAD_LIBS"
				fi
			fi
			GTK_CFLAGS+=" -Wno-deprecated-declarations"
		elif test "$gtk_ver" = "3" ;then
			GTK_CFLAGS+=" -DGDK_DISABLE_DEPRECATED -DGTK_DISABLE_DEPRECATED"
		fi
		AC_DEFINE_UNQUOTED(MAC_INTEGRATION, [$has_mac], [Has GTK MAC integration])
		#The following will substitute the variables name
		#Reexport after updating
		AC_SUBST(GTK_CFLAGS)
		AC_SUBST(GTK_LIBS)
	])

	AM_CONDITIONAL(GUI,    [test x$gtk_ver  != x0])
	AM_CONDITIONAL(CAIRO,  [test x$has_cairo = x1])
	AM_CONDITIONAL(NOTIFY, [test x$has_notify = x1])
	if test "$gtk_ver" != "0" ;then
		echo "GTK_LIBS=$GTK_LIBS"
	fi
])

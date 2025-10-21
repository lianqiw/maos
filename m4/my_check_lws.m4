AC_DEFUN([MY_CHECK_LWS],[
	AC_ARG_ENABLE(libwebsockets, AS_HELP_STRING([--disable-libwebsockets],[Disable libwebsockets]))
	#Check for libwebsocket
	has_lws=no
	if test "$enable_libwebsockets" != "no" ;then
		save_CFLAGS="$CFLAGS"
		m4_ifdef([PKG_PROG_PKG_CONFIG], [PKG_PROG_PKG_CONFIG
			PKG_CHECK_MODULES([LWS], [libwebsockets], [has_lws="yes" CFLAGS="$CFLAGS $LWS_CFLAGS"],[has_lws="no" LWS_LIBS="-lwebsockets"])
		])
		if test "$has_lws" = "no" ;then #then download pre-compiled library and header.
			AC_CHECK_LIB([websockets], [lws_service], [has_lws="yes"], [has_lws="no"], [])
			if test "$has_lws" = "no" ;then #then download pre-compiled library and header.
				if test "$enable_binary" != "no";then
					MY_DOWNLOAD([LIBWEBSOCKETS], [libwebsockets${libsuffix}.tar.bz2], [${prefix}])
					unset ac_cv_lib_websockets_lws_service
					AC_CHECK_LIB([websockets], [lws_service], [has_lws="yes"], [has_lws="no"], [])
				fi
				if test "$has_lws" = "no" ;then #try to compile from source. Needs cmake
					MY_COMPILE([LIBWEBSOCKETS], [libwebsockets.tar.bz2], [-DLWS_WITH_SSL=OFF -DLWS_IPV6=OFF -DLWS_WITHOUT_TESTAPPS=ON -DLWS_WITHOUT_DAEMONIZE=ON -DLWS_WITHOUT_CLIENT=ON -DLWS_WITH_SHARED=ON -DLWS_WITH_STATIC=ON -DLWS_WITH_EXTERNAL_POLL=ON])
					unset ac_cv_lib_websockets_lws_service
					AC_CHECK_LIB([websockets], [lws_service], [has_lws="yes"], [has_lws="no"], [])
				fi
			fi
		fi
		CFLAGS="$CFLAGS $OPENSSL_INCLUDES" #system libwebsockets needs OPENSSL
		AC_CHECK_HEADERS([libwebsockets.h], [], [has_lws="no"])
		if test "$has_lws" = "yes" ;then
			AC_SUBST(LWS_LIBS)
			AC_SUBST(LWS_CFLAGS)
			$1
		else
			$2
		fi
		CFLAGS=${save_CFLAGS}
		echo LWS_LIBS=$LWS_LIBS
		echo LWS_CFLAGS=$LWS_CFLAGS
	fi
	AM_CONDITIONAL(HAS_LWS, [test "$has_lws" = "yes"])
])

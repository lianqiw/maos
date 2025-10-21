AC_DEFUN([MY_CHECK_CMOCKA],[
	if test "$USE_CPP" = "1" ;then
		cmocka_default="no" #the C++ compiler does not work with cmocka
	fi
	AC_ARG_ENABLE(cmocka, AS_HELP_STRING([--enable-cmocka], [Enable test with cmocka]), [enable_cmocka="$enableval"], [enable_cmocka="$cmocka_default"])
	has_cmocka="no"
	if test "$enable_cmocka" != "no" ;then
		save_CFLAGS="$CFLAGS"
		m4_ifdef([PKG_PROG_PKG_CONFIG], [PKG_PROG_PKG_CONFIG
			PKG_CHECK_MODULES([CMOCKA], [cmocka], [has_cmocka="yes" CFLAGS="$CFLAGS $CMOCKA_CFLAGS"],[has_cmocka="no" CMOCKA_LIBS="-lcmocka"])
		])
		if test $has_cmocka = "no";then
			AC_CHECK_LIB([cmocka], [mock_assert], [has_cmocka="yes"], [has_cmocka="no"])
			if test "$has_cmocka" = "no" ;then
				if test "$enable_binary" != "no";then
					MY_DOWNLOAD([CMOCKA], [cmocka${libsuffix}.tar.bz2], [${prefix}])
					unset ac_cv_lib_cmocka_mock_assert
					AC_CHECK_LIB([cmocka], [mock_assert], [has_cmocka="yes"], [has_cmocka="no"])
				fi
				if test "$has_cmocka" = "no" ;then
					MY_COMPILE([CMOCKA], [cmocka.tar.xz])
					unset ac_cv_lib_cmocka_mock_assert
					AC_CHECK_LIB([cmocka], [mock_assert], [has_cmocka="yes"], [has_cmocka="no"])
				fi
			fi
		fi
		AC_CHECK_HEADERS([cmocka.h],[:],[has_cmocka="no"],[
			#include <stdarg.h>
			#include <stddef.h>
			#include <stdint.h>
			#include <setjmp.h>])
	
		if test "$has_cmocka" = "yes" ;then
			AC_SUBST(CMOCKA_LIBS)
			$1
		else :
			$2
		fi
		CFLAGS=${save_CFLAGS}
		echo CMOCKA_LIBS=$CMOCKA_LIBS
		echo CMOCKA_CFLAGS=$CMOCKA_CFLAGS
	fi
	AM_CONDITIONAL(HAS_CMOCKA, [ test x$has_cmocka = xyes ])
])

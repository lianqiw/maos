AC_DEFUN([MY_CHECK_ZSTD],[
	AC_ARG_ENABLE(zstd, AS_HELP_STRING([--disable-zstd], [Disable zstd file compression]))
	have_zstd=no
	if test "$enable_zstd" != no ;then
		m4_ifdef([PKG_PROG_PKG_CONFIG], [
			PKG_CHECK_MODULES([ZSTD], [libzstd >= 1.5],	[have_zstd=yes],[:])
		])
		if test "$have_zstd" = no ;then
			AC_CHECK_LIB([zstd], [ZSTD_compress], [ZSTD_LIBS="-lzstd"
				AC_CHECK_HEADER([zstd.h], [have_zstd=yes ZSTD_CFLAGS=""])
			])
		fi
	fi
	if test "$have_zstd" = yes; then
		AC_DEFINE([HAVE_LIBZSTD], [1], [Define if libzstd is available])
		AC_SUBST(ZSTD_CFLAGS)
		AC_SUBST(ZSTD_LIBS)
		echo ZSTD_LIBS=$ZSTD_LIBS
		echo ZSTD_CFLAGS=$ZSTD_CFLAGS
		AC_MSG_NOTICE([libzstd is found])
		$1
	else
		AC_MSG_NOTICE([libzstd is not found])
		$2
	fi
])

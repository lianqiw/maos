AC_DEFUN([MY_CHECK_CHOLMOD],[
	save_CFLAGS="$CFLAGS"
	has_chol="no"
	m4_ifdef([PKG_PROG_PKG_CONFIG], [
		PKG_CHECK_MODULES([CHOL], [cholmod], [CFLAGS="$CFLAGS $CHOL_CFLAGS" has_chol="yes"],[:])
	])
	AC_ARG_ENABLE(dlong, AS_HELP_STRING([--enable-dlong], [Use sparse matrix with long index]))
	#Check for cholmod
	if test "$enable_binary" != "no";then
		steps="1 2 3"
	else
		steps="1 3"
	fi
	for step in $steps ;do
		case $step in
		1)
			;;
		2)
			MY_DOWNLOAD([CHOLMOD], [cholmod${libsuffix}.tar.bz2], [${prefix}])
			;;
		3)
			if test "$enable_dlong" != "no"; then
				CPPFLAGS="$CPPFLAGS -DDLONG"
			fi
			MY_COMPILE([CHOLMOD], [cholmod.tar.bz2])
			;;
		esac
		unset ac_cv_lib_cholmod_cholmod_factorize
		unset ac_cv_lib_cholmod_cholmod_l_factorize
		if test "$enable_dlong" != "no"; then
			AC_CHECK_LIB([cholmod], [cholmod_l_factorize], [has_chol="yes" CPPFLAGS="$CPPFLAGS -DDLONG"], [has_chol="no"], [$CHOL_LIBS $LAPACK])
		fi
		if test "$has_chol" = "no" ;then : #do the test at once to set CPPFLAGS
			AC_CHECK_LIB([cholmod], [cholmod_factorize], [has_chol="yes"], [has_chol="no"], [$CHOL_LIBS $LAPACK])
		fi
		if test "$has_chol" = "yes" ;then : #do the test at once to set CPPFLAGS
			break;
		fi
	done
	unset ac_cv_header_cholmod_h ac_cv_header_suitesparse_cholmod_h
	AC_CHECK_HEADERS([cholmod.h], [], [AC_CHECK_HEADERS([suitesparse/cholmod.h], [], [has_chol="no"])])
	if test -z "$CHOL_LIBS" ;then
		CHOL_LIBS="-lcholmod"
	fi
	if test "$has_chol" = "yes" ;then
		AC_SUBST(CHOL_CFLAGS)
		AC_SUBST(CHOL_LIBS)
		$1
	else :
		$2
	fi
	echo CHOL_CFLAGS=$CHOL_CFLAGS
	echo CHOL_LIBS=$CHOL_LIBS
	CFLAGS=${save_CFLAGS}
])

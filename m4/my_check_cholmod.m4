AC_DEFUN([MY_CHECK_CHOLMOD],[
	save_CFLAGS="$CFLAGS"
	m4_ifdef([PKG_PROG_PKG_CONFIG], [PKG_PROG_PKG_CONFIG
		PKG_CHECK_MODULES([CHOL], [cholmod], [CFLAGS="$CFLAGS $CHOL_CFLAGS"],[CHOL_LIBS="-lcholmod"])
	])
	
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
			MY_COMPILE([CHOLMOD], [cholmod.tar.bz2])
			;;
		esac
		unset ac_cv_lib_m_cholmod_factorize
		unset ac_cv_lib_m_cholmod_l_factorize
		AC_CHECK_LIB([m], [cholmod_l_factorize], [has_chol="yes" CPPFLAGS="$CPPFLAGS -DDLONG" break], [
			AC_CHECK_LIB([m], [cholmod_factorize], [has_chol="yes" beak], [has_chol="no"], [$CHOL_LIBS $LAPACK])
		], [$CHOL_LIBS $LAPACK])
	done
	unset ac_cv_header_cholmod_h ac_cv_header_suitesparse_cholmod_h
	AC_CHECK_HEADERS([cholmod.h], [], [AC_CHECK_HEADERS([suitesparse/cholmod.h], [], [has_chol="no"])])
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

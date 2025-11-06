AC_DEFUN([MY_CHECK_LAPACK], [
  AC_ARG_ENABLE([mkl],
    AS_HELP_STRING([--enable-mkl], [Prefer MKL]),
    [use_mkl="$enableval"],
    [use_mkl="default"])

	default="default pkg-config default-3"
	if test "x$enable_threadlibx" != xno ; then
   		openblas="openblas-threads openblas-threads-0 openblas-serial openblas-serial-0"
  	else
   		openblas="openblas-serial openblas-serial-0"
  	fi
	try_lapack="$default $openblas"
	# Determine MKL availability based on CPU
	case "$host_os" in
		*linux*)
			if grep -qE "GenuineIntel|AuthenticAMD" /proc/cpuinfo; then
				#"MKL is good for Intel and AMD CPUs"
				if test x$use_mkl = xyes ;then
					try_lapack="mkl mkl-download $try_lapack"
				else
					try_lapack="mkl $try_lapack mkl-download"
				fi
			fi
		;;
	esac

	# Loop over LAPACK candidates
	has_lapack=no
	has_blas=no
	dgemm=dgemm_
	dposv=dposv_

	for try in $try_lapack; do
		echo trying ${try}
		case "${try}" in
			mkl*)
				if test "$try" = "mkl-download"; then
					MY_DOWNLOAD([MKL], [mkl${libsuffix}.tar.bz2], [${prefix}], [1])
				fi

				if test "$host_cpu" = "x86_64" ; then
					LAPACK="-lmkl_intel_lp64"
				else
					LAPACK="-lmkl_intel"
				fi

				if test -z "$OPENMP_CFLAGS" -o "x$enable_threadlib" = "xno" ; then
					LAPACK="$LAPACK -lmkl_sequential"
				else
					if test "$CN" = "ICC" -o "$CN" = "CLANG" ; then
						LAPACK="$LAPACK -lmkl_intel_thread"
					else
						LAPACK="$LAPACK -lmkl_gnu_thread"
					fi
				fi
				LAPACK="$LAPACK -lmkl_core"

				if test "$try" = "mkl-static" ; then
					lapack_mkl_static="yes"
					LAPACK="-Wl,--start-group,$LAPACK,--end-group"
					LAPACK=$(echo "$LAPACK" | tr ' ' ',') #Convert spaces to commas (POSIX-safe)
				fi

				if test "x$enable_threadlib" != "xno" ; then
					LAPACK="$LAPACK $OPENMP_CFLAGS -lpthread -ldl"
				fi
				;;

			blis*)
				if test "$try" = "blis-download"; then
					MY_DOWNLOAD([blis], [blis${libsuffix}.tar.bz2], [${prefix}], [1])
				fi

				if test -z "$OPENMP_CFLAGS" ; then
					LAPACK="-lblis"
				else
					LAPACK="-lblis-mt"
				fi
				LAPACK="$LAPACK -llapack"
				;;

			pkg-config)
				LAPACK=$(pkg-config --libs blas lapack 2>/dev/null)
				;;

			openblas-threads)
				if test -z "$OPENMP_CFLAGS" ; then
					LAPACK="-lopenblasp"
				else
					LAPACK="-lopenblaso"
				fi
				;;

			openblas-threads-0)
				if test -z "$OPENMP_CFLAGS" ; then
					LAPACK="-l:libopenblasp.so.0"
				else
					LAPACK="-l:libopenblaso.so.0"
				fi
				;;

			openblas-serial)
				LAPACK="-lopenblas"
				;;

			openblas-serial-0)
				LAPACK="-l:libopenblas.so.0"
				;;

			default)
				LAPACK="-llapack -lblas"
				;;

			default-3)
				LAPACK="-l:liblapack.so.3 -l:libblas.so.3"
				;;

			*)
				LAPACK=""
				;;
		esac
		if test -n "$LAPACK" ; then
			LAPACK="$LAPACK -lm"
			unset ac_cv_lib_m_${dgemm}
			unset ac_cv_lib_m_${dposv}
			AC_CHECK_LIB([m], [${dposv}], [has_lapack=yes], [has_lapack=no], [$LAPACK])
			AC_CHECK_LIB([m], [${dgemm}], [has_blas=yes], [has_blas=no], [$LAPACK])
			if test "$has_lapack" = "yes" -a "$has_blas" = "yes" ; then
				break;
			fi
		fi
	done
	if test "$has_lapack" = "yes" -a "$has_blas" = "yes" ; then
		$1
		AC_SUBST(LAPACK)
		break;
	else :
		$2
	fi
	echo LAPACK=$LAPACK
])

AC_DEFUN([MY_CHECK_FFT],[
	save_CFLAGS="${CFLAGS}"
	has_fftw="no"
	m4_ifdef([PKG_PROG_PKG_CONFIG], [
		PKG_CHECK_MODULES([FFTW], [fftw3], [has_fftw=yes CFLAGS="$CFLAGS $FFTW_CFLAGS"],[:])
	])
	#Check FFTW header and library. 
	if test "$enable_binary" != "no";then
		steps="1 2 3"
	else
		steps="1 3"
	fi
	for step in $steps; do
		case $step in
		1)
		;;
		2)
			MY_DOWNLOAD([FFTW], [fftw${libsuffix}.tar.bz2], [${prefix}])
		;;
		3) #compile from source
			FLAGS="--enable-threads --disable-fortran "
			if test "$host_cpu" = "aarch64" ;then
				FLAGS="$FLAGS --enable-neon"
				if test "$system" = "apple" ;then
					FLAGS="$FLAGS --build=aarch64-apple-darwin"
				else
					FLAGS="$FLAGS --build=aarch64 "
				fi
			else
				FLAGS="$FLAGS --enable-sse2 --enable-avx --enable-avx2 --enable-avx512"
			fi
			MY_COMPILE([FFTW], [fftw.tar.bz2], [$FLAGS])
			MY_COMPILE([FFTW], [fftw.tar.bz2], [$FLAGS --enable-float])
		;;
		esac
		unset ac_cv_lib_fftw3_fftw_execute
		unset ac_cv_header_fftw3_h
		AC_CHECK_LIB([fftw3], [fftw_execute], [has_fftw="yes"], [has_fftw="no"],[$FFTW_LIBS])
		AC_CHECK_HEADERS([fftw3.h], [:], [has_fftw="no"])
		if test "$has_fftw" = "yes" ;then :
			break;
		fi
	done
	if test -z "$FFTW_LIBS" ;then
		FFTW_LIBS="-lfftw3"
	fi
	if test "$has_fftw" = "yes" ;then :
		$1
	else :
		$2
	fi
	AC_CHECK_LIB([fftw3f], [fftwf_execute],[has_fftwf=1 FFTW_LIBS="$FFTW_LIBS -lfftw3f"] ,[has_fftwf=0],[$FFTW_LIBS])

	if test "$has_fftw" = "yes" -a x$enable_threadlib != xno;then
		#if test "$has_omp" = "1" ;then
		#	THREADS_SUFFIX="_omp"
		#else
		THREADS_SUFFIX="_threads"
		#fi
		FFTW_LIBS_THREADS="-pthread"
		for name in $FFTW_LIBS ;do
			case "$name" in
				*.la)
				FFTW_LIBS_THREADS="$FFTW_LIBS_THREADS ${name/.la/}${THREADS_SUFFIX}.la"
				;;
				*.a)
				FFTW_LIBS_THREADS="$FFTW_LIBS_THREADS ${name/.la/}${THREADS_SUFFIX}.a"
				;;
				*.so)
				FFTW_LIBS_THREADS="$FFTW_LIBS_THREADS ${name/.la/}${THREADS_SUFFIX}.so"
				;;
				-lfft*)
				FFTW_LIBS_THREADS="$FFTW_LIBS_THREADS ${name}${THREADS_SUFFIX}"
				;;
			esac
		done
		AC_CHECK_LIB([fftw3], [fftw_plan_with_nthreads], [fftw_threads=1 FFTW_LIBS="$FFTW_LIBS $FFTW_LIBS_THREADS"], [fftw_threads=0], [$FFTW_LIBS_THREADS $FFTW_LIBS ])
	else
		fftw_threads=0
	fi
	
	AC_DEFINE_UNQUOTED(HAS_FFTW_THREADS, [$fftw_threads], [FFTW has thread support.])
	AC_DEFINE_UNQUOTED(HAS_FFTWF, [$has_fftwf], [FFTW has single precision support.])
	AC_CHECK_LIB([fftw3], [fftw_threads_set_callback], [has_fftw_callback=1], [has_fftw_callback=0], [$FFTW_LIBS]) #since version 3.3.9
	AC_DEFINE_UNQUOTED(HAS_FFTW_CALLBACK, [$has_fftw_callback], [FFT has fftw_threads_set_callback support])
	AS_IF([test x$has_fftwf = 0 -a x$use_double = xno], [AC_MSG_ERROR([FFT single precision is not available for --disable-double])])
	AC_SUBST(FFTW_CFLAGS)
	AC_SUBST(FFTW_LIBS)
	echo FFTW_CFLAGS=$FFTW_CFLAGS
	echo FFTW_LIBS=$FFTW_LIBS
	CFLAGS=${save_CFLAGS}
])

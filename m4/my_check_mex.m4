AC_DEFUN([MY_CHECK_MEX],[
	#Compile mex routines if MATLAB is detected and single is not enabled
	MEXDIR_AUTO=
	if which matlab >/dev/null 2>&1;then
		MEX="$(which matlab)"
		if test -L "$MEX" ;then
			MEX="$(readlink "$MEX")"
		fi
		MEXDIR_AUTO="$(dirname "$MEX")"
		echo MEXDIR_AUTO=$MEXDIR_AUTO
		unset MEX
	else
		for i in /Applications/MATLAB*.app $HOME/Applications/MATLAB*.app /opt/MATLAB/* /usr/local/share/MATLAB/* ;do
			if test -x "$i/bin/mex" ;then
				MEXDIR_AUTO="${i}/bin"
				echo "Matlab found at ${i}"
				break;
			fi
		done
	fi
	AC_ARG_WITH(matlab, AS_HELP_STRING([--with-matlab[=DIR]],[Enable compiling mex routines for matlab]), [MEXDIR=${withval}], [MEXDIR=${MEXDIR_AUTO}])
	if test -z "$MEXDIR" -o "$MEXDIR" = "no" ;then
		with_matlab=no
	else
		with_matlab=yes
	fi
	if test "$enable_single" = "yes" ;then
		if test "$with_matlab" != "no" ;then
			AC_MSG_WARN([MATLAB linking is only possible without --enable-single. Disabled mex.])
			with_matlab=no
		fi
	fi

	if test "$with_matlab" != "no" ;then
		if test ${MEXDIR:0:1} = '~' ;then
			MEXDIR=$HOME:${MEXDIR:1}
		fi
		if test -x "${MEXDIR}/mex" ;then
			MEXDIR="${MEXDIR}"
		elif test -x "${MEXDIR}/bin/mex" ;then
			MEXDIR="${MEXDIR}/bin"
		else
			AC_MSG_ERROR([mex not found at $MEXDIR. Please reconfigure with --with-matlab=DIR to specify MATLAB location or --without-matlab to disable])
		fi
		echo MEXDIR=$MEXDIR
		MEXEXT="$($MEXDIR/mexext)"
		if test -z "$MEXEXT" ;then
			AC_MSG_ERROR([Unable to determine mex extension. reconfigure with --without-matlab to disable mex.])
		fi
		case $MEXEXT in
		mexglx)
			MEXLIB=glnx86
			MEXOPT="$OPTS -fPIC -m32 -fexceptions -D_FILE_OFFSET_BITS=64"
			#MEXOPT="${MEXOPT} -D_BSD_SOURCE -D_POSIX_C_SOURCE=200809L"
			mex_cpu=x86
			;;
		mexa64)
			MEXLIB=glnxa64
			MEXOPT="$OPTS -fPIC -m64 -fexceptions -fno-omit-frame-pointer"
			#MEXOPT="${MEXOPT} -D_BSD_SOURCE -D_POSIX_C_SOURCE=200809L"
			mex_cpu=x86_64
			;;
		mexmaci)
			MEXLIB=maci
			MEXOPT="$OPTS -fPIC -m32 -fexceptions"
			mex_cpu=x86
			;;
		mexmaci64)
			MEXLIB=maci64
			MEXOPT="$OPTS -fPIC -m64 -fexceptions"
			MEXOPT="${MEXOPT}"
			mex_cpu=x86_64
			;;
		mexw64)
			MEXLIB=win64
			MEXOPT=""
			mex_cpu=x86_64
			;;
		mexw32)
			MEXLIB=win32
			MEXOPT=""
			mex_cpu=x86
			;;
		*)
			AC_MSG_ERROR([Unknown mex extension "$MEXEXT"])
			;;
		esac
		if test "$mex_cpu" != "$host_cpu" ;then
			echo mex_cpu=$mex_cpu
			echo host_cpu=$host_cpu
			echo Architecture mismatch, disable mex support.
			mex_lib="no"
		else
			LDMEX="-L${MEXDIR}/${MEXLIB} -Wl,-rpath,${MEXDIR}/${MEXLIB} "
			AC_CHECK_LIB([mx], [mxGetPr], [mex_lib="yes"], [mex_lib="no"], [$LDMEX -shared])
		fi
		if test x$mex_lib = xyes ;then
			MEXOPT+=" -I${MEXDIR}/../extern/include -DMATLAB_MEX_FILE $OPENMP_CFLAGS"
			#2021-05-24: Linking with mwfftw3 crashes matlab. (mwfftw3 contains thread support)
			LIBMEX="$CHOL_LIBS -lmwlapack -lmwblas -lmx -lmex -lmat -lut -lstdc++ $LDEXE $FFTW_LIBS"
			LDMEX+="-no-fast-install -static -Xcompiler -shared"
		else
			with_matlab="no"
			AC_MSG_NOTICE([mex library test failed. Disable mex])
		fi
		AC_SUBST(MEXEXT)
		AC_SUBST(LIBMEX)
		AC_SUBST(LDMEX)
		AC_SUBST(MEXOPT)
	fi
	AM_CONDITIONAL(USE_MEX, [test "$with_matlab" != "no"])
	if test "$with_matlab" != "no" ;then
		echo "MEXDIR=$MEXDIR"
		echo "LIBMEX=$LIBMEX"
		echo "LDMEX=$LDMEX"
		echo "MEXOPT=$MEXOPT"
	fi
])

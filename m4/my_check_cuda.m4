AC_DEFUN([MY_CHECK_CUDA],[
	if test "${CN}"="CLANG" -a "$system" = "apple" ;then
		CUDAOPT="-stdlib=libc++"
	else
		CUDAOPT= #CFLAGS for CUDA
	fi

	#AUTO detect libraries in the system
	#Try to discover CUDA install in the system

	NVCC="$(which nvcc 2>/dev/null)"
	if test -z "$NVCC" ;then
		NVCC="/usr/local/cuda/bin/nvcc"
	fi
	if test -x "$NVCC" ;then
		if test -L "$NVCC" ;then
			NVCC="$(readlink "$NVCC")"
		fi
		CUDADIR_AUTO="$(dirname "$(dirname "${NVCC}")")"
	fi
	AC_ARG_ENABLE(cuda-double, AS_HELP_STRING([--enable-cuda-double],[Enable cuda double precision]), [cuda_double="$enableval"],[cuda_double="no"])
	AC_ARG_ENABLE(nvcc, AS_HELP_STRING([--enable-nvcc],[Always use nvcc to compile cu files]))
	AC_ARG_WITH(cuda, AS_HELP_STRING([--with-cuda[=DIR]],[With cuda support]), [CUDADIR=${withval}], [CUDADIR=])
	AC_ARG_WITH(ccbin, AS_HELP_STRING([--with-ccbin=CC,],[Specify CCBIN for nvcc]), [CCBIN=${withval}], [CCBIN=${CC}])

	if test -z "$CUDADIR" ;then #default value
		CUDADIR=$CUDADIR_AUTO
	elif test -n "$CUDADIR" ;then #user supplied options
		if test "$CUDADIR" != "no" ;then #user specified
			user_enabled_cuda=yes
		else
			user_enabled_cuda=no
		fi
	fi
	
	if test -z "$CUDADIR" -o "$CUDADIR" = "no" ;then
		with_cuda=no
	else
		with_cuda=yes
	fi
	if test "$with_cuda" != "no" ;then #test CUDA availability
		if test "$enable_single" = "yes" -a "$cuda_double" = "yes" ;then
			AC_MSG_ERROR([cuda_double is only possible without --enable-single])
		fi
		#if test "$enable_allstatic" = "yes"; then
		#AC_MSG_ERROR("all-static conflicts with cuda support. Please disable all-static")
		#fi
		AC_MSG_NOTICE([Checking CUDA])
		#Check for nvcc executable
		if test ${CUDADIR:0:1} = '~' ;then
			CUDADIR=$HOME${CUDADIR:1}
		fi
		NVCC="$CUDADIR/bin/nvcc"
		if test ! -x "$NVCC" -a -x "$CUDADIR/nvcc"; then
			NVCC="$CUDADIR/nvcc"
			CUDADIR="$(dirname $CUDADIR)"
		fi
		if test ! -x "$NVCC" ;then
			AC_MSG_NOTICE([nvcc not found])
			with_cuda="no"
		else
			cudaver="$($NVCC -V |grep release |cut -d ',' -f 2 |cut -d ' ' -f 3)"
			cudaver="${cudaver/.}"
			echo CUDA toolkit version is $cudaver
			if test "$cudaver" -lt 90 ;then
				cuarch=20 #minimum arch version supported.
			elif test "$cudaver" -lt 110 ;then
				cuarch=30
			elif test "$cudaver" -lt 120 ;then
				cuarch=35
			else
				cuarch=50
			fi
			if test "$cudaver" -lt 60 ;then
				cuarchmax=35 #maximum arch version supported.
			elif test "$cudaver" -lt 80 ;then
				cuarchmax=50
			elif test "$cudaver" -lt 90 ;then
				cuarchmax=60
			elif test "$cudaver" -lt 100 ;then
				cuarchmax=70
			elif test "$cudaver" -lt 110 ;then
				cuarchmax=75
			elif test "$cudaver" -lt 120 ;then
				cuarchmax=80
			else
				cuarchmax=90
			fi
			devarch=`nvidia-smi  --query-gpu=compute_cap --format=csv,noheader,nounits |sed  's/\.//g' |sort -n |head -1`
			if test -n "$devarch" ; then
				if test $devarch -lt $cuarch ;then
					AC_MSG_ERROR([CUDA version supports compute capability $cuarch - $cuarchmax, but device is at $devarch])
				elif test "$enable_release" != "yes" ;then
					if test $devarch -gt $cuarchmax ;then
						cuarch=$cuarchmax
					else
						cuarch=$devarch
					fi
				fi
			fi
		fi
	fi
	if test "$with_cuda" != "no" ;then
		#Check cor cuda library
		if test -f "${CUDADIR}/lib${libsuffix}/libcudart.$ldsuffix" ;then
			CUDA_L="${CUDADIR}/lib${libsuffix}"
		elif test -f "${CUDADIR}/lib/libcudart.$ldsuffix" ;then
			CUDA_L="${CUDADIR}/lib"
		elif test -f "${CUDADIR}/lib/x64/cudart.lib" ;then
			CUDA_L="${CUDADIR}/lib/x64"
		else
			CUDA_L=""
		fi
		if test -d "$CUDA_L" ;then
			LDCUDA="-L$CUDA_L -Wl,-rpath,$CUDA_L"
		else
			LDCUDA=
			AC_MSG_NOTICE([CUDA toolkit not found])
			with_cuda="no"
		fi
		AC_CHECK_LIB([cudart], [cudaSetDevice], [], [AC_MSG_NOTICE([libcudart not found]);with_cuda="no"], [$LIBCUDA $LDCUDA])

		#Check for cuda header
		if test -f "${CUDADIR}/include/cuda.h" ;then
			CUDAOPT+=" -I${CUDADIR}/include"
		else
			AC_CHECK_HEADERS([cuda.h], [] ,[AC_MSG_NOTICE([Header not found]);with_cuda="no"])
			CUDAOPT=
		fi
		if test "$with_cuda" != "no" ;then
			LIBCUDA+=" -lcurand -lcusolver -lcusparse -lcufft -lcublas -lcudart -lstdc++"
			if test "$enable_debug" = "yes" ;then
				CUDAOPT+=" -O0 -g"
			else
				CUDAOPT+=" -O3 -g"
			fi
			CUDAOPT+=" -DHAVE_CONFIG_H -I$BUILD_DIR "

			#CCBINOPT is what NVCC passes to C++ compiler
			CCBINOPT="$OPTScommon $OPTSprof $OPTS_CC $CPPFLAGS $OPENMP_CFLAGS ${CFLAGS/-std=c11/} -Wno-unused-parameter"
			CCBINOPT="${CCBINOPT/-fpermissive/}"
			case `$CCBIN --version 2>&1` in
			*clang*|*DPC++*)
				CCN=CLANG
				;;
			*icc*) #Old intel compiler
				CCN=ICC
				;;
			*)
				CCN=GCC
				CCBINOPT+=" -Wno-unused-value -D__STRICT_ANSI__ "
			esac
		fi
	fi
	if test "$with_cuda" != "no"  ;then
		if test $CN = CLANG -a x${enable_nvcc} != xyes ;then #use clang to compile GPU code. Which is as fast as nvcc but more tolerant.
			NVCC=$CC
			AC_MSG_NOTICE([Use $CC with cuda in directory $CUDADIR])
			if test $CC_VERSION -ge 50 ;then
				CUDAOPT+=" --std=c++17" #fix __float128 error. see https://discourse.llvm.org/t/error-float128-is-not-supported-on-this-target/72397
			fi
			CUDAOPT+=" --cuda-gpu-arch=sm_${cuarch} $CCBINOPT --cuda-path=$CUDADIR"
			CUDAOPT+=" --no-cuda-version-check -Wno-unknown-cuda-version -Wno-format-nonliteral"
			CCBINOPT= #no need for this
		else
			AC_MSG_NOTICE([Use nvcc with cuda in directory $CUDADIR])
			#NVCC+=" $CFLAGS_CPU"
			if test -n "$CCBIN" -a "$with_ccbin" != "no" ;then
				CCBIN="${CCBIN%% *}"
				NVCC+=" -ccbin ${CCBIN%% *}"
			fi
			if test "$cudaver" -lt 70 ;then
				#CUDA 4.0 supports compute_20
				#CUDA 7.0 supports C++11.
				AC_MSG_ERROR([Require at least CUDA toolkit 7.0])
			fi

			#virtual GPUs are used here to allow JIT compilation during runtime
			#set env CUDA_CACHE_MAXSIZE to a large number to cache JIT
			CUDAOPT+=" -lineinfo "
			CUDAOPT+=" -Wno-deprecated-gpu-targets"
			#To specify multiple architectures, use
			#CUDAOPT+=" -gencode arch=compute_20 -gencode arch=compute_30"

			if test -n "$CCBINOPT" ;then
				CCBINOPT="--compiler-options \"$CCBINOPT\""
			fi
		fi
		if test "$cuda_double" = "yes" ;then
			cuda_double=1
		else
			cuda_double=0		
		fi
		AC_SUBST(CUDAOPT)
		AC_SUBST(CCBINOPT)
		AC_SUBST(LIBCUDA)
		AC_SUBST(LDCUDA)
		AC_SUBST(NVCC)
	else
		cuda_double=0
		if test "$user_enabled_cuda" = "yes" ;then
			AC_MSG_NOTICE([Specified directory $CUDADIR does not contain CUDA toolkit.])
		elif test "$user_enabled_cuda" = "no" ;then
			AC_MSG_NOTICE([User disabled CUDA.])
		else
			AC_MSG_NOTICE([CUDA toolkit is not found.])
		fi
		cudaver=
	fi
	AC_DEFINE_UNQUOTED(CUDA_DOUBLE, [$cuda_double], [Use double precision mode in CUDA])
	AC_DEFINE_UNQUOTED(USE_CUDA, [${cudaver:-0}], [CUDA version])
	AM_CONDITIONAL(CUDA_DOUBLE, [test "$cuda_double" = "1" ])
	AM_CONDITIONAL(USE_CUDA, [test "$with_cuda" != "no"])
	if test "$with_cuda" != "no" ;then
		echo "NVCC=$NVCC"
		echo "CUDAOPT=$CUDAOPT"
		#echo "CCBINOPT=$CCBINOPT"
		#echo "LIBCUDA=$LIBCUDA"
		echo "LDCUDA=$LDCUDA"
	fi
])

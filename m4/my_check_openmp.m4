AC_DEFUN([MY_CHECK_OPENMP],[
	save_CFLAGS="$CFLAGS"
	has_omp="no"
	AC_ARG_ENABLE(openmp, AS_HELP_STRING([--enable-openmp], [Enable openmp support]))
	AS_IF([test "$enable_openmp" != "no" ],[
		AS_IF([test "$system" = "apple" ], [
			#MACOS supports openmp with separate omp library
			AC_CHECK_LIB([omp],[omp_set_num_threads],[has_omp=yes],[has_omp=no
				AS_IF([test "$enable_openmp" = "yes"], [ #force enable openmp on macos
					MY_DOWNLOAD([LIBOMP], [libomp${libsuffix}.tar.bz2], [${prefix}])
					unset ac_cv_lib_omp_omp_set_num_threads
					AC_CHECK_LIB([omp],[omp_set_num_threads],[has_omp=yes],[has_omp=no
						#ompurl="https://github.com/llvm/llvm-project/releases/download/llvmorg-13.0.0/openmp-13.0.0.src.tar.xz"
						MY_COMPILE([LIBOMP], ["libomp.tar.xz"], [${prefix} -DOPENMP_LIBDIR_SUFFIX=${libsuffix}])
						unset ac_cv_lib_omp_omp_set_num_threads
						AC_CHECK_LIB([omp],[omp_set_num_threads],[has_omp=yes],[has_omp=no])
					])
				])
			])
			AS_IF([test "$has_omp" = "yes"], [OPENMP_CFLAGS="-Xpreprocessor -fopenmp" OPENMP_LIBS=" -lomp" libtool_fix_openmp=1])
		], [AC_OPENMP])
		AS_IF([test -z "$OPENMP_CFLAGS" -a -n "$OPENMP_CXXFLAGS"], [OPENMP_CFLAGS="$OPENMP_CXXFLAGS"])
		#Don't add OPENMP_CFLAGS to CFLAGS yet to avoid building external library with the flag
		AS_IF([test -z "$OPENMP_CFLAGS"], [has_omp=no], [has_omp=yes])
	], [has_omp=no])
	AS_IF([test "$has_omp" = yes ], [
			AS_IF([test "$INTEL" = "1" ],[OPENMP_CFLAGS=${OPENMP_CFLAGS//-fopenmp/-qopenmp}])
			$1
		],[
			AS_IF([test "$INTEL" = "1" ],[CFLAGS+=" -pthread -qopenmp-simd"],[CFLAGS+=" -pthread -fopenmp-simd"])
			$2
		])
	echo OPENMP_CFLAGS=$OPENMP_CFLAGS
	echo OPENMP_LIBS=$OPENMP_LIBS
	CFLAGS=${save_CFLAGS}
	AM_CONDITIONAL(LIBTOOL_FIX_OPENMP, [test "$libtool_fix_openmp" = "1" -a -n "$LIBTOOL"])
])

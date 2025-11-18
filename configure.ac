AC_INIT(maos, 3.3)
AC_PREREQ(2.57)
AC_CONFIG_AUX_DIR([m4]) #for autoconf
AC_CONFIG_MACRO_DIR([m4]) #for aclocal
AC_PREFIX_DEFAULT(${HOME}/.aos) #do not use variable set in this file.
AM_INIT_AUTOMAKE([1.9.6 foreign -Wall no-define subdir-objects ])
#silent building
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])])
#The URL where misc supportive libaries are located.
BASEURL="https://github.com/lianqiw/files/raw/master/"

test "x$prefix" = xNONE && prefix=$ac_default_prefix #prefix is set too late in configure.ac by default
export PATH+=":${prefix}/bin" #add to PATH for doxygen if downloaded
CPPFLAGS+=" -I${prefix}/include " #use CPPFLAGS is used for preprocessor. Use it for -D, -I flags
SRC_DIR=$(cd "$srcdir" && pwd) #top_srcdir and abs_top_srcdir is empty. no idea why. but srcdir sometimes is the full path.
BUILD_DIR=$(pwd)
TMP_DIR="${BUILD_DIR}/external" #temporary directory for compiling
DOWNLOAD_DIR="${TMP_DIR}/downloads"

if ! which which >/dev/null 2>&1;then
	AC_MSG_ERROR([Please install 'which'])
fi
#Discover default download method
if which curl >/dev/null 2>&1;then
	wget="curl --retry 1 -L -s "
	wgetoutput="--output"
elif which wget >/dev/null 2>&1;then
	wget="wget -t1 --no-check-certificate"
	wgetoutput="-O"
else
	wget="echo Please install wget or curl && false"
fi

#If user does not supply any CFLAGS, set something so that AC_PROG_CC won't put a -g -O2 to CFLAGS
if test -z "$CFLAGS" ;then
	CFLAGS=" "
fi
CFLAGS="$CFLAGS -fPIC"
if test -z "$LIBS" ;then
	LIBS=" "
fi
if test -z "$CC" -a -n "$MAOS_CC" && which "$MAOS_CC" >/dev/null 2>&1;then
	CC="$MAOS_CC"
fi
AC_CONFIG_HEADERS([config.h]) #contains package name.
AC_CANONICAL_HOST #defines host host_os host_cpu host_vendor etc
AC_LANG([C])
AC_TYPE_SIZE_T #it sets CC if empty
AC_PROG_CC([$CC]) #overrides CC
#AC_LANG([C++]) #uncomment this line to use C++ compiler (for compatibility checking)
#The following must be after AC_PROG_CC
AC_PROG_MAKE_SET
AC_PROG_INSTALL
AC_CHECK_HEADERS([unistd.h])
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
LT_INIT([disable-shared, static]) #prefer static internal linking.
AC_C_BIGENDIAN([BIGENDIAN=1], [BIGENDIAN=0], AC_MSG_ERROR([Unknown Endianness]))
AC_DEFINE_UNQUOTED(BIGENDIAN, [$BIGENDIAN], [Endianness])
AC_DEFINE_UNQUOTED(SRCDIR,  ["$SRC_DIR"], [SOURCE folder])
AC_DEFINE_UNQUOTED(BUILDDIR,["$BUILD_DIR"],[BUILD folder])

#CFLAGS are for C compiler
#CXXFLAGS are for C++ compiler
LDOPTS=
LDEXE=
OPTS=

case "x$CC $CFLAGS" in
	x*++*)
		USE_CPP=1
		;;
	*)
		USE_CPP=0
esac
#find the real flavor of CC
case `$CC --version 2>&1` in
	*clang*)
		CN=CLANG
		;;
	*DPC++*) #New intel compiler with clang backend
	  	CN=CLANG
		LDEXE="$LDEXE -limf"
		INTEL=1
		;;
	*icpc*|*icc*) #Old intel compiler
		CN=ICC
		;;
	*g++*|*gcc*)
		CN=GCC
		;;
	*)
		CN=GCC
		AC_MSG_ERROR([Compiler $CC is unknown])
esac

AC_DEFINE_UNQUOTED(USE_CPP, [$USE_CPP], [Compile using C++ compiler])
AM_CONDITIONAL(GCC, [test x$CN = xGCC])
AM_CONDITIONAL(ICC, [test x$CN = xICC])
AM_CONDITIONAL(CLANG, [test x$CN = xCLANG])
CC_VERSION=$($CC -dumpversion)
AC_DEFINE_UNQUOTED(COMPILER, ["$CN (v$CC_VERSION)"], [Compiler name and version])

#Common configure options.
AC_ARG_ENABLE(debug, AS_HELP_STRING([--enable-debug],[Turn on debug]))
AC_ARG_ENABLE(gcov, AS_HELP_STRING([--enable-gcov],[Enable gcov for code coverage analysis.]))
AC_ARG_ENABLE(double, AS_HELP_STRING([--disable-double],[Use single precision for CPU and cuda.]))
AC_ARG_ENABLE(single, AS_HELP_STRING([--enable-single],[Use single precision for CPU and cuda.]))
AC_ARG_ENABLE(scheduler, AS_HELP_STRING([--disable-scheduler],[Disable built in scheduler]))
AC_ARG_ENABLE(release, AS_HELP_STRING([--enable-release], [Compile for enhanced compatibility]))
AC_ARG_ENABLE(doxygen, AS_HELP_STRING([--enable-doxygen], [Force enable doxygen. It will download doxygen binary if not available in the system]))
AC_ARG_ENABLE(binary, AS_HELP_STRING([--disable-binary], [Disable downloading precompiled library]))
AC_ARG_ENABLE(download, AS_HELP_STRING([--disable-download], [Disable downloading library]))
AC_ARG_ENABLE(openmp, AS_HELP_STRING([--enable-openmp], [Enable openmp support]))
AC_ARG_ENABLE(threadlib, AS_HELP_STRING([--enable-threadlib], [Enable threading library for blas and fftw.]))
#echo enable_openmp="$enable_openmp"
#echo openmp_default="$openmp_default"

if test "$enable_double" = "no" -o "$enable_single" = "yes" ;then
	enable_single=1
else
	enable_single=0
fi
AC_DEFINE_UNQUOTED(CPU_SINGLE, $enable_single, [Single precision CPU mode])
AM_CONDITIONAL(GCOV, [ test x$enable_gcov = "xyes" ])
if test x$enable_gcov = "xyes" ;then
	OPTSprof+="  -fprofile-arcs -ftest-coverage --coverage"
	enable_debug="yes"
fi
if test "$enable_release" = "yes" ;then
	enable_scheduler="no"
fi
AM_CONDITIONAL(USE_SCHEDULER, [test "$enable_scheduler" != "no"])
if test "$enable_scheduler" = "no" ;then
	AC_DEFINE(MAOS_DISABLE_SCHEDULER, [1], [Disable built-in scheduler])
	enable_lws="no"
fi
#Specifying both -O3 and -g in mac produces internal error in gcc 4.2
#-all-static will produce a static exe if all libraries are available in static linking mode.
#use -Wno-error=xxxx to disable certain warnings from being treated as errors.
#use -Wno-xxxx to disable certain warnings entirely.
OPTScommon=

#OPTScommon+=" -Werror -Wno-error=unknown-pragmas -Wformat=2 -Wno-error=deprecated-declarations -Wno-error=narrowing"

COMMON_ALL="-Wall -Wshadow -Wextra -Wno-missing-braces -Wno-missing-field-initializers "
COMMON_ICC="-wd193 -wd1418 -wd2259 -wd1372 -wd981 -wd1292 -wd10237 -wd858 -wd1478 -wd1786 -wd3180 -diag-disable remark -vec-report=3"
COMMON_CLANG=" -Wno-pass-failed -Wno-sign-compare -Wno-gnu-folding-constant"  #CLANG has false positive in detecting sign comparison.
#Set optimize flags for different compilers
OPTIM_ALL="-g -O3 -mtune=native "

if test "$INTEL" = "1" ;then
	OPTIM_CLANG+=" -fp-model=precise"
fi

DEBUG_ALL="-g3 -O0" #-og
DEBUG_GCC="-ggdb -fno-omit-frame-pointer -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls " #-fsanitize=address"
if test "$USE_CPP" = "1" ;then
	COMMON_GCC+=" -Wno-changes-meaning -Wno-class-memaccess -Wno-error=use-after-free -Wno-format-truncation -Wno-error=int-in-bool-context -fpermissive"
	COMMON_CLANG+=" -Wno-vla-cxx-extension  -Wno-deprecated -Wno-address-of-temporary"
	COMMON_ALL+=" -Wno-deprecated-declarations " #do not emit warning
fi
#inline definition:

#in GCC, extern inline will use the inline and not emit function definition.
#in GCC, inline will cause use the inline, but also emit function deminition, causing multiple definitions if put in header.
#in C99, extern inline will cause use the inline, but also emit function deminition, causing multiple definitions if put in header. This is equivalent to GCC inline
#in C99, inline will use the inline and not emit function definition. This is equivalent to GCC extern inline
#GCC mode can be forced on using __attribute__((gnu_inline))
#always_inline causes inline to be done in debug mode.

if test "x$enable_debug" = "xyes" ;then
	AC_DEFINE(DEBUG,[1],[Compile maos with debug information])
	CPPFLAGS+=" -DDEBUG_TRACE"
	ccprefix=DEBUG
else
	AC_DEFINE(DEBUG,[0],[Compile maos with optimization])
	ccprefix=OPTIM
fi
AM_CONDITIONAL(DEBUG, [test "x$enable_debug" = "xyes"])
optname1=COMMON_ALL
optname2=COMMON_${CN}
optname3=${ccprefix}_ALL
optname4=${ccprefix}_${CN}
#echo OPTS_CC from "${optname1} ${optname2} ${optname3} ${optname4}"
OPTS_CC="${!optname1} ${!optname2} ${!optname3} ${!optname4}"
#echo OPTS_CC=$OPTS_CC
unset optname1 optname2 optname3 optname4

#Check operating system.
AM_CONDITIONAL(LINUX,  false)
AM_CONDITIONAL(BSD,    false)
AM_CONDITIONAL(DARWIN, false)
AM_CONDITIONAL(WIN,    false)

#echo host_os=$host_os
case "$host_os" in
	*linux*)
		AM_CONDITIONAL(LINUX, true)
		LDEXE="$LDEXE -lrt" #need real time library.
		#RUNPATH is not used for dlopen or secondary dependencies. RPATH is. disable-new-dtags uses RPATH instead of RUNPATH
		#RUNPATH is searched after LD_LIBRARY_PATH while RPATH is searched before LD_LIBRARY_PATH. Both are before the default paths.
		LDOPTS+=" -Wl,--no-as-needed,--disable-new-dtags" #as-needed is default in some system (ubuntu) and cause mkl library checking to fail
		CPPFLAGS+=" -D_DEFAULT_SOURCE" #this is needed to compile libwebsockets and to use certain symbols
		case ${host_cpu} in
			aarch64|arm64)
				libsuffix="arm64"
			;;
			arm|arm32)
				libsuffix="arm32"
			;;
			x86_64|amd64)
				libsuffix="64"
			;;
			*86|*32)
				libsuffix="32"
			;;
			*)
				libsuffix="${host_cpu}"
			;;
		esac
		
		system="linux"
		ldsuffix="so"
		LDSO="-Wl,--no-undefined" #force an error if symbols are not resolved.
	;;
	*darwin*)
		AM_CONDITIONAL(DARWIN, true)
		system="apple"
		ldsuffix="dylib"
		if test "$enable_release" = "yes" ;then
			CFLAGS+=" -mmacosx-version-min=10.13" #_open_memstream is only available on macos 10.13 and later
		fi
		if test "$host_cpu" = "aarch64" ;then
			libsuffix=a64 #arm mac
		else
			libsuffix=fat #intel mac
		fi
		LDSO="-Wl,-undefined,error" #force an error if symbols are not resolved.
	;;
	*bsd*)
		AM_CONDITIONAL(BSD, true)
		system="bsd"
		ldsuffix="so"
	;;
	*cygwin*|*mingw*)
		AM_CONDITIONAL(WIN, true)
		system="win"
		ldsuffix="dll"
	;;
	*)
		AC_MSG_ERROR([Unrecognized OS])
esac

case $host_cpu in
	*64)
	;;
	*86|*32|arm)
		CPPFLAGS+=" -D_FILE_OFFSET_BITS=64" #to fix error about stat() in 32bit.
	;;
	*)
		AC_MSG_ERROR([Unknown host_cpu])
	;;
esac

AC_DEFINE_UNQUOTED(LDSUFFIX, ["$ldsuffix"], [Suffix of dynamic load library])
libdir="${prefix}/lib${libsuffix}" #external library installation location
mkdir -p "${libdir}"
export LD_LIBRARY_PATH+=":${libdir}" #avoid using rpath in configure since it will be hard coded into 3rd library during compilation
LDFLAGS+=" -L${libdir} "
m4_ifdef([PKG_PROG_PKG_CONFIG], [PKG_PROG_PKG_CONFIG])
#compile LIBNAME srcfile.tar.gz [compile options]
AS_IF([test "$enable_openmp" != "no" ],[
	AS_IF([test "$system" = "apple" ], [
		#MACOS supports openmp with separate omp library
		AC_CHECK_LIB([omp],[omp_set_num_threads],[has_omp=1],[has_omp=0
			AS_IF([test "$enable_openmp" = "yes"], [ #force enable openmp on macos
				MY_DOWNLOAD([LIBOMP], [libomp${libsuffix}.tar.bz2], [${prefix}])
				unset ac_cv_lib_omp_omp_set_num_threads
				AC_CHECK_LIB([omp],[omp_set_num_threads],[has_omp=1],[has_omp=0
					#ompurl="https://github.com/llvm/llvm-project/releases/download/llvmorg-13.0.0/openmp-13.0.0.src.tar.xz"
					MY_COMPILE([LIBOMP], ["libomp.tar.xz"], [${prefix} -DOPENMP_LIBDIR_SUFFIX=${libsuffix}])
					unset ac_cv_lib_omp_omp_set_num_threads
					AC_CHECK_LIB([omp],[omp_set_num_threads],[has_omp=1],[has_omp=0])
				])
			])
		])
		AS_IF([test "$has_omp" = "1"], [OPENMP_CFLAGS="-Xpreprocessor -fopenmp" LDFLAGS+=" -lomp" libtool_fix_openmp=1])
	], [AC_OPENMP])
	AS_IF([test -z "$OPENMP_CFLAGS" -a -n "$OPENMP_CXXFLAGS"], [OPENMP_CFLAGS="$OPENMP_CXXFLAGS"])
	#Don't add OPENMP_CFLAGS to CFLAGS yet to avoid building external library with the flag
	AS_IF([test -z "$OPENMP_CFLAGS"], [has_omp=0], [has_omp=1])
], [has_omp=0])
AS_IF([test "$has_omp" = 1], 
	[AS_IF([test "$INTEL" = "1" ],[OPENMP_CFLAGS=${OPENMP_CFLAGS//-fopenmp/-qopenmp}])],
	[AS_IF([test "$INTEL" = "1" ],[CFLAGS+=" -pthread -qopenmp-simd"],[CFLAGS+=" -pthread -fopenmp-simd"])])
AM_CONDITIONAL(LIBTOOL_FIX_OPENMP, [test "$libtool_fix_openmp" = "1" -a -n "$LIBTOOL"])
#Check for zlib
AC_CHECK_LIB([z],[gzopen],[:],[
	MY_COMPILE([ZLIB], [https://www.zlib.net/current/zlib.tar.gz])
	unset ac_cv_lib_z_gzopen
	AC_CHECK_LIB([z],[gzopen],[:],[AC_MSG_ERROR([Please install zlib.])])
])

#Must be before FFT and cholmod checking
MY_CHECK_LAPACK(
	[AC_MSG_NOTICE([Found blas/lapack])],
	[AC_MSG_ERROR([Please install blas/lapack])]
)
AX_CHECK_OPENSSL(
	[AC_DEFINE([HAVE_OPENSSL], [1], [Define if OpenSSL is available])],
	[AC_MSG_NOTICE([OpenSSL library not found])]
)
MY_CHECK_FFT([:],[AC_MSG_ERROR([FFTW library is not found])]) #must be before CHECK_MEX
MY_CHECK_CHOLMOD([:],[AC_MSG_ERROR([CHOLMOD library is not found])])
MY_CHECK_LWS([AC_DEFINE(HAS_LWS, [1], [libwebsockets is not found])],[:])
MY_CHECK_CMOCKA([AC_DEFINE(HAS_CMOCKA, [1], [Enable unit testing with cmocka])],[:])
MY_CHECK_MEX()
MY_CHECK_CUDA()
MY_CHECK_GTK()
MY_CHECK_ZSTD()
#Doxygen Features
if test "$enable_doxygen" = "yes" ;then
	if ! which doxygen ;then
		MY_DOWNLOAD([DOXYGEN], [doxygen${libsuffix}.tar.bz2], [${prefix}], [1])
	fi

	if ! which pdflatex ;then
		$wget "https://yihui.org/gh/tinytex/tools/install-unx.sh" $wgetoutput - | sh
	fi
fi
DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)
DX_INIT_DOXYGEN($PACKAGE_NAME, doxygen.cfg)
AC_CONFIG_FILES( [Makefile
	sys/Makefile
	math/Makefile
	lib/Makefile
	cuda/math/Makefile
	cuda/sim/Makefile
	cuda/recon/Makefile
	cuda/test/Makefile
	cuda/Makefile
	maos/Makefile
	skyc/Makefile
	tools/Makefile
	test/Makefile
	mex/Makefile
])
OPTS+=" $OPTScommon $OPTS_CC $OPTSprof "
#put rpath here to avoid using it when compiling 3rd party library
LDFLAGS+=" $LDOPTS -Wl,-rpath,${libdir}" #to run artifacts from github workflow avoid using rpath
CFLAGS+=" $FFTW_CFLAGS $OPTS $OPENMP_CFLAGS" #do not replace 
LDEXE+=" -lz -lm -ldl"
LDEXE_FULL="$FFTW_LIBS $CHOL_LIBS $LAPACK $LDEXE" #Make sure FFTW_LIBS is before LAPACK to avoid MKL providing FFTW functions.
#echo LIBS=$LIBS
LIBS="" #somehow LIBS contain -fftw3 -lnvidia-smlA
AC_SUBST(ldsuffix)
AC_SUBST(CN)
AS_UNSET(OPTS)
AC_SUBST(LDEXE)
AC_SUBST(LDEXE_FULL)
AC_SUBST(SRC_DIR)
AC_SUBST(BUILD_DIR)
AC_SUBST(LIBS)
AC_SUBST(LDSO) #to link .so with all symbols resolved
AC_OUTPUT
echo SRC_DIR=$SRC_DIR
echo BUILD_DIR=$BUILD_DIR
echo BIGENDIAN=$BIGENDIAN
echo CC=$CC
echo USE_CPP=$USE_CPP
echo CPPFLAGS=$CPPFLAGS
echo CFLAGS=$CFLAGS
echo LDEXE=$LDEXE
echo LDFLAGS=$LDFLAGS

yesno=("no" "yes")
FP=("FP64" "FP32")
echo
echo Compiling options:
echo
echo "host_os:         ${host_os}"
echo "host_cpu:        ${host_cpu}"
echo "CC:              ${CN} ($CC_VERSION)"
echo "Debug:           ${enable_debug:-no}"
echo "Scheduler:       ${enable_scheduler:-yes}"
echo "CMOCKA:          ${has_cmocka}"
echo "Coverage:        ${enable_gcov:-no}"
echo "OpenMP:          ${yesno[[$has_omp]]}"
echo "FFTW Threads:    ${yesno[[$fftw_threads]]}"
echo "Precision:       ${FP[[$enable_single]]}"
echo "CUDA:            ${cudaver:-no} ${cudaver:+ (${FP[[1-$cuda_double]]}, arch: $cuarch)}"
echo "GTK:             ${gtk_ver}"
echo "MATLAB MEX:      ${with_matlab}"
echo "Doxygen:         ${DX_DOXYGEN:-no}${DX_DOXYGEN:+ (make doc)}"
echo

echo "*************************************************************************"
echo "*** MAOS Documentation can be found at https://lianqiw.github.io/maos ***"
echo "*** MAOS Source can be obtained at https://github.com/lianqiw/maos    ***"
echo "*************************************************************************"

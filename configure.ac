AC_INIT(maos, 2.6.0, lianqiw@gmail.com)
AC_PREREQ(2.57)
AC_CONFIG_AUX_DIR([m4])
AM_INIT_AUTOMAKE([1.9.6 foreign -Wall no-define subdir-objects ])
#silent building
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])]) 
#The URL where misc supportive libaries are located.
BASEURL="https://github.com/lianqiw/files/raw/master/"

#If user does not supply any CFLAGS, set something so that AC_PROG_CC won't put a -g -O2 to CFLAGS
if test -z "$CFLAGS" ;then
    CFLAGS=" "
fi
if test -z "$LIBS" ;then
    LIBS=" "
fi
if test -z "$CC" ;then 
    CC="$MAOS_CC"
fi
AC_CONFIG_MACRO_DIR([m4])
#After enable HEADER, the package name will not appear in CFLAGS. Appear in config.h instead.
AC_CONFIG_HEADERS([config.h])
#defines host host_os host_cpu host_vendor etc
AC_CANONICAL_HOST
AC_LANG([C])
AC_TYPE_SIZE_T #it sets CC if empty
AC_PROG_CC([$CC])
#AM_PROG_CC_C_O #deprecated. Included in AC_PROG_CC.
#AC_PROG_CC_STDC #enable the current c standard, which is c99 circa 2014
AC_PROG_CC_C99  #enable C99
#AC_PROG_CXX([$CXX])
AC_OPENMP
#The following must be after AC_PROG_CC
AC_PROG_MAKE_SET
AC_PROG_INSTALL
AC_DISABLE_SHARED
AC_ENABLE_STATIC
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
AC_PROG_LIBTOOL

AC_C_BIGENDIAN([BIGENDIAN=1], [BIGENDIAN=0], AC_MSG_ERROR("Unknown Endianness"))
AC_DEFINE_UNQUOTED(BIGENDIAN, [$BIGENDIAN], ["Endianness"])

#top_srcdir and abs_top_srcdir is empty. no idea why. but srcdir sometimes is the full path.
SRC_DIR=$(cd "$srcdir" && pwd)
BUILD_DIR=$(pwd)
#AC_DEFINE_UNQUOTED means expand the $ variable
AC_DEFINE_UNQUOTED(SRCDIR,  ["$SRC_DIR"], ["Define SRC folder"])
AC_DEFINE_UNQUOTED(BUILDDIR,["$BUILD_DIR"],["Define BUILD folder"])

die(){
    shift
    AC_MSG_ERROR(["Error:$@"])
}

USE_CPP=0
#find the real flavor of CC
CCv=`$CC -v 2>&1`
if grep clang++ <<< $CCv ;then
    CN=CLANG
    USE_CPP=1
elif grep clang <<< $CCv ;then
    CN=CLANG
elif grep icpc <<< $CCv ;then
    CN=ICC
    USE_CPP=1
elif grep icc <<< $CCv ;then
    CN=ICC
elif grep g++ <<< $CCv ;then
    CN=GCC
    USE_CPP=1
elif grep gcc <<< $CCv ;then
    CN=GCC
else
    CN=GCC
    AC_MSG_ERROR([Compiler $CC is unknown])
fi
AC_DEFINE_UNQUOTED(USE_CPP, [$USE_CPP], ["Whether MAOS is using C++ compiler"])
#canonicalize the host_cpu name
case "$host_cpu" in
    *86)
	host_cpu=x86
	;;
    *64)
	host_cpu=x86_64
	;;
esac
#AUTO detect libraries in the system
#Try to discover CUDA install in the system
if ! which which >/dev/null 2>&1;then
    die "Please install 'which'"
fi
if test -n "`which nvcc 2>/dev/null`" ;then
    NVCC="$(which nvcc)"
    if test -L "$NVCC" ;then
	NVCC="$(readlink "$NVCC")"
    fi
else
    for cuda in /etc/alternatives/cuda /usr/local/cuda /usr/local/cuda-* /opt/cuda /opt/cuda-* /opt/cuda/default /opt/cuda/cuda-*;do
		if test -x "${cuda}/bin/nvcc" ;then
			NVCC="${cuda}/bin/nvcc"
			break;
		fi
    done
    unset nvcc
fi
if test -n "$NVCC" ;then
    CUDADIR_AUTO="$(dirname "$(dirname "${NVCC}")")"
    unset NVCC
fi

#Compile mex routines if MATLAB is detected and single is not enabled
MEXDIR_AUTO=
if test -n "`which matlab 2>/dev/null`" ;then
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

#Common configure options.
AC_ARG_ENABLE(allstatic, AS_HELP_STRING([--enable-allstatic], [Compile statically]))
AC_ARG_ENABLE(debug, AS_HELP_STRING([--enable-debug],[Turn on debug]))
AC_ARG_ENABLE(gcov, AS_HELP_STRING([--enable-gcov],[Enable gcov for code coverage analysis.]))
AC_ARG_ENABLE(single, AS_HELP_STRING([--enable-single],[Use single precision for CPU and cuda.]))
AC_ARG_ENABLE(cuda-double, AS_HELP_STRING([--enable-cuda-double],[Enable cuda double precision]), [cuda_double="$enableval"],[cuda_double="no"])
AC_ARG_ENABLE(scheduler, AS_HELP_STRING([--disable-scheduler],["Disable built in scheduler"]))
AC_ARG_ENABLE(websockets, AS_HELP_STRING([--disable-websockets],[Disable websockets]))
AC_ARG_ENABLE(mkl, AS_HELP_STRING([--enable-mkl], ["Prefer to use MKL"]), [use_mkl="$enableval"], [use_mkl="default"])
AC_ARG_WITH(cuda, AS_HELP_STRING([--with-cuda[=DIR]],[With cuda support]), [CUDADIR=${withval}], [CUDADIR="$CUDADIR_AUTO"])
AC_ARG_WITH(ccbin, AS_HELP_STRING([--with-ccbin=CC,],[Specify CCBIN for nvcc]), [CCBIN=${withval}], [CCBIN=${CC}])
AC_ARG_WITH(gtk,  AS_HELP_STRING([--with-gtk[=DIR]],[Enable gtk support]), [GTKDIR=${withval}], [GTKDIR=])
AC_ARG_ENABLE(gtk-3, AS_HELP_STRING([--enable-gtk-3],[Using GTK+-3]), [use_gtk_3="$enableval"], [use_gtk_3="yes"])
AC_ARG_ENABLE(gtk-4, AS_HELP_STRING([--enable-gtk-4],[Using GTK4]), [use_gtk_4="$enableval"], [use_gtk_4="no"])
AC_ARG_WITH(matlab, AS_HELP_STRING([--with-matlab[=DIR]],[Enable compiling mex routines for matlab]), [MEXDIR=${withval}], [MEXDIR=${MEXDIR_AUTO}])
AC_ARG_ENABLE(release, AS_HELP_STRING([--enable-release], [Compile for broad compatibility]))
AC_ARG_ENABLE(doxygen, AS_HELP_STRING([--enable-doxygen], [Download doxygen binary if not available in the systemy]))
AC_ARG_ENABLE(cmocka, AS_HELP_STRING([--enable-cmocka], [Enable test with cmocka]))
AC_ARG_ENABLE(binary, AS_HELP_STRING([--disable-binary], [Disable downlong precompiled library]))
if test -z "$MEXDIR" -o "$MEXDIR" = "no" ;then
    with_matlab=no
else
    with_matlab=yes
fi
if test -z "$CUDADIR" -o "$CUDADIR" = "no" ;then
    with_cuda=no
else
    with_cuda=yes
fi
if test "$enable_single" = "yes" ;then
    if test "$cuda_double" = "yes" ;then
        AC_MSG_ERROR(["cuda_double is only possible without --enable-single"])
    fi
    if test "$with_matlab" != "no" ;then
		AC_MSG_WARN(["MATLAB linking is only possible without --enable-single. Disabled mex."])
		with_matlab=no
    fi
    enable_single=1
else
    enable_single=0
fi
AC_DEFINE_UNQUOTED(CPU_SINGLE, $enable_single, "Use single precision mode")

if test x$enable_gcov = "xyes" ;then
    OPTSprof+="  -fprofile-arcs -ftest-coverage --coverage"
    debug="yes"
fi
if test "$enable_release" = "yes" ;then
    enable_scheduler="no"
    enable_websockets="no"
	#enable_allstatic="yes"
fi
#Specifying both -O3 and -g in mac produces internal error in gcc 4.2
#-all-static will produce a static exe if all libraries are available in static linking mode.
OPTScommon=
if test "${USER}" = "lianqiw" ;then
    OPTScommon+=" -Werror -Wno-error=unknown-pragmas -Wno-unknown-pragmas -Wformat=2"
    #OPTScommon+=" -Wno-unknown-pragmas -Wno-unused-function -fno-common -Wno-missing-field-initializers -Wno-format-security -Wno-vla " 
fi
LDEXE=
LDRT=
OPTS=
COMMON_ALL="-Wall -Wshadow -Wextra " 
#Set optimize flags for different compilers
#-xO Enables SSE3, SSE2 and SSE instruction sets optimizations for non-Intel CPUs 
COMMON_ICC="-wd193 -wd1418 -wd2259 -wd1372 -wd981 -wd1292 -wd10237 -wd858 -wd1478 -wd1786 -wd3180 -diag-disable remark -vec-report=3"
COMMON_GCC+=" -Wno-missing-braces -Wno-missing-field-initializers" #gcc-4.8 error
COMMON_CLANG+=" -Wno-sign-compare" #CLANG has false positive in detecting sign comparison.

OPTIM_ALL="-g -O3 -mtune=native "
OPTIM_GCC="-ftree-vectorize -ffast-math -fno-finite-math-only "

DEBUG_ALL="-g3 -Og -DDEBUG_TRACE"
DEBUG_GCC="-ggdb -fno-omit-frame-pointer -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls " #-fsanitize=address"

#inline definition:

#in GCC, extern inline will use the inline and not emit function definition.
#in GCC, inline will cause use the inline, but also emit function deminition, causing multiple definitions if put in header.
#in C99, extern inline will cause use the inline, but also emit function deminition, causing multiple definitions if put in header. This is equivalent to GCC inline
#in C99, inline will use the inline and not emit function definition. This is equivalent to GCC extern inline 
#GCC mode can be forced on using __attribute__((gnu_inline))
#always_inline causes inline to be done in debug mode.
AM_CONDITIONAL(GCC, [test x$CN = xGCC])
AM_CONDITIONAL(ICC, [test x$CN = xICC])
AM_CONDITIONAL(CLANG, [test x$CN = xCLANG])
case "$CN" in
    ICC)
	CC_VERSION=$($CC -dumpversion)
	if test -n "$CC_VERSION" ;then
	    major=${CC_VERSION%%.*}
	    minor=${CC_VERSION#*.}
	    minor=${minor%%.*}
	    if test $major -lt 11 ;then
			AC_MSG_ERROR([OpenMP is not supported by ICC $CC_VERSION])
			unset OPENMP_CFLAGS
	    fi
	fi
	;;
    GCC)
	#Test GCC version
	CC_VERSION=$($CC -dumpversion)
	if test -n "$CC_VERSION" ;then
	    major=${CC_VERSION%%.*}
	    minor=${CC_VERSION#*.}
	    if test $minor = $CC_VERSION ;then
			minor=0
	    else
			minor=${minor%%.*}
	    fi
	    echo major=$major minor=$minor
	    if test "$major" -lt 4 ;then
			AC_MSG_ERROR([At least version 4 is required when compiling with gcc.])
	    fi

	    #GCC introduced openmp task in 4.4
	    if test $major -le 4 -a $minor -le 4 ;then
			AC_MSG_WARN([OpenMP is not supported by GCC $CC_VERSION])
			unset OPENMP_CFLAGS
	    fi
	fi
	;;
    CLANG)
		CC_VERSION=$($CC -dumpversion)
		if test "$system" = "apple" ;then
			CFLAGS+=" -stdlib=libc++"
		fi
	;;
    *)
	;;	
esac

AC_DEFINE_UNQUOTED(COMPILER, ["$CN ($CC_VERSION)"], ["Compiler"])
if test "x$enable_debug" = "xyes" ;then
    AM_CONDITIONAL(DEBUG,true)
    AC_DEFINE(DEBUG,[1],["Debug the code"])
    ccprefix=DEBUG
else
    AM_CONDITIONAL(DEBUG,false)
    AC_DEFINE(DEBUG,[0],["Do not debug"])
    ccprefix=OPTIM
fi
optname1=COMMON_ALL
optname2=COMMON_${CN}
optname3=${ccprefix}_ALL
optname4=${ccprefix}_${CN}
echo OPTScc from "${optname1} ${optname2} ${optname3} ${optname4}"
OPTScc="${!optname1} ${!optname2} ${!optname3} ${!optname4}"
echo OPTScc=$OPTScc
unset optname1 optname2 optname3 optname4

#Used for embedding conf folder in maos exe
AC_SUBST(objcopy_output)
AC_SUBST(objcopy_arch)

#Check operating system.
AM_CONDITIONAL(LINUX,  false)
AM_CONDITIONAL(BSD,    false)
AM_CONDITIONAL(DARWIN, false)
AM_CONDITIONAL(WIN,    false)
mkl_avail="no"
echo host_os=$host_os
case "$host_os" in
    *linux*)
		AM_CONDITIONAL(LINUX, true)
		LDRT=" -lrt" #need real time library.
		LDOPTS+=" -Wl,--no-as-needed" #as-needed is default in some system (ubuntu) and cause mkl library checking to fail
		CFLAGS+=" -fPIC"
		#CFLAGS_EXTRA="-g -D_XOPEN_SOURCE=700 -D_BSD_SOURCE -D_DEFAULT_SOURCE -D_POSIX_C_SOURCE=200809L"
		if test "$enable_release" = "yes" ;then
			CFLAGS+=" -D_GLIBC_2_10_SOURCE" #Does not seem to be effective
		fi
		if test "${host_cpu}" = "x86_64" ;then
			libsuffix="64"
		else
			libsuffix="32"
		fi
		mkl_avail="yes"
		system="linux"
		ldsuffix="so"
		LDSO="-Wl,--no-undefined" #force an error if symbols are not resolved.
	;;
    *darwin*)
		AM_CONDITIONAL(DARWIN, true)
		CFLAGS+=" -fPIC"
		CFLAGS_EXTRA=
		mkl_avail="no"
		libsuffix=fat
		system="apple"
		ldsuffix="dylib"
		#for compatibility with OS X >=10.6. Don't use 10.5 which will render mex files not usable.
		if test "$enable_release" = "yes" ;then
			CFLAGS+=" -mmacosx-version-min=10.7"
		fi
		if test "$CN" = "CLANG" ;then
			CFLAGS+=" -stdlib=libc++"
		fi
		LDSO="-Wl,-undefined,error" #force an error if symbols are not resolved.
		#MACOS supports openmp with separate omp library
		#AC_CHECK_LIB([omp],[omp_set_num_threads],[OPENMP_CFLAGS="-Xpreprocessor -fopenmp" LDFLAGS+="-lomp" libtool_fix_openmp=1],[])
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
    none)
		AC_MSG_NOTICE([OS is none])
	;;
    *)
	echo host_os=$host_os
	AC_MSG_ERROR([Unrecognized OS])
esac
AM_CONDITIONAL(LIBTOOL_FIX_OPENMP, [test "$libtool_fix_openmp" = "1"])
if test -z "$OPENMP_CFLAGS" ;then
    CFLAGS+=" -pthread" #enables pthreads feature and links in the library
    has_openmp=0
else
    CFLAGS+=" $OPENMP_CFLAGS"
    has_openmp=1
fi
AC_DEFINE_UNQUOTED(HAS_OPENMP,[$has_openmp],["Defined to 1 if Openmp is enabled"])

AC_DEFINE_UNQUOTED(LDSUFFIX, ["$ldsuffix"], "Suffix of dynamic load library")
PRE_DIR="${HOME}/.aos"          #prefix for installation
TMP_DIR="${BUILD_DIR}/external" #temporary directory for compiling
LIB_DIR="${PRE_DIR}/lib${libsuffix}" #external library
LDFLAGS+=" -L${LIB_DIR} -Wl,-rpath,${LIB_DIR}" #-avoid-version
CPPFLAGS+=" -I${PRE_DIR}/include " #use CPPFLAGS so that the preprocessor can also use it.
PATH+=":${PRE_DIR}/bin"
DOWNLOAD_DIR="${TMP_DIR}/downloads"
OPTS+=" $CFLAGS_EXTRA "
if test "$system" = "win" ;then
    LDFLAGS+=" -Wl,-rpath,${PREDIR}/bin"
fi
if test "$host_cpu" = "x86_64" ;then
    CFLAGS+=" -DDLONG "
    mflag="-m64"
else
    CFLAGS+=" -D_FILE_OFFSET_BITS=64" #to fix error about stat() in 32bit.
    mflag="-m32"
fi
if ! echo $CFLAGS | grep -q -- "$mflag" ;then #don't repeat mflag
    CFLAGS+=" $mflag"
fi
unset mflag

#Discover default download method
if test -x "`which curl 2>/dev/null`" ;then
    wget="curl --retry 1 -L -s "
    wgetoutput="--output"
elif test -x "`which wget 2>/dev/null`" ;then
    wget="wget -t1 --no-check-certificate"
    wgetoutput="-O"
else
    wget="die wget/curl not found ; "
fi


AM_CONDITIONAL(USE_SCHEDULER, true)
if test "$enable_scheduler" = "no" ;then
    AC_DEFINE(MAOS_DISABLE_SCHEDULER, [1], "Disable built in scheduler")
    AM_CONDITIONAL(USE_SCHEDULER, false)
    enable_websockets="no"
fi
if test "$host_cpu" = "x86_64" ;then 
    chol_fun=cholmod_l_factorize
    AM_CONDITIONAL(USE_LONG, true)
else
    chol_fun=cholmod_factorize
    AM_CONDITIONAL(USE_LONG, false)
fi

LDEXE+="$LDRT " #-static" #show that we prefer static linking.
if test "x$enable_allstatic" = "xyes" ;then
    LDEXE_MAOS="-all-static" #only all static linking for maos
    AM_CONDITIONAL(USE_STATIC, true)
    AC_DEFINE(USE_STATIC, [1], "We are compling maos statically")
	mkl="mkl-static mkl" #prefer static
else
    AM_CONDITIONAL(USE_STATIC, false)
    AC_DEFINE(USE_STATIC, [0], "We are compling maos dynamically")
	mkl="mkl mkl-static" #prefer dynamic
fi

#call convention of fortran libraries
#AC_F77_FUNC(dgemm)
#AC_F77_FUNC(dposv)
#download LIBNANE file.tar.bz2 [outdir]
extract(){
    fn="$1"
    OUTDIR="$2"
    if test -z "$OUTDIR" ;then
	AC_MSG_ERROR([Please specify third argument])
	return 1
    else
	echo Extracting "${DOWNLOAD_DIR}/${fn}" to "${OUTDIR}"
	mkdir -p "${OUTDIR}" || die "failed to mkdir"
	if ! tar -xf "${DOWNLOAD_DIR}/${fn}" -C "${OUTDIR}" ;then
	    rm -rf ${DOWNLOAD_DIR}/${fn}
	    AC_MSG_NOTICE([File downloaded is not correct.])
	    return 1
	else
	    return 0
	fi
    fi
}
download(){
    fn="$2"
    OUTDIR="$3"
	if test "$enable_binary" = "no" -a -z "$4" ;then
		echo "Download pre-built binary is disabled"
		return
	fi
    mkdir -p ${DOWNLOAD_DIR} || die "failed to mkdir"
    if test -f ${DOWNLOAD_DIR}/${fn} ;then
		echo "Reusing downloaded ${DOWNLOAD_DIR}/${fn}"
		extract "${fn}" "$OUTDIR"
    fi
    if ! test -f ${DOWNLOAD_DIR}/${fn} ;then
		echo $wget "${BASEURL}/${fn}" "$wgetoutput" "${DOWNLOAD_DIR}/${fn}" 
		$wget "${BASEURL}/${fn}" "$wgetoutput" "${DOWNLOAD_DIR}/${fn}"
		if ! extract "$fn" "$OUTDIR" ;then
			AC_MSG_NOTICE([Failed to download the file from ${BASEURL}/${fn}])
			rm -rf ${DOWNLOAD_DIR}/${fn}
		fi
    fi
}
#compile LIBNAME srcfile.tar.gz [compile options]
compile(){
    echo Download and compile $1
    name=$1
    file=$2
    shift
    shift
    rm -rf ${TMP_DIR}/$name
    download $name $file ${TMP_DIR}/$name 1
    fnlog=${TMP_DIR}/compile.log
    echo "Compiling $name with CFLAGS=$CFLAGS $@. Check $fnlog for details"
    echo CFLAGS=$CFLAGS > $fnlog
    if (cd ${TMP_DIR}/$name/* && export CFLAGS="$CFLAGS -O3 -Wno-error"  && \
		((./configure --libdir=${LIB_DIR} $@ --prefix=${PRE_DIR} && make -j4 install) ||
		 (mkdir build && cd build && cmake .. -DCMAKE_INSTALL_PREFIX:PATH=${PRE_DIR} -DCMAKE_INSTALL_LIBDIR=lib${libsuffix} -DLIB_SUFFIX=${libsuffix} $@ && make -j4 install)))	>> $fnlog 2>&1 ;then
		echo "Compiled $name successfully."
    else
		echo "Failed to compile $name. Please see $fnlog for details"
		cat $fnlog
    fi
}
function set_lapack(){
    type=$1
    case $type in
	mkl*)
	    if test "$type" = "mkl-download" -a "$mkl_avail" = "yes";then
			download MKL "mkl${libsuffix}.tar.bz2" ${PRE_DIR}
	    fi
	    #First use MKL if available. 

	    if test "$host_cpu" = "x86_64" ;then
			LAPACK="-lmkl_intel_lp64"
	    else
			LAPACK="-lmkl_intel"
	    fi

	    LAPACK+=" -lmkl_core"

	    if test -z "$OPENMP_CFLAGS"; then
			LAPACK+=" -lmkl_sequential -lpthread"
	    else 
			if test "$CN" = "ICC" -o "$CN" = "CLANG" ;then
				LAPACK+=" -lmkl_intel_thread" #OPENMP_CFLAGS was already in CFLAGS
			else #GCC or others. Notice that clang can work with either thread model.
				LAPACK+=" -lmkl_gnu_thread"
			fi
	    fi
	    if test "$type" = "mkl-static" ;then 
			LAPACK="-Wl,--start-group,${LAPACK},--end-group"
			LAPACK=${LAPACK// /,}
			LAPACK=${LAPACK//,,/,}
	    fi
	    ;;
	pkg-config)
	    LAPACK=`pkg-config --libs blas lapack 2>/dev/null`
	    ;;
	openblas-openmp)
	    LAPACK="-lopenblaso"
	    ;;
	openblas-threads)
	    LAPACK="-lopenblasp"
	    ;;
	openblas)
	    LAPACK="-lopenblas"
	    ;;
	default)
	    LAPACK="-llapack -lblas"
		;;
	default-3)
		LAPACK="-l:liblapack.so.3 -l:libblas.so.3 " #notice that this method is not supported in macos (clang).
		;;
	*)
	    LAPACK=
	    ;;
    esac
    echo trying LAPACK=$LAPACK
}
has_lapack=no
base_lapack="pkg-config default openblas-openmp openblas-threads openblas default-3"
case $use_mkl in
    yes) #prefer MKL
		trys="$mkl mkl mkl-download $base_lapack"
		;;
    no) #disable MKL
		trys="$base_lapack"
		;;
    default) #default order
		trys="$mkl $base_lapack mkl-download"
		;;
esac

for try in ${trys}; do
    if test "$has_lapack" = "no" -o "$has_blas" = "no";then
		unset ac_cv_lib_m_dgemm_
		unset ac_cv_lib_m_dposv_
		set_lapack ${try}
		if test -n "$LAPACK" ;then
			#the trailing -lm is necessary for ubuntu-16.04 with mkl
			AC_CHECK_LIB([m],[dposv_],[has_lapack=yes],[has_lapack=no],[$LAPACK -lm])
			AC_CHECK_LIB([m],[dgemm_],[has_blas=yes],[has_blas=no],[$LAPACK -lm])
		fi
    fi
done

if test "$has_lapack" = "no" -o "$has_blas" = "no" ;then
    AC_MSG_ERROR([Please install blas/lapack])
fi

#Check for zlib
AC_CHECK_LIB([z],[gzopen],[has_z='yes'],[has_z='no'])
if test "$has_z" = "no" ;then
    compile ZLIB zlib.tar.bz2 #--static
    unset ac_cv_lib_z_gzopen
    AC_CHECK_LIB([z],[gzopen],[has_z='yes'],[AC_MSG_ERROR([Invalid zlib])])
fi

#Check FFTW header and library. FFTW3 interface in MKL performs as good as FFTW3
AC_CHECK_HEADERS([fftw3.h],[
	AC_CHECK_LIB([fftw3], [fftw_execute], [has_fftw="yes" LDFFT="-lfftw3"],[
		AC_CHECK_LIB([m], [fftw_execute], [has_fftw="builtin"], [has_fftw="no"], [$LAPACK -lm])
	])
], [has_fftw='no'])
if test "$has_fftw" = "builtin" -a "$with_matlab" != "no" ;then
	has_fftw=no #requires FFT for matlab
fi
if test "$has_fftw" = "builtin";then
    AC_MSG_NOTICE([FFTW is already provided by $LAPACK])
	AC_CHECK_LIB([m], [fftwf_execute],[has_fftwf=1] ,[has_fftwf=0])
    fftw_threads=1
else
    if test "$has_fftw" = "no" ;then #fftw does not exist in system.
		#download pre-compiled
		download FFTW "fftw${libsuffix}.tar.bz2" ${PRE_DIR}
		unset ac_cv_lib_fftw3_fftw_execute
		AC_CHECK_LIB([fftw3], [fftw_execute],[has_fftw='yes' LDFFT="-lfftw3"] ,[has_fftw='no'])

		if test "$has_fftw" = "no" ;then #fftw does not exist in system.
			#compile from source
			FLAGS="--enable-threads --enable-static --enable-shared"
			if test -n "$OPENMP_CFLAGS" ;then
				FLAGS+=" --enable-openmp"
			fi
			FLAGS+=" --disable-fortran --enable-sse2 --enable-avx --enable-avx2"
			compile FFTW fftw.tar.bz2 $FLAGS  
			compile FFTW fftw.tar.bz2 $FLAGS --enable-float 
			unset ac_cv_lib_fftw3_fftw_execute
			AC_CHECK_LIB([fftw3], [fftw_execute],[has_fftw='yes' LDFFT="-lfftw3"] ,[has_fftw='no'])
			AS_IF([test "$has_fftw" = "no"], AC_MSG_ERROR(["FFTW library or header is not available."]))
		fi
	fi
	AC_CHECK_LIB([fftw3f], [fftwf_execute],[has_fftwf=1 LDFFT+=" -lfftw3f"] ,[has_fftwf=0])
	
	if test -z "$OPENMP_CFLAGS" ;then
		THREADS_SUFFIX="_threads"
	else
		THREADS_SUFFIX="_omp"
	fi
	LDFFT_THREADS=
	for name in $LDFFT ;do
		case "$name" in
		*.a)
			LDFFT_THREADS+=" ${name/.a/}${THREADS_SUFFIX}.a"
			;;
		*.so)
			LDFFT_THREADS+=" ${name/.so/}${THREADS_SUFFIX}.so"
			;;
		*)
			LDFFT_THREADS+=" $name${THREADS_SUFFIX}"
		esac
	done
	AC_CHECK_LIB([m], [fftw_plan_with_nthreads], [fftw_threads=1 LDFFT+=" $LDFFT_THREADS"], [fftw_threads=0], [$LDFFT_THREADS $LDFFT $OPENMP_CFLAGS])
	LDFFT_MEX=$LDFFT 
fi
AC_DEFINE_UNQUOTED(HAS_FFTW_THREADS, [$fftw_threads], "Has FFTW3 threads support.")
AC_DEFINE_UNQUOTED(HAS_FFTWF, [$has_fftwf], "FFTW single precision exists")
AS_IF([test x$fftwf = 0 -a x$use_double = xno], [AC_MSG_ERROR(["FFT single precision is not available for --disable-double"])])
#Check for cholmod
AC_CHECK_LIB([cholmod], [${chol_fun}], [has_chol="yes"], [has_chol="no"], [$LAPACK -lm])
AC_CHECK_HEADERS([cholmod.h], [], [AC_CHECK_HEADERS([suitesparse/cholmod.h], [], [has_chol="no"])])
if test "$has_chol" = "no" ;then
	download CHOLMOD "cholmod${libsuffix}.tar.bz2" ${PRE_DIR}
	unset ac_cv_lib_cholmod_${chol_fun}
	AC_CHECK_LIB([cholmod], [${chol_fun}], [has_chol="yes"], [has_chol="no"], [$LAPACK -lm])

	if test "$has_chol" = "no" ;then
		compile CHOLMOD "cholmod.tar.bz2"
		unset ac_cv_lib_cholmod_${chol_fun}
		AC_CHECK_LIB([cholmod], [${chol_fun}], [has_chol="yes"], [AC_MSG_ERROR([Suitesparse library is not available])], [$LAPACK -lm])
	fi
	unset ac_cv_header_cholmod_h ac_cv_header_suitesparse_cholmod_h
	AC_CHECK_HEADERS([cholmod.h], [], [AC_CHECK_HEADERS([suitesparse/cholmod.h], [], [AC_MSG_ERROR([Suitesparse header is not available])])])
fi
LDCHOL="-lcholmod"

#Check for libwebsocket
AM_CONDITIONAL(HAS_LWS, false)
if test "$enable_websockets" != "no" ;then
    AC_CHECK_LIB([websockets], [lws_service], [has_lws="yes"], [has_lws="no"], [])
    AC_CHECK_HEADERS([libwebsockets.h], [], [has_lws="no"])
    if test "$has_lws" = "no" ;then #then download pre-compiled library and header.
		download LIBWEBSOCKETS "libwebsockets${libsuffix}.tar.bz2" ${PRE_DIR}
		unset ac_cv_lib_websockets_lws_service
		AC_CHECK_LIB([websockets], [lws_service], [has_lws="yes"], [has_lws="no"], [])
    
		if test "$has_lws" = "no" ;then #try to compile from source. Needs cmake
			compile LIBWEBSOCKETS libwebsockets.tar.bz2 -DLWS_WITH_SSL=OFF -DLWS_IPV6=OFF -DLWS_WITHOUT_TESTAPPS=ON -DLWS_WITHOUT_DAEMONIZE=ON -DLWS_WITHOUT_CLIENT=ON -DLWS_WITH_SHARED=OFF -DLWS_WITH_STATIC=ON
			unset ac_cv_lib_websockets_lws_service
			AC_CHECK_LIB([websockets], [lws_service], [has_lws="yes"], [has_lws="no"], [])
		fi
		
		unset ac_cv_header_libwebsockets_h
		AC_CHECK_HEADERS([libwebsockets.h], [], [AC_MSG_ERROR([libwebsocket header is not available])])
	fi
	if test "$has_lws" != "no" ;then
		AM_CONDITIONAL(HAS_LWS, true)
		AC_DEFINE(HAS_LWS, [1], "libwebsocket is available")
	else
		AC_MSG_NOTICE(["libwebsocket library or header is not available"])
		enable_websockets="no"
	fi
fi
if test "$enable_cmocka" != "no" ;then
	AC_CHECK_LIB([cmocka], [mock_assert], [enable_cmocka="yes"], [enable_cmocka="no"])
	AC_CHECK_HEADERS([cmocka.h],[],[enable_cmocka="no"],[
		#include <stdarg.h>
		#include <stddef.h>
		#include <stdint.h>
		#include <setjmp.h>])
	if test "$enable_cmocka" = "no" ;then
		download CMOCKA cmocka${libsuffix}.tar.bz2 ${PRE_DIR}
		unset ac_cv_lib_cmocka_mock_assert
		AC_CHECK_LIB([cmocka], [mock_assert], [enable_cmocka="yes"], [enable_cmocka="no"])
	
		if test "$enable_cmocka" = "no" ;then
			compile CMOCKA cmocka.tar.xz
			unset ac_cv_lib_cmocka_mock_assert
			AC_CHECK_LIB([cmocka], [mock_assert], [enable_cmocka="yes"], [enable_cmocka="no"])
		fi
	fi
	if test "$enable_cmocka" != "no" ;then
		AC_DEFINE(USE_CMOCKA, [1], "Enable unit testing with cmocka")
	fi
fi
AM_CONDITIONAL(USE_CMOCKA, [ test x$enable_cmocka != xno ])
PKG_PROG_PKG_CONFIG
gtk_ver="0"
AS_IF([ test -n "$PKG_CONFIG"],[#use AS_IF instead of if to avoid error regarding pkgconfig not available.
    #Check availability of gtk and libnotify. If available, we will build drawdaemon and the monitor.
    #cairo 1.2 comes with gtk+-2.10
    #make sure there are spaces before and after >=
	if test -n "$GTKDIR" ;then #User specified GTK directory
		if test -d "$GTKDIR/pkgconfig" ;then
			GTKDIR="$GTKDIR/../"
		fi
		if test -d "$GTKDIR/lib/pkgconfig" ;then
			export PKG_CONFIG_PATH="$PKG_CONFIG_PATH:$GTKDIR/lib/pkgconfig"
		fi
		TT=-I$GTKDIR/include
		GTK_CFLAGS="$TT $TT/cairo $TT/atk-1.0 $TT/gtk-${gtk_ver}.0 $TT/glib-2.0 $TT/pango-1.0 $TT/gdk-pixbuf-2.0 -I$GTKDIR/lib/glib-2.0/include -I${GTKDIR}/lib/gtk-2.0/include"
		GTK_LIBS="-L$GTKDIR/lib/ -Wl,-rpath,$GTKDIR/lib/"
		if test -d "$GTKDIR/bin" ;then
			GTK_LIBS="$GTK_LIBS -L$GTKDIR/bin -Wl,-rpath,GTKDIR/bin"
		fi
		gtk_ver="2"
    else
		AS_IF([test "$use_gtk_4" != "no" ], [PKG_CHECK_MODULES([GTK], gtk4,               [gtk_ver="4"], [gtk_ver="0"])])
		AS_IF([test "$gtk_ver" = "0" -a "$use_gtk_3" != "no" ], [PKG_CHECK_MODULES([GTK], gtk+-3.0,           [gtk_ver="3"], [gtk_ver="0"])])
		AS_IF([test "$gtk_ver" = "0" ],     [PKG_CHECK_MODULES([GTK], gtk+-2.0 >= 2.2,    [gtk_ver="2"], [gtk_ver="0"])])
    fi
    AS_IF([test "$gtk_ver" != "0"], [
	      PKG_CHECK_MODULES(DRAW,   cairo >= 1.4,    [has_cairo="1"], [has_cairo="0"])
	      PKG_CHECK_MODULES(NOTIFY, libnotify >= 0.1,[has_notify="1"], [has_notify="0"])
	      PKG_CHECK_MODULES(MAC,    gtk-mac-integration-gtk${gtk_ver}, [has_mac="1"], [has_mac="0"])
	      ])
 
    if test "$gtk_ver" = "2" ;then
		if ! echo "$GTK_LIBS" | grep -q gthread ;then
			#OLD GLIB has GTHREAD separated
			PKG_CHECK_MODULES(GTHREAD, gthread-2.0 >= 2.0, [has_gthread="yes"], [has_gthread="no"])
			if test -n "$GTHREAD_LIBS" ;then
			GTK_LIBS+=" $GTHREAD_LIBS"
			fi
		fi
		GTK_CFLAGS+=" -Wno-deprecated-declarations"
    fi
    AC_DEFINE_UNQUOTED(MAC_INTEGRATION, [$has_mac], ["Mac integration"])
    #The following will substitute the variables name 
    #Reexport after updating
    AC_SUBST(GTK_CFLAGS)
    AC_SUBST(GTK_LIBS)
])

AM_CONDITIONAL(GUI,  test x$gtk_ver  != x0)
AM_CONDITIONAL(CAIRO, test x$has_cairo = x1)
AM_CONDITIONAL(NOTIFY, test x$has_notify = x1)

AC_CHECK_HEADERS([unistd.h])
OPTS+=" $OPTScommon $OPTScc $OPTSprof "

if test "$with_matlab" != "no" ;then
    if test ${MEXDIR:0:1} = '~' ;then
		MEXDIR=$HOME:${MEXDIR:1}
    fi
    if test -x "${MEXDIR}/mex" ;then
		MEXDIR="${MEXDIR}"
    elif test -x "${MEXDIR}/bin/mex" ;then
		MEXDIR="${MEXDIR}/bin"
    else
		AC_MSG_ERROR([mex not found at $MEXDIR. please reconfigure with --with-matlab=DIR to specify MATLAB location or --without-matlab to disable])
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
	    #MEXOPT="${MEXOPT} -DDLONG -D_BSD_SOURCE -D_POSIX_C_SOURCE=200809L"
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
	    MEXOPT="${MEXOPT} -DDLONG"
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
	AC_CHECK_LIB([mx], [mxGetPr], [mex_lib="yes"], [mex_lib="no"], [$LDMEX])
    fi
    if test x$mex_lib = xyes ;then
    	MEXOPT="${MEXOPT} -I${MEXDIR}/../extern/include -DMATLAB_MEX_FILE $OPENMP_CFLAGS" 
		#2021-05-24: Linking with mwfftw3 crashes matlab. (mwfftw3 contains thread support)
		LIBMEX=" -lmwlapack -lmwblas -lmx -lmex -lmat -lut -lstdc++ -lz -lm -ldl $LDRT $LDFFT_MEX"
		LDMEX+="-no-fast-install -static -Xcompiler -shared"
		AM_CONDITIONAL(USE_MEX, true)
    else
		with_matlab="no"
		AC_MSG_NOTICE([mex library test failed. Disable mex])
		AM_CONDITIONAL(USE_MEX, false)
    fi
    AC_SUBST(MEXEXT)
    AC_SUBST(LIBMEX)
    AC_SUBST(LDMEX)
    AC_SUBST(MEXOPT)
else
    AM_CONDITIONAL(USE_MEX, false)
fi

if test "$with_cuda" != "no" ;then
    if test -n "$with_cuda" ;then
		user_enable_cuda=yes
    fi
    #if test "$enable_allstatic" = "yes"; then
	#AC_MSG_ERROR("all-static conflicts with cuda support. Please disable all-static")
    #fi
    AC_MSG_NOTICE([Enabling CUDA])
    #Check for nvcc executable
    if test ${CUDADIR:0:1} = '~' ;then
		CUDADIR=$HOME:${CUDADIR:1}
    fi
    NVCC="$CUDADIR/bin/nvcc"
    if test ! -x "$NVCC" ;then
		NVCC="$CUDADIR/nvcc"
		CUDADIR="$(dirname $CUDADIR)"
		if test ! -x "$NVCC" ;then
            AC_MSG_NOTICE([nvcc not found])
	    	with_cuda="no"
		fi
    fi
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
		LIBCUDA="-L$CUDA_L -Wl,-rpath,$CUDA_L"
    else
		LIBCUDA=
		AC_MSG_NOTICE([CUDA library not found])
		with_cuda="no"
    fi
    AC_CHECK_LIB([cudart], [cudaSetDevice], [], [AC_MSG_NOTICE([libcudart not found]);with_cuda="no"], [$LIBCUDA])

    #Check for cuda header
    if test -f "${CUDADIR}/include/cuda.h" ;then
		CUDAOPT="-I${CUDADIR}/include"
    else
		AC_CHECK_HEADERS([cuda.h], [] ,[AC_MSG_NOTICE([Header not found]);with_cuda="no"])
		CUDAOPT=
    fi
    if test "$with_cuda" = "no" ;then
		if test "$user_enable_cuda" = "yes" ;then
			AC_MSG_ERROR([Specified cuda DIR is not valid])
		else
			AC_MSG_NOTICE(["Cannot compile with cuda. Please use --with-cuda=DIR to specify CUDA location or --without-cuda"])
		fi
    fi
fi

if test "$with_cuda" != "no" ;then
    AC_MSG_NOTICE(["Use cuda in directory $CUDADIR"])
    cudaver="$($NVCC -V |grep release |cut -d ',' -f 2 |cut -d ' ' -f 3)"
    cudaver="${cudaver/.}"
    echo CUDA toolkit is $cudaver
    NVCC+=" $CFLAGS_CPU"
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
    if test "$cudaver" -lt 90 ;then
		cuarch=20 #<cuda-8.0 supports compute_20
    elif test "$cudaver" -lt 110 ;then
		cuarch=30 #compute_30 is the minimum version
	else
		cuarch=50 #compute_50 is the minimum version
    fi
	CUDAOPT+=" -arch=compute_$cuarch" 
    #To specify multiple architectures, use
    #CUDAOPT+=" -gencode arch=compute_20 -gencode arch=compute_30"
    if test "$enable_debug" = "yes" ;then
        CUDAOPT+=" -O0"
    else
	CUDAOPT+=" -O3"
    fi
    CUDAOPT+=" -g -lineinfo -DHAVE_CONFIG_H -I.. "
    #CUDAOPT+="-Xptxas=-v -dlcm=cg"
    if test "$host_cpu" = "x86_64" ;then
        CUDAOPT+=" -DDLONG"
    fi
    #CUDAOPT+=" -std=c++11"
    CUDAOPT+=" -Wno-deprecated-gpu-targets"
    CCBINOPT="$OPTScommon $OPTSprof -fPIC $CPPFLAGS $CFLAGS -Wno-unused-parameter"
    CCBINOPT="${CCBINOPT/-fpermissive/}"
    LIBCUDA+=" -lcurand -lcusparse -lcufft -lcublas -lcudart"
    if $CCBIN -v 2>&1 |grep clang ;then
		CCN=CLANG
		LIBCUDA+=" -lstdc++"
    elif $CCBIN -v 2>&1 |grep icc ;then
	 	CCN=ICC
	 	LIBCUDA+=" -lstdc++"
    else
		CCN=GCC
		CCBINOPT+=" -Wno-unused-value -D__STRICT_ANSI__ $"
		LIBCUDA+=" -lstdc++"
    fi
    optname=${ccprefix}_${CCN}
    CCBINOPT+=" ${!optname}"
    if test -n "$OPENMP_CFLAGS";then
        CCBINOPT+=" $OPENMP_CFLAGS"
    fi
else
    cudaver="no"
fi
if test "$with_cuda" != "no" ;then
    if test "$cuda_double" = "yes" ;then
		AC_DEFINE(CUDA_DOUBLE, [1], "Use double precision mode in cuda")
		AM_CONDITIONAL(CUDA_DOUBLE, true)
    else
		AC_DEFINE(CUDA_DOUBLE, [0], "Use single precision mode in cuda")
		AM_CONDITIONAL(CUDA_DOUBLE, false)
    fi
    AM_CONDITIONAL(USE_CUDA, true)
    AC_DEFINE(USE_CUDA, [1], "Has Cuda")
    AC_SUBST(CUDAOPT)
    AC_SUBST(CCBINOPT)
    AC_SUBST(LIBCUDA)
    AC_SUBST(NVCC)
else
    AM_CONDITIONAL(USE_CUDA, false)
    AM_CONDITIONAL(CUDA_DOUBLE, false)
    AC_DEFINE(USE_CUDA, [0], "Don't have Cuda")
fi

#Doxygen Features
if test "$enable_doxygen" = "yes" ;then
    if ! which doxygen ;then
        download DOXYGEN doxygen${libsuffix}.tar.bz2 ${PRE_DIR}
    fi
    PATH+=":${HOME}/bin"
    if ! which latex ;then
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
	tools/Makefile
	maos/Makefile 
	skyc/Makefile
	bin/Makefile
	mex/Makefile
	cuda/Makefile
	])
if test -d "${srcdir}/test" ;then
    AC_CONFIG_FILES([test/Makefile])
fi

LDOPTS+=" -no-fast-install -avoid-version"
LDEXE+=" -lz -lm -ldl"
LDEXE_FULL="$LDFFT $LDCHOL $LAPACK $LDEXE" #when require libaos
echo LIBS=$LIBS
LIBS="" #somehow LIBS contain -fftw3 -lnvidia-smlA
AC_SUBST(ldsuffix)
AC_SUBST(CN)
AC_SUBST(LDAOS)
AC_SUBST(OPTS)
AC_SUBST(LDCHOL)
AC_SUBST(LDFFT)
AC_SUBST(LAPACK)
AC_SUBST(LDEXE)
AC_SUBST(LDEXE_FULL)
AC_SUBST(LDEXE_MAOS)
AC_SUBST(SRC_DIR)
AC_SUBST(BUILD_DIR)
AC_SUBST(LIBS)
AC_SUBST(LDOPTS)
AC_SUBST(LDSO)
AC_OUTPUT
echo SRC_DIR=$SRC_DIR
echo BUILD_DIR=$BUILD_DIR
echo host_os=$host_os
echo BIGENDIAN=$BIGENDIAN
echo enable_allstatic=$enable_allstatic

echo CN=$CN
echo CFLAGS=$CFLAGS
echo CPPFLAGS=$CPPFLAGS
echo LDFLAGS=$LDFLAGS
echo OPTS=$OPTS
echo LDOPTS=$LDOPTS
echo LDEXE_FULL=$LDEXE_FULL
if test "$gtk_ver" != "0" ;then
    echo DRAW_LIBS=$DRAW_LIBS
    echo GTK_LIBS=$GTK_LIBS
fi
if test "$with_cuda" != "no" ;then
    echo NVCC=$NVCC
    echo CUDAOPT=$CUDAOPT
    echo LIBCUDA=$LIBCUDA
fi
if test "$with_matlab" != "no" ;then
    echo MEXDIR=$MEXDIR
    echo LIBMEX=$LIBMEX
    echo MEXOPT=$MEXOPT
fi
echo
echo Compiling options:
echo 
echo "CC:              $CC"
echo "OpenMP:          $has_openmp"
echo "Single:          $enable_single"
echo "CUDA:            $cudaver (double: $cuda_double, arch: $cuarch)"
echo "Drawdaemon:      $has_cairo"
echo "Monitor:         $((gtk_ver>0)) (libnotify: $has_notify)"
echo "GTK version:     ${gtk_ver}"

echo "Websocket:       ${enable_websockets:-yes}"
echo "MATLAB MEX:      $with_matlab"

echo "Doxygen:         ${DX_DOXYGEN:+yes (make doc)}"
echo "CMOCKA:          $enable_cmocka"
echo

echo "*************************************************************************"
echo "*** MAOS Documentation can be found at http://lianqiw.github.io/maos ***"
echo "*** MAOS Source can be obtained at http://github.com/lianqiw/maos     ***"
echo "*************************************************************************"

AC_PREREQ(2.57)
AC_INIT(maos, 1.8.1, lianqiw@gmail.com)
AM_INIT_AUTOMAKE([1.9.6 gnu -Wall no-define subdir-objects])
#silent building
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])]) 
#The URL where misc supportive libaries are located.
BASEURLS="https://sourceforge.net/projects/lianqiw/files/external"
BASEURLS+=" http://project.tmt.org/~lianqiw/maos/external"
#If user does not supply any CFLAGS, set something so that AC_PROG_CC won't put a -g -O2 to CFLAGS
if test -z "$CFLAGS" ;then
    CFLAGS=" "
fi
if test -z "$LIBS" ;then
    LIBS=" "
fi
AC_CONFIG_MACRO_DIR([m4])
#After enable HEADER, the package name will not appear in CFLAGS. Appear in config.h instead.
AC_CONFIG_HEADER([config.h])
#defines host host_os host_cpu host_vendor etc
AC_CANONICAL_HOST 
AC_LANG([C])
AC_TYPE_SIZE_T #it sets CC if empty
AC_PROG_CC([$CC])
AM_PROG_CC_C_O
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
#Doxygen Features
DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)
DX_INIT_DOXYGEN($PACKAGE_NAME, doxygen.cfg)
if test x$DX_FLAG_doc = x0 ;then
    AM_CONDITIONAL(DO_DOXYGEN, false)
else
    AM_CONDITIONAL(DO_DOXYGEN, true)
fi

AC_C_BIGENDIAN([BIGENDIAN=1], [BIGENDIAN=0], AC_MSG_ERROR("Unknown Endianness"))
AC_DEFINE_UNQUOTED(BIGENDIAN, [$BIGENDIAN], ["Endianness"])

#top_srcdir and abs_top_srcdir is empty. no idea why. but srcdir sometimes is the full path.
SRC_DIR=$(cd "$srcdir" && pwd)
BUILD_DIR=$(pwd)
#AC_DEFINE_UNQUOTED means expand the $ variable
AC_DEFINE_UNQUOTED(SRCDIR,  ["$SRC_DIR"], ["Define SRC folder"])
AC_DEFINE_UNQUOTED(BUILDDIR,["$BUILD_DIR"],["Define BUILD folder"])
#Specifying both -O3 and -g in mac produces internal error in gcc 4.2
#-all-static will produce a static exe if all libraries are available in static linking mode.
OPTScommon="-Wno-unknown-pragmas -Wno-unused-function -fno-common" 
if test "${USER}" = "lianqiw" -a "$CC" = "gcc" ;then
    OPTScommon="-Wshadow -Wall -Werror $OPTScommon"
fi
LDEXE="-no-fast-install -ldl"
if test -z "$OPENMP_CFLAGS" ;then
    LDEXE+=" -lpthread"
fi
OPTSextra=""

#Enable the user to force 32 or 64 bit program
AC_ARG_ENABLE(64,  AS_HELP_STRING([--enable-64],[Force 64 (disable=32) bit executables]))

if test "$enable_64" = "yes" ;then
    host_cpu=x86_64
elif test "$enable_64" = "no" ;then
    host_cpu=i386
else
    host_cpu=$(uname -m)
fi
case ${host_cpu} in 
    *64)
	host_cpu=x86_64
	CFLAGS_CPU="-m64"
	CFLAGS_ARCH="-m64 -DDLONG"
        LDFLAGS+=" -m64"
        CDEFINE="__x86_64__"
	AC_DEFINE(X86_64, [1], "X86_64 Machine")
	objcopy_output=elf64-x86-64
	objcopy_arch=i386:x86-64
	echo "CPU is 64 bit"
	;;
    *86)
	host_cpu=i386
	CFLAGS_CPU="-m32"
	CFLAGS_ARCH="-m32 -DDINT -D_FILE_OFFSET_BITS=64"
	LDFLAGS+=" -m32"
        CDEFINE="__i386__"
	AC_DEFINE(X86, [1], "I386 Machine")
	objcopy_output=elf32-i386
	objcopy_arch=i386
	echo "CPU is 32 bit"
	;;
    *)
        echo host_cpu=$host_cpu
        AC_MSG_ERROR([Unrecognized host_cpu])
	;;
esac   
#Used for embedding conf folder in maos exe
AC_SUBST(objcopy_output)
AC_SUBST(objcopy_arch)

#Check operating system.
AM_CONDITIONAL(LINUX,  false)
AM_CONDITIONAL(BSD,    false)
AM_CONDITIONAL(DARWIN, false)
AM_CONDITIONAL(WIN,    false)
LDLIB="-avoid-version"
mkl_avail="no"
case "$host_os" in
    *linux*)
	AM_CONDITIONAL(LINUX, true)
	LDEXE+=" -lrt" #need real time library.
	CFLAGS_EXTRA="-g -D_BSD_SOURCE -D_POSIX_C_SOURCE=200112L -fPIC"
	CFLAGS_ARCH+=" -D_GLIBC_2_5_SOURCE"
	if test "${host_cpu}" = "x86_64" ;then
	    libsuffix="64"
	else
	    libsuffix="32"
	fi
	mkl_avail="yes"
	system="linux"
	ldsuffix="so"
	;;
    *darwin*)
	AM_CONDITIONAL(DARWIN, true)
	CFLAGS_EXTRA="-fPIC"
	mkl_avail="no"
	libsuffix=fat
	system="apple"
	ldsuffix="dylib"
	#for compatibility with OS X >=10.6. Don't use 10.5 which will render mex files not usable.
	CFLAGS_ARCH+=" -mmacosx-version-min=10.6"
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
AC_SUBST(ldsuffix)
PRE_DIR="${HOME}/.aos/"          #prefix for installation
TMP_DIR="${BUILD_DIR}/external/" #temporary directory for compiling
LIB_DIR="${PRE_DIR}/lib${libsuffix}" #external library
LDFLAGS+=" -L${LIB_DIR} -Wl,-rpath,${LIB_DIR}" #-avoid-version
DOWNLOAD_DIR="${TMP_DIR}/downloads"
CFLAGS="$CFLAGS_EXTRA $CFLAGS_ARCH -I${PRE_DIR}/include"

die(){
    shift
    echo "Error:$@"
    exit
}
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
    for nvcc in "/opt/cuda/default/bin/nvcc" "/opt/cuda/bin/nvcc" "/usr/local/cuda/" ;do
	if test -x "$nvcc" ;then
	    NVCC="$nvcc"
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
if test -n "`which mex 2>/dev/null`" ;then
    MEX="$(which mex)"
    if test -L "$MEX" ;then
	MEX="$(readlink "$MEX")"
    fi
    MEXDIR_AUTO="$(dirname "$MEX")"
elif test -n "`ls -d /Applications/MATLAB*.app 2>/dev/null`" ;then
    for i in "$(ls -d /Applications/MATLAB*.app)" ;do
	if test -x "$i/bin/mex" ;then
       	    MEXDIR_AUTO="${i}/bin"
	fi
    done
fi

#Discover default download method
if test -x "`which wget 2>/dev/null`" ;then
    wget="wget -t1 --no-check-certificate"
    wgetoutput="-O"
elif test -x "`which curl 2>/dev/null`" ;then
    wget="curl --retry 1 -L --progress-bar"
    wgetoutput="--output"
else
    wget="die wget/curl not found ; "
fi
#Common configure options.
AC_ARG_ENABLE(all-static, AS_HELP_STRING([--diable-all-static], [Compile statically]), [use_static="$enableval"],[use_static="no"])
AC_ARG_ENABLE(debug, AS_HELP_STRING([--enable-debug],[Turn on debug]), [debug="$enableval"], [debug="no"])
AC_ARG_ENABLE(oprof, AS_HELP_STRING([--enable-oprof],[Enable oprof for code profiling.]),[use_oprof="$enableval"],[use_oprf="no"])
AC_ARG_ENABLE(cuda-double, AS_HELP_STRING([--enable-cuda-double],[Enable cuda double precision]), [cuda_double="$enableval"],[cuda_double="no"])
AC_ARG_ENABLE(scheduler,AS_HELP_STRING([--disable-scheduler],["Disable built in scheduler"]), [use_scheduler="$enableval"], [use_scheduler="yes"])
AC_ARG_ENABLE(websockets,AS_HELP_STRING([--disable-websockets],[Disable websockets (default: enabled)]), [use_websockets="$enableval"], [use_websockets="yes"])
AC_ARG_WITH(cuda, AS_HELP_STRING([--with-cuda[=DIR]],[With cuda support]), [CUDADIR=${withval}], [CUDADIR="$CUDADIR_AUTO"])
AC_ARG_WITH(ccbin, AS_HELP_STRING([--with-ccbin=CC,],[Specify CCBIN for nvcc]), [CCBIN=${withval}], [CCBIN=${CC}])
AC_ARG_WITH(gtk,  AS_HELP_STRING([--with-gtk[=DIR]],[Enable gtk support]), [GTKDIR=${withval}], [GTKDIR=])
AC_ARG_ENABLE(gtk-3, AS_HELP_STRING([--enable-gtk-3],[Using GTK+-3]), [use_gtk_3="$enableval"], [use_gtk_3="no"])
AC_ARG_WITH(matlab, AS_HELP_STRING([--with-matlab[=DIR]],[Enable compiling mex routines for matlab]), [MEXDIR=${withval}], [MEXDIR=${MEXDIR_AUTO}])
if test -z "$MEXDIR" ;then
    with_matlab=no
fi

AM_CONDITIONAL(USE_SCHEDULER, true)
if test "$use_scheduler" = "no" ;then
    AC_DEFINE(MAOS_DISABLE_SCHEDULER, [1], "Disable built in scheduler")
    AM_CONDITIONAL(USE_SCHEDULER, false)
    use_websockets="no"
fi
if test "$host_cpu" = "x86_64" ;then 
    chol_fun=cholmod_l_factorize
    AM_CONDITIONAL(USE_LONG, true)
else
    chol_fun=cholmod_factorize
    AM_CONDITIONAL(USE_LONG, false)
fi

if test "x$use_static" = "xyes" ;then
    LDEXE+=" -all-static" #Enable all static build.
    AM_CONDITIONAL(USE_STATIC, true)
    AC_DEFINE(USE_STATIC, [1], "We are compling maos statically")
    #(cd ${SRC_DIR} && tar zcf ${BUILD_DIR}/config.tar.gz config/)
else
    LDEXE+=" -static" #show that we prefer static linking.
    AM_CONDITIONAL(USE_STATIC, false)
    AC_DEFINE(USE_STATIC, [0], "We are compling maos dynamically")
fi

if test x$use_oprof = "xyes" ;then
    OPTSprof+=" -g -fno-omit-frame-pointer -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls"
fi
#Set optimize flags for different compilers
#-xO Enables SSE3, SSE2 and SSE instruction sets optimizations for non-Intel CPUs 
COMMON_ICC="-wd193 -wd1418 -wd2259 -wd1372 -wd981 -wd1292 -wd10237 -wd858 -wd1478 -wd1786 -wd3180 -diag-disable remark -vec-report=3"
COMMON_GCC="-fms-extensions -static-libgcc"
COMMON_GCC+=" -Wno-deprecated-declarations" 
COMMON_CLANG="-Wno-deprecated-declarations" #"-Qunused-arguments"
OPTIM_ALL="-O3 "
OPTIM_ICC="$OPTIM_ALL -funroll-loops " #-static-intel"
OPTIM_GCC="$OPTIM_ALL -funroll-loops -fprefetch-loop-arrays -fno-math-errno"
OPTIM_CLANG="$OPTIM_ALL"
DEBUG_ALL="-g -O0 -DDEBUG_TRACE"
DEBUG_ICC="$DEBUG_ALL"
DEBUG_GCC="$DEBUG_ALL -ggdb"
DEBUG_CLANG="$DEBUG_ALL"

#inline definition:

#in GCC, extern inline will use the inline and not emit function definition.
#in GCC, inline will cause use the inline, but also emit function deminition, causing multiple definitions if put in header.
#in C99, extern inline will cause use the inline, but also emit function deminition, causing multiple definitions if put in header. This is equivalent to GCC inline
#in C99, inline will use the inline and not emit function definition. This is equivalent to GCC extern inline 
#GCC mode can be forced on using __attribute__((gnu_inline))
#always_inline causes inline to be done in debug mode.
CN=""
cc_full=`which ${CC%% *}`
while test -L $cc_full ;do
    cc_full=`readlink $cc_full`
done
case "${cc_full}" in	
    *icc*)
	CN=icc
	COMMON="$COMMON_ICC"
	OPTIM="$OPTIM_ICC"
	DEBUG="$DEBUG_ICC"
	#Test ICC version
	CC_VERSION=$($CC -dumpversion)
	if test -n "$CC_VERSION" ;then
	    major=${CC_VERSION%%.*}
	    minor=${CC_VERSION%.*}
	    minor="${minor#4.}"
	    #if test $major -ge 13 ;then
	#	OPTIM+=" -axCORE-AVX2,SSE4.2"
	 #   elif test $major -ge 11 ;then
	#	OPTIM+=" -axSSE4.2"
	 #   fi
	    if test $major -lt 11 ;then
		unset OPENMP_CFLAGS
	    fi
	fi
	;;
    *gcc*|*g++*)
	if test -n "`$CC -v 2>&1 |grep clang`" ;then
            CN=clang
            OPTIM="$OPTIM_CLANG"
	    COMMON="$COMMON_CLANG"
	    DEBUG="$DEBUG_CLANG"
	else
            CN=gcc
	    OPTIM="$OPTIM_GCC"
	    COMMON="$COMMON_GCC"
	    DEBUG="$DEBUG_GCC"
	fi
	#Test GCC version
	CC_VERSION=$($CC -dumpversion)
	if test -n "$CC_VERSION" ;then
	    major="${CC_VERSION%%.*}"
	    minor="${CC_VERSION%.*}"
	    minor="${minor#4.}"
	    echo major=$major minor=$minor
	    if test "$major" -lt 4 ;then
		AC_MSG_ERROR([At least version 4 is required when compiling with gcc.])
	    fi
	    if test "$major" -ge 4 -a "$minor" -ge 2;then
		#echo "Enabling -mtune=native"
		#OPTIM="$OPTIM -mtune=native" #introduced in GCC 4.2
		OPTIM+=" -ftree-vectorize " #introduced in GCC 4
	    fi
	    #GCC4.3 changed meaning of inline to C99
	    if test $major -eq 4 -a $minor -lt 5;then
		gnu_inline="yes"
	    fi
	    #GCC introduced openmp task in 4.4
	    if test $major -lt 4 -o $minor -lt 4 ;then
		unset OPENMP_CFLAGS
	    fi
	fi
	;;
    *clang*)
        CN=clang
	OPTIM="$OPTIM_CLANG"
	COMMON="$COMMON_CLANG"
	DEBUG="$DEBUG_CLANG"
	;;
    *) #Other compilers
     	OPTIM="-O3"
        COMMON=""
	DEBUG="-O0"
	;;	
esac
AC_SUBST(CN)
if test "$gnu_inline" = "yes" ;then
    INLINE="extern __inline __attribute__((gnu_inline, always_inline))" #C99
else
    INLINE="inline __attribute__((always_inline))" #C99
fi
AC_DEFINE_UNQUOTED(INLINE, $INLINE, [Portable definition of inline])
if test "x$debug" = "xyes" ;then
    AM_CONDITIONAL(DEBUG,true)
    AC_DEFINE(DEBUG,[1],["Debug the code"])
    OPTScc="$DEBUG $COMMON"
else
    AM_CONDITIONAL(DEBUG,false)
    AC_DEFINE(DEBUG,[0],["Do not debug"])
    OPTScc="$OPTIM $COMMON"
fi


#call convention of fortran libraries
#AC_F77_FUNC(dgemm)
#AC_F77_FUNC(dposv)
#download LIBNANE file.tar.bz2 [outdir]
download(){
    fn="$2"
    if test -n "$3" ;then
	OUTDIR="$3"
    else
	AC_MSG_ERROR([Please specify third argument])
    fi
    mkdir -p ${DOWNLOAD_DIR} || die "failed to mkdir"
    retry=0
    while test ${retry} -le 2 ;do
	if test ! -f ${DOWNLOAD_DIR}/${fn} ;then
	    echo Downloading $1 from the MAOS website to ${DOWNLOAD_DIR}/${fn}
	    success=0
	    for BASEURL in $BASEURLS ; do
		echo $wget "${BASEURL}/${fn}" "$wgetoutput" "${DOWNLOAD_DIR}/${fn}" 
		if $wget "${BASEURL}/${fn}" "$wgetoutput" "${DOWNLOAD_DIR}/${fn}" ;then
		    break
		else
		    rm -rf ${DOWNLOAD_DIR}/${fn}
		    AC_MSG_NOTICE([Failed to download the file from ${BASEURL}/${fn}])
		fi
	    done
	    if test ! -f ${DOWNLOAD_DIR}/${fn} ;then
		AC_MSG_ERROR(["Please download it manually and put in ${DOWNLOAD_DIR}/${fn}"])
	    fi
	fi
	if test -f ${DOWNLOAD_DIR}/${fn} ;then
	    echo Extracting "${DOWNLOAD_DIR}/${fn}" to "${OUTDIR}"
	    mkdir -p "${OUTDIR}" || die "failed to mkdir"
	    if ! tar -xf "${DOWNLOAD_DIR}/${fn}" -C "${OUTDIR}" ;then
		rm -rf ${DOWNLOAD_DIR}/${fn}
		retry=$((retry+1))
		AC_MSG_NOTICE([File download is not correct, retry])
	    else
		break
	    fi
	fi
    done
}
#compile LIBNAME srcfile.tar.gz [compile options]
compile(){
    echo Download and compile $1
    name=$1
    file=$2
    shift
    shift
    rm -rf ${TMP_DIR}/$name
    download $name $file ${TMP_DIR}/$name
    fnlog=${TMP_DIR}/compile.log
    echo "Compiling $name. Check $fnlog for details"
    if (cd ${TMP_DIR}/$name/* && export CFLAGS="$CFLAGS_ARCH -O3 -fPIC"  && \
	(./configure --libdir=${LIB_DIR} $@ --prefix=${PRE_DIR} > $fnlog || \
	cmake -DCMAKE_INSTALL_PREFIX:PATH=${PRE_DIR} $@ > $fnlog) && \
	make -j4 install >> $fnlog) ;then
	echo "Compiled $name successfully."
    else
	echo "Failed to compile $name. Please see $fnlog for details"
    fi
}
function set_mkl_lapack(){
    #First use MKL if available. 
    LAPACK="-lmkl_core"
    if test -n "$OPENMP_CFLAGS" ;then
	if test "$CN" = "icc" ;then
	    LAPACK+=" -lmkl_intel_thread $OPENMP_CFLAGS"
	elif test "$CN" = "gcc" ;then
	    LAPACK+=" -lmkl_gnu_thread $OPENMP_CFLAGS"
	else
	    LAPACK+=" -lmkl_sequential -lpthread"
	fi
    else
	LAPACK+=" -lmkl_sequential -lpthread"
    fi
    if test "$host_cpu" = "x86_64" ;then
	LAPACK+=" -lmkl_intel_lp64"
    else
	LAPACK+=" -lmkl_intel"
    fi
    echo trying LAPACK=$LAPACK
}
set_mkl_lapack
AC_CHECK_LIB([m],[dposv_],[has_lapack=yes],[has_lapack=no],[$LAPACK])
if test "$has_lapack" = "no" ;then
    LAPACK="-lblas -llapack"
    echo trying LAPACK=$LAPACK

    unset nset ac_cv_lib_m_dgemm_
    unset nset ac_cv_lib_m_dposv_
    AC_CHECK_LIB([m],[dgemm_],[has_blas=yes],[has_blas=no],[$LAPACK])
    AC_CHECK_LIB([m],[dposv_],[has_lapack=yes],[has_lapack=no],[$LAPACK])
fi
if test "$has_blas" = "no" -o "$has_lapack" = "no" ;then
    unset nset ac_cv_lib_m_dgemm_
    unset nset ac_cv_lib_m_dposv_
    #Check for LAPACK library in the system
    if test -n "`which pkg-config 2>/dev/null`" ;then
	LAPACK=`pkg-config --libs blas lapack 2>/dev/null`
	if test -n "$LAPACK" ;then
	    AC_CHECK_LIB([m],[dgemm_],[has_blas=yes],[has_blas=no],[$LAPACK])
	    AC_CHECK_LIB([m],[dposv_],[has_lapack=yes],[has_lapack=no],[$LAPACK])
	fi
    fi
fi
#Some system has liblapack.so.? without symbol links to liblapack.so
if test "$has_blas" = "no" ;then
    fn_blas=`ls /usr/lib${libsuffix}/libblas.so.? 2>/dev/null`
    if test -n "${fn_blas}" ;then
	ln -sf "${fn_blas%% *}" "${LIB_DIR}/libblas.so"
    fi
    unset nset ac_cv_lib_m_dgemm_
    AC_CHECK_LIB([m],[dgemm_],[has_blas=yes],[has_blas=no],[$LAPACK])
fi
if test "$has_lapack" = "no" ;then
    fn_lapack=`ls /usr/lib${libsuffix}/liblapack.so.? 2>/dev/null`
    if test -n "${fn_lapack}" ;then
	ln -sf "${fn_lapack%% *}" "${LIB_DIR}/liblapack.so"
    fi
    unset nset ac_cv_lib_m_dposv_
    AC_CHECK_LIB([m],[dposv_],[has_lapack=yes],[has_lapack=no],[$LAPACK])
fi
if test "$has_lapack" = "no" -o "$has_blas" = "no" ;then
    if test "$mkl_avail" = "yes" ;then
	download MKL "mkl${libsuffix}.tar.bz2" ${PRE_DIR}
	set_mkl_lapack
	unset nset ac_cv_lib_m_dposv_
	AC_CHECK_LIB([m],[dposv_],[],[AC_MSG_ERROR([Please install blas/lapack])],[$LAPACK])
    else
	AC_MSG_ERROR([Please install blas/lapack])
    fi
fi

AC_CHECK_HEADERS([numa.h], [LDFLAGS+=" -lnuma"], [])
#Check for zlib
AC_CHECK_LIB([z],[gzopen],[],[has_z='no'])
if test "$has_z" = "no" ;then
    compile ZLIB zlib.tar.bz2 --static
    unset ac_cv_lib_z_gzopen
    AC_CHECK_LIB([z],[gzopen],[],[AC_MSG_ERROR([Invalid zlib])])
fi
LDZ="-lz"
#Check for FFTW3 header and library.
has_fftw='yes'
AC_CHECK_LIB([fftw3], [fftw_execute],[has_fftw="yes"],[has_fftw='no'])
AC_CHECK_HEADERS([fftw3.h],[], [has_fftw='no'])
if test "$has_fftw" = "no" ;then #fftw does not exist in system.
    compile FFTW fftw.tar.bz2 --enable-threads --enable-openmp --enable-shared
    unset ac_cv_lib_fftw3_fftw_execute
    #compile FFTW fftw.tar.bz2 --enable-float 
    AC_CHECK_LIB([fftw3], [fftw_execute],[],[AC_MSG_ERROR([Invalid FFTW library])])
fi
if test "$with_matlab" != "no" ;then
    #cannot compile with single floating fft when compiling mex routines
    has_fftwf=0
else
    AC_CHECK_LIB([fftw3f], [fftwf_execute],[has_fftwf=1],[has_fftwf=0])
fi
AC_DEFINE_UNQUOTED(HAS_FFTWF, [$has_fftwf], "FFTW single precision exists")

AC_DEFINE(USE_FFTW_THREADS, [0], "Disable FFTW3 threads support.")
if test x$has_fftwf = x1 ;then
    LDFFT="-lfftw3 -lfftw3f"
else
    LDFFT="-lfftw3"
fi

#Check for cholmod
AC_CHECK_LIB([cholmod], [${chol_fun}], [], [has_chol="no"], [$LAPACK -lm])
AC_CHECK_HEADERS([cholmod.h], [], [has_chol="no"])
if test "$has_chol" = "no" ;then
    compile CHOLMOD cholmod.tar.bz2
    unset ac_cv_lib_cholmod_${chol_fun}
    AC_CHECK_LIB([cholmod], [${chol_fun}], [], [AC_MSG_ERROR([Invalid CHOLMOD library])], [$LAPACK -lm])
fi
LDCHOL="-lcholmod"

#Check for libwebsocket
AM_CONDITIONAL(HAS_LWS, false)
if test "$use_websockets" = "yes" ;then
    AC_CHECK_LIB([websockets], [libwebsocket_service], [has_lws="yes"], [has_lws="no"], [])
    AC_CHECK_HEADERS([libwebsockets.h], [], [has_lws="no"])
    if test "$has_lws" = "no" ;then #first try to compile from source. Needs cmake
	compile LIBWEBSOCKETS libwebsockets.tar.bz2 -DLIB_SUFFIX=${libsuffix} -DLWS_WITH_SSL=OFF -DLWS_IPV6=OFF -DLWS_WITHOUT_TESTAPPS=ON -DLWS_WITHOUT_DAEMONIZE=ON -DLWS_WITHOUT_CLIENT=ON -DLWS_WITHOUT_DEBUG=ON
	unset ac_cv_lib_websockets_libwebsocket_service
	AC_CHECK_LIB([websockets], [libwebsocket_service], [has_lws="yes"], [has_lws="no"], [])
    fi
    if test "$has_lws" = "no" ;then #then download pre-compiled library and header.
	download LIBWEBSOCKETS "libwebsockets${libsuffix}.tar.bz2" ${PRE_DIR}
	unset ac_cv_lib_websockets_libwebsocket_service
	AC_CHECK_LIB([websockets], [libwebsocket_service], [has_lws="yes"], [has_lws="no"], [])
    fi
    if test "$has_lws" != "no" ;then
	AM_CONDITIONAL(HAS_LWS, true)
	AC_DEFINE(HAS_LWS, [1], "libwebsocket is available")
    else
	AC_MSG_NOTICE(["libwebsocket library or header is not available"])
    fi
fi
AM_CONDITIONAL(DRAW,false) #for drawdaemon
AM_CONDITIONAL(MON, false)  #for parital monitor without notify
AM_CONDITIONAL(MONN,false) #for full monitor with notify
if test "$with_gtk" != "no" ;then
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
    fi
    PKG_CHECK_MODULES(GTK3, gtk+-3.0,           [gtk_3="yes"], [gtk_3="no"])
    PKG_CHECK_MODULES(GTK2, gtk+-2.0 >= 2.2,    [gtk_2="yes"], [gtk_2="no"])

    #if test x$($PKG_CONFIG --variable=target gdk-2.0) = xquartz ;then
#	macint=gtk-mac-integration
    #fi
    if test "$gtk_3" = "yes" -a "$use_gtk_3" != "no" ;then
	PKG_CHECK_MODULES(DRAW, gtk+-3.0 >= 0.0 cairo >= 0.0,     [comp_draw="yes"],[comp_draw="no"])
	PKG_CHECK_MODULES(MON,  gtk+-3.0 >= 0.0,                  [comp_mon="yes"], [comp_mon="no"])
	PKG_CHECK_MODULES(MONN, gtk+-3.0 >= 0.0 libnotify >= 0.0, [comp_monn="yes"],[comp_monn="no"])
	gtk_ver=3
    else
	PKG_CHECK_MODULES(DRAW, gtk+-2.0 >= 2.2 $macint cairo >= 1.4,     [comp_draw="yes"],[comp_draw="no"])
	PKG_CHECK_MODULES(MON,  gtk+-2.0 >= 2.2 $macint ,                 [comp_mon="yes"], [comp_mon="no"])
	PKG_CHECK_MODULES(MONN, gtk+-2.0 >= 2.2 $macint libnotify >= 0.1, [comp_monn="yes"],[comp_monn="no"])
	gtk_ver=2
    fi
    if test -n "$GTKDIR" ;then
	TT=-I$GTKDIR/include
	GTK_CFLAGS="$TT $TT/cairo $TT/atk-1.0 $TT/gtk-${gtk_ver}.0 $TT/glib-2.0 $TT/pango-1.0 $TT/gdk-pixbuf-2.0 -I$GTKDIR/lib/glib-2.0/include -I${GTKDIR}/lib/gtk-2.0/include"
	GTK_LIBS="-L$GTKDIR/lib/ -Wl,-rpath,$GTKDIR/lib/"
	if test -d "$GTKDIR/bin" ;then
	    GTK_LIBS="$GTK_LIBS -L$GTKDIR/bin -Wl,-rpath,GTKDIR/bin"
	fi
    else
	GTK_CFLAGS=""
	GTK_LIBS=""
    fi

    AM_CONDITIONAL(DRAW, test x$comp_draw = xyes)
    AM_CONDITIONAL(MON,  test x$comp_mon  = xyes)
    AM_CONDITIONAL(MONN, test x$comp_monn = xyes)
    if test "$comp_mon" = "yes" ;then
	if ! echo "$MON_LIBS" | grep -q gthread ;then
	    #OLD GLIB has GTHREAD separated
	    PKG_CHECK_MODULES(GTHREAD, gthread-2.0 >= 2.0, [has_gthread="yes"], [has_gthread="no"])
	    if test -n "$GTHREAD_LIBS" ;then
		DRAW_LIBS="$DRAW_LIBS $GTHREAD_LIBS"
		MON_LIBS="$MON_LIBS $GTHREAD_LIBS"
		MONN_LIBS="$MONN_LIBS $GTHREAD_LIBS"
	    fi
	fi
    fi

    #The following will substitute the variables name 
    AC_SUBST(DRAW_CFLAGS)
    AC_SUBST(DRAW_LIBS)
    AC_SUBST(MON_CFLAGS)
    AC_SUBST(MON_LIBS)
    AC_SUBST(MONN_LIBS)
    AC_SUBST(MONN_CFLAGS)
    AC_SUBST(GTK_CFLAGS)
    AC_SUBST(GTK_LIBS)
fi

AC_CHECK_HEADERS([unistd.h])
LDOPTS="$LDCHOL $LDFFT $LAPACK $LDZ -lm $LDEXE"
OPTS="$OPTScommon $OPTScc $OPTSprof $OPTSextra $OPENMP_CFLAGS"
AC_SUBST(LDAOS)
AC_SUBST(OPTS)
AC_SUBST(LDOPTS)
AC_SUBST(OPTScommon)
AC_SUBST(SRC_DIR)
AC_SUBST(BUILD_DIR)

if test "$with_matlab" != "no" ;then
    if test -x "${MEXDIR}/mex" ;then
	MEXDIR="${MEXDIR}"
    elif test -x "${MEXDIR}/bin/mex" ;then
	MEXDIR="${MEXDIR}/bin"
    else
	AC_MSG_ERROR([mex not found. please reconfigure with --with-matlab=DIR to specify MATLAB location or --without-mex to disable])
    fi
    MEXEXT="$($MEXDIR/mexext)"
    if test -z "$MEXEXT" ;then
	AC_MSG_ERROR([Unable to determine mex extension. reconfigure with --disable-mex to disable mex.])
    fi
    case $MEXEXT in
   	mexglx)
	    MEXLIB=glnx86
	    MEXOPT="$OPTS -fPIC -m32 -fexceptions -D_FILE_OFFSET_BITS=64"
	    MEXOPT="${MEXOPT} -D_BSD_SOURCE -D_POSIX_C_SOURCE=200112L"
	    mex_cpu=i386
	    ;;
   	mexa64)
	    MEXLIB=glnxa64
	    MEXOPT="$OPTS -fPIC -m64 -fexceptions -fno-omit-frame-pointer"
	    MEXOPT="${MEXOPT} -DDLONG -D_BSD_SOURCE -D_POSIX_C_SOURCE=200112L"
	    mex_cpu=x86_64
	    ;;
	mexmaci)
	    MEXLIB=maci
	    MEXOPT="$OPTS -fPIC -m32 -fexceptions" 
	    mex_cpu=i386
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
	    mex_cpu=i386
	    ;;
	*)
	    AC_MSG_ERROR([Unknown mex extension "$MEXEXT"])
	    ;;
    esac
    if test "$mex_cpu" != "$host_cpu" ;then
	mex_lib="no"
    else       
	LIBMEX="-L${MEXDIR}/${MEXLIB} -Wl,-rpath,${MEXDIR}/${MEXLIB} "
	AC_CHECK_LIB([mx], [mxGetPr], [mex_lib="yes"], [mex_lib="no"], [$LIBMEX])
    fi
    if test x$mex_lib = xyes ;then
    	MEXOPT="${MEXOPT} -I${MEXDIR}/../extern/include -DMATLAB_MEX_FILE $OPENMP_CFLAGS" 
	LIBMEX+=" -lmwlapack -lmwblas -lmx -lmex -lmat -lut -lstdc++ -lz -lm -ldl"
	LDMEX="-no-fast-install -static -Xcompiler -shared"
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

if test -z "$CUDADIR" ;then
    with_cuda="no"
fi
if test "$with_cuda" != "no" ;then
    if test -n "$with_cuda" ;then
	user_enable_cuda=yes
    fi
    if test "$use_static" = "yes"; then
	AC_MSG_ERROR("all-static conflicts with cuda support. Please disable all-static")
    fi
    AC_MSG_NOTICE([Enabling CUDA])
    #Check for nvcc executable
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
            AC_MSG_NOTICE(["Cannot compile with cuda. Please use --with-cuda=DIR in configure to specify CUDA location or disable cuda use --without-cuda"])
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
    if test "$cudaver" -lt 40 ;then
	AC_MSG_ERROR([Require at least CUDA toolkit 4.0 ])
    fi
    if test "$cudaver" -ge 40 ;then
	CUDAOPT+=" -gencode arch=compute_20,code=sm_20"
    fi
    if test "$cudaver" -ge 42 ;then
	CUDAOPT+=" -gencode arch=compute_30,code=sm_30"
    fi
    if test "$cudaver" -ge 65 ;then
        CUDAOPT="$CUDAOPT -gencode arch=compute_50,code=sm_50"
    fi
    #if test -n "$OPENMP_CFLAGS";then
    #    CUDAOPT="$CUDAOPT -Xcompiler $OPENMP_CFLAGS"
    #fi
    if test "$debug" = "yes" ;then
        CUDAOPT+=" -O0"
    else
	CUDAOPT+=" -O3"
    fi
    CUDAOPT+=" -g -lineinfo -DHAVE_CONFIG_H -I.. "
    #CUDAOPT+="-Xptxas=-v -dlcm=cg"
    if test "$host_cpu" = "x86_64" ;then
        CUDAOPT+=" -DDLONG"
    fi
    CCBINOPT="$OPTS -fPIC "
    case "$CCBIN" in 
        *gcc*)
	CCBINOPT+=" -Wno-unused-value -D__STRICT_ANSI__"
  	;;			  
    esac
    LIBCUDA+=" -lcudart -lcurand -lcusparse -lcufft -lcublas -lstdc++"
else
    with_cuda="no"
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

AC_CONFIG_FILES( [Makefile
	sys/Makefile
	math/Makefile
	lib/Makefile
	tools/Makefile
	maos/Makefile 
	skyc/Makefile
	mex/Makefile
	cuda/Makefile
	])
if test -d "${srcdir}/test" ;then
    AC_CONFIG_FILES([test/Makefile])
fi
LDFLAGS+=" -no-fast-install"
LIBS="" #somehow LIBS contain -fftw3 -lnvidia-sml
AC_SUBST(LIBS)
AC_OUTPUT
echo SRC_DIR=$SRC_DIR
echo BUILD_DIR=$BUILD_DIR
echo host_os=$host_os
echo BIGENDIAN=$BIGENDIAN
echo use_static=$use_static

echo CC=$CC
echo CFLAGS=$CFLAGS
echo LDFLAGS=$LDFLAGS
echo OPTS=$OPTS
echo LDOPTS=$LDOPTS

if test "$with_gtk" != "no" ;then
    echo Compile Drawdaemon=$comp_draw
    echo Compile Monitor=$comp_mon
    echo Compile Monitor with Libnotify=$comp_monn
    echo DRAW_LIBS=$DRAW_LIBS
    echo GTK_LIBS=$GTK_LIBS
    echo Using GTK+-${gtk_ver}
else
    echo "Without GTK"
fi
echo Compile websocket server=$has_lws
if test "$with_cuda" != "no" ;then
    echo NVCC=$NVCC
    echo CUDAOPT=$CUDAOPT
    echo LIBCUDA=$LIBCUDA
else
    echo "Without CUDA support"
fi
if test "$with_matlab" != "no" ;then
    echo MEXDIR=$MEXDIR
    echo LIBMEX=$LIBMEX
    echo MEXOPT=$MEXOPT
else
    echo "Without matlab support"
fi

echo "*************************************************************************"
echo "*** MAOS Documentation can be found at http://lianqiw.github.com/maos ***"
echo "*** MAOS Source can be obtained at http://github.com/lianqiw/maos     ***"
echo "*************************************************************************"

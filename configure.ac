AC_PREREQ(2.57)
AC_INIT(maos, 1.7.0, lianqiw@gmail.com)
AM_INIT_AUTOMAKE([1.9.6 gnu -Wall no-define subdir-objects])
#silent building
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])]) 
#The URL where misc supportive libaries are located.
BASEURLS="http://project.tmt.org/~lianqiw/maos/"
#If user does not supply any CFLAGS, set something so that AC_PROG_CC won't put a -g -O2 to CFLAGS
if test -z "$CFLAGS" ;then
   CFLAGS=" "
fi
if test -z "$LIBS" ;then
   LIBS=" "
fi
AC_LANG([C])
AC_CONFIG_MACRO_DIR([m4])
#After enable HEADER, the package name will not appear in CFLAGS. Appear in config.h instead.
AC_CONFIG_HEADER([config.h])
#defines host host_os host_cpu host_vendor etc
AC_CANONICAL_HOST 
#In darwin it reports wrong host_cpu (i386). The following is correct in Linux and SL.
if which uname 2>/dev/null;then
   host_cpu=$(uname -m)
fi
AC_C_BIGENDIAN([BIGENDIAN=1], [BIGENDIAN=0], AC_MSG_ERROR("Unknown Endianness"))
AC_DEFINE_UNQUOTED(BIGENDIAN, [$BIGENDIAN], ["Endianness"])

#top_srcdir and abs_top_srcdir is empty. no idea why. but srcdir sometimes is the full path.
SRC_DIR=$(cd "$srcdir" && pwd)
BUILD_DIR=$(pwd)
CFLAGS="$CFLAGS -I${SRC_DIR}/external/include"
LDFLAGS="$LDFLAGS -lm -L${SRC_DIR}/external/lib -Wl,-rpath,${SRC_DIR}/external/lib"
#AC_DEFINE_UNQUOTED means expand the $ variable
AC_DEFINE_UNQUOTED(SRCDIR,  ["$SRC_DIR"], ["Define SRC folder"])
AC_DEFINE_UNQUOTED(BUILDDIR,["$BUILD_DIR"],["Define BUILD folder"])
#Specifying both -O3 and -g in mac produces internal error in gcc 4.2
#-all-static will produce a static exe.
OPTScommon="-Wno-unknown-pragmas -Wno-unused-function" 
if test "${USER}" = "lianqiw" ;then
   OPTScommon="-Wshadow -Wall -Werror $OPTScommon"
fi
LDEXE="-lm -no-fast-install -ldl -lpthread"
OPTSextra=""

#Compile CUDA branch if CUDA toolkit is detected
use_cuda_default="no"
CUDADIR_AUTO=
if test -n "`which nvcc 2>/dev/null`" ;then
   use_cuda_default="yes"
   NVCC="$(which nvcc)"
   if test -L "$NVCC" ;then
      NVCC="$(readlink "$NVCC")"
   fi
   CUDADIR_AUTO="$(dirname "$(dirname "${NVCC}")")"
   unset NVCC
elif test -x "/opt/cuda/bin/nvcc" ;then
   use_cuda_default="yes"
   CUDADIR_AUTO="/opt/cuda"
elif test -x "/usr/local/cuda/bin/nvcc" ;then
   use_cuda_default="yes"
   CUDADIR_AUTO="/usr/local/cuda/"
fi


#Enable the user to force 32 or 64 bit program
AC_ARG_ENABLE(32,  [  --enable-32   Compile 32 bit executables], [use_32="$enableval"], [use_32="no"])
AC_ARG_ENABLE(64,  [  --enable-64   Compile 64 bit executables], [use_64="$enableval"], [use_64="no"])
cross_compile="no"
if test "$use_64" = "yes" ;then
   case $host_cpu in
      *86)
         cross_compile="yes"
      ;;
   esac
   host_cpu=x86_64
elif test "$use_32" = "yes" ;then
   case $host_cpu in
      *64)
         cross_compile="yes"
      ;;
   esac
   host_cpu=i386
fi
case ${host_cpu} in 
   *64)
	host_cpu=x86_64
        cpu="64"
        CFLAGS="$CFLAGS -m64"
        LDFLAGS="$LDFLAGS -m64"
        CDEFINE="__x86_64__"
	AC_DEFINE(X86_64, [1], "X86_64 Machine")
	objcopy_output=elf64-x86-64
	objcopy_arch=i386:x86-64
	echo "CPU is 64 bit"
   ;;
   *86)
	host_cpu=i386
        cpu="32"
	CFLAGS="$CFLAGS -m32"
	LDFLAGS="$LDFLAGS -m32"
        CDEFINE="__i386__"
	AC_DEFINE(X86, [1], "I386 Machine")
	objcopy_output=elf32-i386
	objcopy_arch=i386
	echo "CPU is 32 bit"
   ;;
   *)
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

win32=no
use_mkl="no"
mkl_avail="no"
case "$host_os" in
*linux*)
 AM_CONDITIONAL(LINUX, true)
 LDEXE="$LDEXE -lrt" #need real time library.
 OPTSextra="-g -D_BSD_SOURCE -D_POSIX_C_SOURCE=200112L"
 mkl_avail="yes"
 libsuffix=${cpu}
 system="linux"
 ldsuffix="so"
 ;;
*darwin*)
 AM_CONDITIONAL(DARWIN, true)
 mkl_avail="no"
 libsuffix=fat
 system="apple"
 ldsuffix="dylib"
 #set CFLAGS=-mmacosx-version-min=10.6 before calling configure for compatibility with OS X >=10.6. Don't use 10.5 which will render mex files not usable.
 ;;
*bsd*)
 AM_CONDITIONAL(BSD, true)
 system="bsd"
 ldsuffix="so"
 ;;
*cygwin*|*mingw*)
 AM_CONDITIONAL(WIN, true)
 LAPACK="-lblas -llapack"
 win32=yes
 system="win"
 ldsuffix="dll"
 ;;
*)
AC_MSG_ERROR([Unrecognized OS])
esac
AC_SUBST(ldsuffix)

if test -n "`which wget 2>/dev/null`" ;then
   wget="wget -t1 --no-check-certificate"
   wgetoutput="-O"
elif test -n "`which curl 2>/dev/null`" ;then
   wget="curl --retry 1 -L --progress-bar"
   wgetoutput="--output"
else
   wget="(echo wget/curl not found ; exit)"
fi
#Common configure options.
AC_ARG_ENABLE(all-static,[  --diable-all-static Compile statically], [use_static="$enableval"],[use_static="no"])
AC_ARG_ENABLE(debug, [  --enable-debug  Turn on debug], [debug="$enableval"], [debug="no"])
AC_ARG_ENABLE(mem,   [  --enable-mem    Turn on memory management],[use_mem="$enableval"],[use_mem="no"])
AC_ARG_ENABLE(oprof, [  --enable-oprof  Enable oprof for code profiling.],[use_oprof="$enableval"],[use_oprf="no"])
AC_ARG_ENABLE(gnu-inline, [  --enable-gnu-inline  Force GNU type of inline ], [gnu_inline="$enableval"], [gnu_inline="no"])

AC_ARG_ENABLE(cuda,  [  --enable-cuda   Enable Cuda], [use_cuda="$enableval"], [use_cuda="$use_cuda_default"])
AC_ARG_ENABLE(cuda-double, [  --enable-cuda-double   Enable cuda double precision], [cuda_double="$enableval"],[cuda_double="no"])
AC_ARG_WITH(cuda, --with-cuda=DIR, [CUDADIR=${withval} use_cuda="yes"], [CUDADIR="$CUDADIR_AUTO"])
AC_ARG_WITH(ccbin, --with-ccbin=CC, [CCBIN=${withval}])
AC_ARG_ENABLE(scheduler,[  --disable-scheduler  Disable built in scheduler (default: enabled)], [use_scheduler="$enableval"], [use_scheduler="yes"])
AC_ARG_ENABLE(websockets,[  --disable-websockets  Disable websockets (default: enabled)], [use_websockets="$enableval"], [use_websockets="yes"])
AC_ARG_WITH(gtk,  --with-gtk=DIR, [GTKDIR=${withval}], [GTKDIR=""])

AM_CONDITIONAL(USE_SCHEDULER, true)
if test "$use_scheduler" = "no" ;then
   AC_DEFINE(MAOS_DISABLE_SCHEDULER, [1], "Disable built in scheduler")
   AM_CONDITIONAL(USE_SCHEDULER, false)
fi
if test "$cpu" = "64" ;then 
   OPTScommon="$OPTScommon -DDLONG"
   chol_fun=cholmod_l_factorize
   AM_CONDITIONAL(USE_LONG, true)
else
   chol_fun=cholmod_factorize
   AM_CONDITIONAL(USE_LONG, false)
fi

#Compile mex routines if MATLAB is detected and single is not enabled
use_mex_default="no"
MEXDIR_AUTO=
if test -n "`which mex 2>/dev/null`" ;then
   use_mex_default="yes"
   MEX="$(which mex)"
   if test -L "$MEX" ;then
      MEX="$(readlink "$MEX")"
   fi
   MEXDIR_AUTO="$(dirname "$MEX")"
elif test -n "`ls -d /Applications/MATLAB*.app 2>/dev/null`" ;then
   for i in "$(ls -d /Applications/MATLAB*.app)" ;do
       if test -x "$i/bin/mex" ;then
       	  MEXDIR_AUTO="${i}/bin"
          use_mex_default="yes"
       fi
   done
fi
AC_ARG_ENABLE(mex,   [  --enable-mex    Enable compiling of matlab mex routines, requires mex to be in the PATH],[use_mex="$enableval"],[use_mex="$use_mex_default"])
AC_ARG_WITH(matlab, --with-matlab=DIR, [MEXDIR=${withval} use_mex="yes" ], [MEXDIR=${MEXDIR_AUTO}])

if test "x$use_static" = "xyes" ;then
   LDEXE="$LDEXE -all-static" #Enable all static build.
   AM_CONDITIONAL(USE_STATIC, true)
   AC_DEFINE(USE_STATIC, [1], "We are compling maos statically")
   #(cd ${SRC_DIR} && tar zcf ${BUILD_DIR}/config.tar.gz config/)
else
   LDEXE="$LDEXE -static" #show that we prefer static linking.
   AM_CONDITIONAL(USE_STATIC, false)
   AC_DEFINE(USE_STATIC, [0], "We are compling maos dynamically")
fi

if test x$use_oprof = "xyes" ;then
   OPTSprof="$OPTSprof -g -fno-omit-frame-pointer -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls"
fi
#Set optimize flags for different compilers
#-xO Enables SSE3, SSE2 and SSE instruction sets optimizations for non-Intel CPUs 
#-wd1292 : disables warning "error #1292: attribute "__alloc_size__" ignored"
COMMON_ICC="-wd193 -wd1418 -wd2259 -wd1372 -wd-wd981 -wd1292 -wd10237 -diag-disable remark -vec-report=3"
OPTIM_ICC="-O3 -funroll-loops -DHAVE_SSE2=1 -msse2 " #-static-intel"
OPTIM_GCC="-O3 -funroll-loops -DHAVE_SSE2=1 -msse2 -msse -fprefetch-loop-arrays -fno-math-errno -ffinite-math-only -static-libgcc -Wno-deprecated-declarations" #-Mcache_align"
OPTIM_CLANG="-O3 -DHAVE_SSE2=1 -msse2 -msse -Qunused-arguments"
OPTIM_OPENCC="-O3"
OPT_DEBUG="-g -ggdb -O0 -DDEBUG_TRACE"
if test "x$use_mem" = "xyes" ;then
    AC_DEFINE(USE_MEM,[1],["Use memory management"])
    debug="yes"
else
    AC_DEFINE(USE_MEM,[0],["Do not use memory management"])	
fi
AC_PROG_CC([$CC])
AC_PROG_CXX([$CXX])
AC_OPENMP
AC_TYPE_SIZE_T #it sets CC if empty
#inline definition:

#in GCC, extern inline will use the inline and not emit function definition.
#in GCC, inline will cause use the inline, but also emit function deminition, causing multiple definitions if put in header.
#in C99, extern inline will cause use the inline, but also emit function deminition, causing multiple definitions if put in header. This is equivalent to GCC inline
#in C99, inline will use the inline and not emit function definition. This is equivalent to GCC extern inline 
#GCC mode can be forced on using __attribute__((gnu_inline))
#always_inline causes inline to be done in debug mode.
cc_full=`which ${CC}`
while test -L $cc_full ;do
   cc_full=`readlink $cc_full`
done
case "${cc_full}" in	
  *icc*)
    OPTIM="$OPTIM_ICC $COMMON_ICC"
    OPT_DEBUG="$OPT_DEBUG $COMMON_ICC"
    which icc 2>/dev/null|| AC_MSG_ERROR([icc is not found in path.])
    #Test ICC version
    CC_VERSION=$($CC -dumpversion)
    if test -n "$CC_VERSION" ;then
	major=${CC_VERSION%%.*}
	minor=${CC_VERSION%.*}
	minor="${minor#4.}"
	if test $major -ge 12 ;then
	    OPTIM="$OPTIM -axSSE4.2"
        elif test $major -ge 11 ;then
	    OPTIM="$OPTIM -axSSE4.2"
        elif test $major -ge 10 ;then
	    OPTIM="$OPTIM -axS"
	else
	    echo "icc version is old"
	fi
    fi
    ;;
  *gcc*)
    if test -n "`$CC -v 2>&1 |grep clang`" ;then
        OPTIM="$OPTIM_CLANG"
    else
	OPTIM="$OPTIM_GCC"
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
	    echo "Enabling -mtune=native"
	    OPTIM="$OPTIM -mtune=native" #introduced in GCC 4.2
	    OPTIM="$OPTIM -ftree-vectorize " #introduced in GCC 4
	fi
	 #GCC4.3 changed meaning of inline to C99
	if test $major -eq 4 -a $minor -lt 5;then
	    gnu_inline="yes"
	fi
    fi
    ;;
  *clang*)
	OPTIM="$OPTIM_CLANG"
    ;;
  *) #Other compilers
    OPTIM="-O3"
    ;;
esac
if test "$gnu_inline" = "yes" ;then
INLINE="extern __inline __attribute__((gnu_inline, always_inline))" #C99
else
INLINE="inline __attribute__((always_inline))" #C99
fi
AC_DEFINE_UNQUOTED(INLINE, $INLINE, [Portable definition of inline])
if test "x$debug" = "xyes" ;then
    AM_CONDITIONAL(DEBUG,true)
    AC_DEFINE(DEBUG,[1],["Debug the code"])
    OPTScc="$OPT_DEBUG"
else
    AM_CONDITIONAL(DEBUG,false)
    AC_DEFINE(DEBUG,[0],["Do not debug"])
    OPTScc="$OPTIM"
fi
AM_PROG_AR
AC_PROG_CC_STDC
AM_PROG_CC_C_O
AC_PROG_CC_C99 #enable C99.

#call convention of fortran libraries
#AC_F77_FUNC(dgemm)
#AC_F77_FUNC(dposv)
#download LIBNANE file.tar.bz2 [outdir]
download(){
   echo The $1 library is not available in the system
   echo Will download one from the website.
   fn="$2"
 if test ! -f /tmp/${fn} ;then
   success=0
   for BASEURL in $BASEURLS ; do
      if $wget "${BASEURL}/${fn}" "$wgetoutput" "/tmp/${fn}" ;then
      	 success=1
	 break;
      else
         AC_MSG_NOTICE([Failed to download the file from ${BASEURL}/${fn}])
      fi
   done
   if test $success -eq 0 ;then
      AC_MSG_ERROR([Please download it manually and put in /tmp/${fn}])
   fi
 fi
   if test -n "$3" ;then
      OUTDIR="$3"
   else
      AC_MSG_ERROR([Please specify third argument])
   fi
   if test ! -d "${OUTDIR}" ;then
      mkdir -p "${OUTDIR}" || exit 1
   fi
   tar -xf "/tmp/${fn}" -C "${OUTDIR}"  || exit 1
}
#compile LIBNAME srcfile.tar.gz [compile options]
compile(){
   echo Download and compile $1
   name=$1
   file=$2
   shift
   shift
   opts="$@"
   if test ! -d ${SRC_DIR}/external/$name ;then
      download $name $file ${SRC_DIR}/external/$name
   fi
   fnlog=${SRC_DIR}/external/compile.log
   if (cd ${SRC_DIR}/external/$name/* && \
      (./configure "$@" --prefix=${SRC_DIR}/external > $fnlog || \
      cmake -DCMAKE_INSTALL_PREFIX:PATH=${SRC_DIR}/external "$@" > $fnlog) && \
      make -j4 install >> $fnlog) ;then
      echo "Compiled $name successfully."
   else
      echo "Failed to compile $name. Please see $fnlog for details"
      exit
   fi
}

#Check for LAPACK library in the system
if test -z "$LAPACK" ;then
   LAPACK=`pkg-config --libs blas lapack`
   if test -z "$LAPACK" ;then
       LAPACK="-lblas -llapack"
   fi
fi
AC_CHECK_LIB([blas],  [dgemm_],[has_lapack=yes],[has_lapack=no],[$LAPACK])
AC_CHECK_LIB([lapack],[dposv_],[],[has_lapack=no],[$LAPACK])
if test "$has_lapack" = "no" ;then
   if test "$mkl_avail" = "yes" ;then
       mklname=mkl${libsuffix}
       download MKL "${mklname}.tar.bz2" /tmp/MKL
       mv /tmp/MKL/*/*.* ${SRC_DIR}/external/lib
       unset ac_cv_lib_blas_dgemm_
       AC_CHECK_LIB([blas],  [dgemm_],[],[AC_MSG_ERROR([Please install blas/lapack])],[$LAPACK])
   else
       AC_MSG_ERROR([Please install blas/lapack])
   fi
fi

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
   compile FFTW fftw.tar.bz2
   unset ac_cv_lib_fftw3_fftw_execute
   #compile FFTW fftw.tar.bz2 --enable-float 
   AC_CHECK_LIB([fftw3], [fftw_execute],[],[AC_MSG_ERROR([Invalid FFTW library])])
fi
LDFFT="-lfftw3"
#AC_CHECK_LIB([fftw3f], [fftwf_execute],[has_fftwf=1],[has_fftwf=0],[$LDFFT])
#AC_DEFINE_UNQUOTED(HAS_FFTWF, [$has_fftwf], "FFTW single precision exists)

AC_ARG_ENABLE(fftw-threads, [ --disable-fftw-threads Disable FFTW threads], use_fftw_threads="$enableval", use_fftw_threads="no")
if test "$use_fftw_threads" = "yes" ;then
     AC_CHECK_LIB([fftw3_threads], [fftw_init_threads], [has_fftw_threads="yes"],[has_fftw_threads="no"],[$LDFFT])
     AC_CHECK_LIB([fftw3_threads], [omp_get_num_threads], [has_fftw_omp="yes"], [has_fftw_omp="no"],[$LDFFT])
else
     has_fftw_threads="no"
fi
if test "x$has_fftw_threads" = "xyes" -a "x$has_fftw_omp" = "xno" ;then
    AC_DEFINE(USE_FFTW_THREADS, [1], "We have FFTW3 threads support")
    if test x$has_fftwf = x1 ;then
        LDFFT="$LDFFT -lfftw3_threads -lfftw3 -lfftw3f_threads -lfftw3f"
    else
        LDFFT="$LDFFT -lfftw3_threads -lfftw3"
    fi
else
    AC_DEFINE(USE_FFTW_THREADS, [0], "Disable FFTW3 threads support.")
    if test x$has_fftwf = x1 ;then
        LDFFT="$LDFFT -lfftw3 -lfftw3f"
    else
        LDFFT="$LDFFT -lfftw3"
    fi
fi

#Check for cholmod
AC_CHECK_LIB([cholmod], [${chol_fun}], [], [has_chol="no"], [$LAPACK])
AC_CHECK_HEADERS([cholmod.h], [], [has_chol="no"])
if test "$has_chol" = "no" ;then
   compile CHOLMOD cholmod.tar.bz2
   unset ac_cv_lib_cholmod_${chol_fun}
   AC_CHECK_LIB([cholmod], [${chol_fun}], [], [AC_MSG_ERROR([Invalid CHOLMOD library])], [$LAPACK])
fi
LDCHOL="-lcholmod"

#Check for libwebsocket
AM_CONDITIONAL(HAS_LWS, false)
if test "$use_websockets" = "yes" ;then
AC_CHECK_LIB([websockets], [libwebsocket_service], [has_lws="yes"], [has_lws="no"], [])
AC_CHECK_HEADERS([libwebsockets.h], [], [has_lws="no"])
if test "$has_lws" = "no" ;then
   compile LIBWEBSOCKETS libwebsockets.tar.bz2 -DLWS_WITH_SSL=OFF -DLWS_IPV6=OFF -DLWS_WITHOUT_CLIENT -DLWS_WITHOUT_SERVER
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

#Check availability of gtk and libnotify. If available, we will build drawdaemon and the monitor.
AM_CONDITIONAL(DRAW,true) #for drawdaemon
AM_CONDITIONAL(MON,true)  #for parital monitor without notify
AM_CONDITIONAL(MONN,true) #for full monitor with notify
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
if test "$gtk_2" = "yes" -a  "$gtk_3" = "yes" ;then 
   #both gtk2 & gtk3 exist in system. let use choose
   AC_ARG_ENABLE(gtk-3, [  --enable-gtk-3  Using GTK+-3 ], [gtk_3="$enableval"], [gtk_3="no"])
fi
#The following will define MON_CFLAGS and MON_LIBS with their values in Makefile.in,
#so that you can use $(DRAW_CFLAGS) in Makefile.am to represent the values of the shell variable here.
if test x$($PKG_CONFIG --variable=target gdk-2.0) = xquartz ;then
   macint=gtk-mac-integration
fi
if test "$gtk_3" = "yes" ;then
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

AC_CHECK_HEADERS([unistd.h])
LDLIB=""
LDOPTS="$LDCHOL $LDFFT $LAPACK $LDZ -lm"
OPTS="$OPTScommon $OPTScc $OPTSprof $OPTSextra $OPENMP_CFLAGS"
AC_SUBST(LDAOS)
AC_SUBST(OPTS)
AC_SUBST(LDOPTS)
AC_SUBST(LDEXE)
AC_SUBST(LDLIB)
AC_SUBST(OPTScommon)
AC_SUBST(SRC_DIR)
AC_SUBST(BUILD_DIR)
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
#The following must be after AC_PROG_CC
AC_PROG_MAKE_SET
AC_PROG_INSTALL
AC_ENABLE_SHARED
AC_ENABLE_STATIC
AC_PROG_LIBTOOL

if test "x${use_mex}" = "xyes" ; then
   AC_MSG_NOTICE([MEXDIR=${MEXDIR}])
   if test -n "$MEXDIR" ;then
       if test -x "${MEXDIR}/mex" ;then
           MEXDIR="${MEXDIR}"
       elif test -x "${MEXDIR}/bin/mex" ;then
           MEXDIR="${MEXDIR}/bin"
       else
           AC_MSG_ERROR([mex not found. please reconfigure with --with-matlab=DIR to specify MATLAB location or --disable-mex to disable])
       fi
   else
       AC_MSG_ERROR([mex not found. please reconfigure with --with-matlab=DIR to specify MATLAB location or --disable-mex to disable])
   fi
   MEXEXT="$($MEXDIR/mexext)"
   if test -z "$MEXEXT" ;then
       AC_MSG_ERROR([Unable to determine mex extension. reconfigure with --disable-mex to disable mex.])
   fi
   case $MEXEXT in
   	mexglx)
		MEXLIB=glnx86
		MEXOPT="-O3 -DNDEBUG -fPIC -pthread -m32 -fexceptions -D_FILE_OFFSET_BITS=64"
		MEXOPT="${MEXOPT} -D_BSD_SOURCE -D_POSIX_C_SOURCE=200112L"
		mex_cpu=i386
	;;
   	mexa64)
		MEXLIB=glnxa64
		MEXOPT="-O3 -DNDEBUG -fPIC -pthread -m64 -fexceptions -fno-omit-frame-pointer"
		MEXOPT="${MEXOPT} -DDLONG -D_BSD_SOURCE -D_POSIX_C_SOURCE=200112L"
		mex_cpu=x86_64
	;;
	mexmaci)
		MEXLIB=maci
		MEXOPT="-O3 -DNDEBUG -fno-common -m32 -fexceptions" #omitted -no-cpp-precomp
		mex_cpu=i386
	;;
	mexmaci64)
		MEXLIB=maci64
		MEXOPT="-O3 -DNDEBUG -fno-common -fexceptions -arch x86_64" #omitted -no-cpp-precomp
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
       LDMEX="-L${MEXDIR}/${MEXLIB} -Wl,-rpath,${MEXDIR}/${MEXLIB} -lmx -lmex -lmat -lz"
       AC_CHECK_LIB([mex], [mxGetPr], [mex_lib="yes"], [mex_lib="no"], [$LDMEX])
   fi
   if test x$mex_lib = xyes ;then
       AM_CONDITIONAL(USE_MEX, true)
   else
       AC_MSG_NOTICE([mex library test failed. Disable mex])
       AM_CONDITIONAL(USE_MEX, false)
   fi
   MEXOPT="${MEXOPT} -I${MEXDIR}/../extern/include -DMATLAB_MEX_FILE"
   AC_SUBST(MEXEXT)
   AC_SUBST(LDMEX)
   AC_SUBST(MEXOPT)
else
   AM_CONDITIONAL(USE_MEX, false)
fi

if test "$use_cuda" = "yes" ;then
   if test "$use_static" = "yes"; then
      AC_MSG_ERROR("all-static conflicts with cuda support. Please disable all-static")
   fi
 
   AC_MSG_NOTICE([Enabling CUDA])
   LDCUDA=
   if test -n "$CUDADIR" ;then
       NVCC="$CUDADIR/bin/nvcc"
       if test ! -x "$NVCC" ;then
       	  NVCC="$CUDADIR/nvcc"
	  CUDADIR="$(dirname $CUDADIR)"
       fi
   else
       AC_MSG_ERROR([CUDADIR is empty. Please use --with-cuda=DIR to specify CUDA location or --disable-cuda to disable cuda])
   fi
   if test ! -x "$NVCC" ;then
       AC_MSG_ERROR([nvcc not found])
   fi
   AC_MSG_NOTICE([Trying cuda in directory $CUDADIR])
   if test -f "${CUDADIR}/lib${cpu}/libcudart.$ldsuffix" ;then
       CUDA_L="${CUDADIR}/lib${cpu}"
   elif test -f "${CUDADIR}/lib/libcudart.$ldsuffix" ;then
       CUDA_L="${CUDADIR}/lib"
   fi
   if test -d "$CUDA_L" ;then
       LDCUDA="-L$CUDA_L -Wl,-rpath,$CUDA_L"
   else
       LDCUDA=
   fi
   if test -f "${CUDADIR}/include/cuda.h" ;then
       CUDAOPT="-I${CUDADIR}/include"
   else
       AC_CHECK_HEADERS([cuda.h], [] ,[AC_MSG_ERROR([Please use --with-cuda=DIR to specify CUDA location or --disable-cuda to disable cuda])])
       CUDAOPT=
   fi
   cudaver="$($NVCC -V |grep release |cut -d ',' -f 2 |cut -d ' ' -f 3)"
   cudaver="${cudaver/.}"
   echo CUDA toolkit is $cudaver
   if test -z "$CCBIN" ;then
       CCBIN="$CC"
       NVCC="$NVCC -ccbin ${CCBIN%% *}"
   fi  
   
   if test "$cudaver" -lt 40 ;then
       AC_MSG_ERROR([Require at least CUDA toolkit 4.0 ])
   fi
   if test "$cudaver" -ge 40 ;then
       CUDAOPT="$CUDAOPT -gencode arch=compute_20,code=sm_20"
   fi
   if test "$cudaver" -ge 42 ;then
       CUDAOPT="$CUDAOPT -gencode arch=compute_30,code=sm_30"
   fi
   #if test "$cudaver" -ge 50 ;then
   #    CUDAOPT="$CUDAOPT -gencode arch=compute_35,code=sm_35"
   #fi
   if test -n "$OPENMP_CFLAGS";then
       CUDAOPT="$CUDAOPT -Xcompiler $OPENMP_CFLAGS"
   fi
   echo LDCUDA="$LDCUDA"
   echo CUDAOPT="$CUDAOPT"
   AC_CHECK_LIB([cudart], [cudaSetDevice], [cuda_lib="yes"], [cuda_lib="no"], [$LDCUDA])
   if test x$cuda_lib = xno ;then
       AC_MSG_NOTICE([libcudart not found or currect architecture. Please use --with-cuda=DIR to specify CUDA location])
       use_cuda="no"
   fi
   LDCUDA="$LDCUDA -lcudart -lcurand -lcusparse -lcufft -lcublas -lstdc++"
fi
if test "$use_cuda" = "yes" ;then
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
   AC_SUBST(LDCUDA)
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

LIBS="" #somehow LIBS contain -fftw3 -lnvidia-sml
AC_SUBST(LIBS)
AC_OUTPUT
echo SRC_DIR=$SRC_DIR
echo BUILD_DIR=$BUILD_DIR
echo host_os=$host_os
echo use_mex=$use_mex
echo use_cuda=$use_cuda
echo use_static=$use_static
echo Cross compile=$cross_compile
echo Compile Drawdaemon=$comp_draw
echo Compile Monitor=$comp_mon
echo Compile Monitor with Libnotify=$comp_monn
echo DRAW_LIBS=$DRAW_LIBS
echo GTK_LIBS=$GTK_LIBS
echo Using GTK+-3=$gtk_3
echo CC=$CC
echo CFLAGS=$CFLAGS
echo LDFLAGS=$LDFLAGS
echo OPTS=$OPTS
echo NVCC=$NVCC
echo CUDAOPT=$CUDAOPT
echo LDCUDA=$LDCUDA
echo LDLIB=$LDLIB
echo LDEXE=$LDEXE
echo LDOPTS=$LDOPTS
echo BIGENDIAN=$BIGENDIAN
echo "*************************************************************************"
echo "*** MAOS Documentation can be found at http://lianqiw.github.com/maos ***"
echo "*** MAOS Source can be obtained at http://github.com/lianqiw/maos     ***"
echo "*************************************************************************"

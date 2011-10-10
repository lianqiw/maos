AC_PREREQ(2.57)
AC_INIT(maos, 0.9.6, lianqiw@gmail.com)
AM_INIT_AUTOMAKE([1.9.6 gnu -Wall no-define])
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])]) #silent building #disable this if automake is too old.
BASEURL=https://github.com/downloads/lianqiw/maos
#If user does not supply any CFLAGS, set something so that AC_PROG_CC won't put a -g -O2 to CFLAGS
if test -z "$CFLAGS" ;then
   CFLAGS=" "
fi
if test -z "$LIBS" ;then
   LIBS=" "
fi
AC_LANG([C])
#After enable HEADER, the package name will not appear in CFLAGS. Appear in config.h instead.
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADER([config.h])
AC_CANONICAL_HOST #defines host host_os host_cpu host_vendor etc

AC_TYPE_SIZE_T  

#Check operating system.
AM_CONDITIONAL(LINUX,  false)
AM_CONDITIONAL(BSD,    false)
AM_CONDITIONAL(DARWIN, false)
AM_CONDITIONAL(WIN,    false)
#top_srcdir and abs_top_srcdir is empty. no idea why. but srcdir sometimes is the full path.
SRC_DIR=$(cd $srcdir && pwd)
BUILD_DIR=$(pwd)
if test -d ${HOME} ;then
LIBDIR=${HOME}/.aos/lib
else
LIBDIR=${SRC_DIR}/external
fi
AC_DEFINE_UNQUOTED(SRCDIR,  ["$SRC_DIR"], ["Define SRC folder"])
AC_DEFINE_UNQUOTED(BUILDDIR,["$BUILD_DIR"],["Define BUILD folder"])
AC_DEFINE_UNQUOTED(BASEURL, ["$BASEURL"], ["The location to download files"])
#Specifying both -O3 and -g in mac produces internal error in gcc 4.2
OPTScommon=" -Wshadow -Wall " #-all-static will produce a static exe.
if test "${USER}" = "lianqiw" ;then
   OPTScommon="$OPTScommon -Werror"
fi
LDSYS="-lz -lm"
OPTSextra=""
use_shm="no"
use_static="no"

if test -n "$(icc -dumpversion)" ;then	
   use_icc_default="yes"
else
   echo "Intel Compiler is not found."
   use_icc_default="no"
fi
if test -n "$(which mex)" ;then
   use_mex_default="yes"
else
   use_mex_default="no"
fi
if test -n "$(nvcc -V)" ;then
   use_cuda_default="yes"
else
   use_cuda_default="no"
fi
AC_ARG_ENABLE(cuda,  [  --enable-cuda   Enable Cuda], [use_cuda="$enableval"], [use_cuda="$use_cuda_default"])
if test "x$use_cuda" = "xyes" ;then
   case "$host_os" in
       *darwin*)
       #In MAC OS, gcc defaults to 64 bit, but nvcc is 32 bit. force both to 64 bit.
       host_cpu=x86_64;
       CFLAGS="$CFLAGS -m64"
       LDFLAGS="$LDFLAGS -m64"
       ;;
   esac
fi
case "$host_os" in  #compile 32 exe in 64 bit.
   *linux*)
   case "$host_cpu" in
      *64)
      AC_ARG_ENABLE(32,   [  --enable-32   Compile 32 bit executables], [use_32="$enableval"], [use_32="no"])
      if test "$use_32" = "yes" ;then
         CFLAGS="$CFLAGS -Xcompiler -m32"
         CHOST="i686-pc-linux-gnu"
         CDEFINE="__i386__"
         HOST="i686-linux"
         BUILD="i686-linux"
         use_icc_default="no"
         host_cpu=i386
         use_mex_default="no" #disable mex
      fi
      ;;
   esac
   ;;
esac
echo host_cpu=$host_cpu
case "$host_cpu" in 
    *86)
        cpu="32"
	AC_DEFINE(X86, [1], "I386 Machine")
	objcopy_output=elf32-i386
	objcopy_arch=i386
	echo "CPU is 32 bit"
	;;
    *64)
        cpu="64"
	AC_DEFINE(X86_64, [1], "X86_64 Machine")
	objcopy_output=elf64-x86-64
	objcopy_arch=i386:x86-64
	echo "CPU is 64 bit"
	;;
    *)
        AC_MSG_ERROR([Unrecognized host_cpu])
esac
AC_SUBST(objcopy_output)
AC_SUBST(objcopy_arch)

mklextra=_redistributed

win32=no
use_mkl="no"
try_mkl="no"
mkl_avail="no"
case "$host_os" in
*linux*)
 AM_CONDITIONAL(LINUX, true)
 LDSYS="$LDSYS -lrt" #real time library.
 OPTSextra="-g -D_BSD_SOURCE -D_POSIX_C_SOURCE=200112L"
 try_mkl="yes"
 use_static_default="no"
 wget="wget --no-check-certificate"
 wgetoutput="-O"
 mkldep="-lguide"
 mklname=mkl${cpu}
 system="linux"
 ldsuffix="so"
 LDFLAGS="$LDFLAGS -Wl,--as-needed"
 ;;
*darwin*)
 AM_CONDITIONAL(DARWIN, true)
 try_mkl="yes"
 use_static_default="no"
 wget="curl -L --progress-bar"
 wgetoutput="--output"
 mkldep="-liomp5"
 mklname=mklfat
 system="apple"
 ldsuffix="dylib"
 ;;
*bsd*)
 AM_CONDITIONAL(BSD, true)
 system="bsd"
 ldsuffix="so"
 ;;
*cygwin*|*mingw*)
 AM_CONDITIONAL(WIN, true)
 LAPACK="-lblas -llapack -lgfortran"
 win32=yes
 use_static_default=yes #blas/lapack are linked statically, but others are still dynamic. ship those
 system="win"
 ldsuffix="dll"
 ;;
*)
AC_MSG_ERROR([Unrecognized OS])
esac
AC_SUBST(ldsuffix)
mkl_lib="${LIBDIR}/${mklname}/libmkl.la"
if test "$try_mkl" = "yes" ;then
    #Enable System V Interface definition (IPC,SHM)
    mkl_avail="yes" #available for download.
    if test -f "$mkl_lib" ;then
       use_mkl_default="yes"
    else
       use_mkl_default="no"
    fi
    if test -d /dev/shm ;then
       use_shm_default="yes"
    else
       use_shm_default="no"
    fi
    #We only have MKL redistributable for linux and MAC.
    AC_ARG_ENABLE(mkl,   [  --enable-mkl   Enable MKL blas/lapack], [use_mkl="$enableval"], [use_mkl="$use_mkl_default"])
    AC_ARG_ENABLE(shm,   [  --disable-shm  Disable posix shared memory usage], [use_shm="$enableval"],[use_shm="$use_shm_default"])
fi

#Common configure options.
AC_ARG_ENABLE(all-static,[  --diable-all-static Compile statically], [use_static="$enableval"],[use_static=$use_static_default])
AC_ARG_ENABLE(debug, [  --enable-debug  Turn on debug], [debug="$enableval"], [debug="no"])
AC_ARG_ENABLE(mem,   [  --enable-mem    Turn on memory management],[use_mem="$enableval"],[use_mem="no"])
AC_ARG_ENABLE(icc,   [  --enable-icc    Using intel compiler], [use_icc="$enableval"], [use_icc=$use_icc_default])
AC_ARG_ENABLE(oprof, [  --enable-oprof  Enable oprof for code profiling.],[use_oprof="$enableval"],[use_oprf="no"])
AC_ARG_ENABLE(opencc,[  --enable-opencc Enable opencc compiling], [use_opencc="$enableval"], [use_opencc="no"])
AC_ARG_ENABLE(long,  [  --enable-long   Use long for sparse matrix index], [use_long="$enableval"], [use_long="yes"])
AC_ARG_ENABLE(duma,  [  --enable-duma   Use duma, the memory debugger],[use_duma="$enableval"],[use_duma="no"])
AC_ARG_ENABLE(mex,   [  --enable-mex    Enable compiling of matlab mex routines, requires mex to be in the PATH],[use_mex="$enableval"],[use_mex="$use_mex_default"])
if test "x$use_static" = "xyes" ;then
   LDEXE="-all-static" #Enable all static build.
   AM_CONDITIONAL(USE_STATIC, true)
   AC_DEFINE(USE_STATIC, [1], "We are compling maos statically")
   #(cd ${SRC_DIR} && tar zcf ${BUILD_DIR}/config.tar.gz config/)
else
   LDEXE=""
   AM_CONDITIONAL(USE_STATIC, false)
   AC_DEFINE(USE_STATIC, [0], "We are compling maos dynamically")
fi

if test "x$use_shm" = "xyes" ;then
   AC_DEFINE(USE_POSIX_SHM, [1], "Enable posix shared memory")
else
   AC_DEFINE(USE_POSIX_SHM, [0], "Disable posix shared memory")
fi 
AC_CHECK_LIB([z],[gzopen],[],[AC_MSG_ERROR([zlib is not found])])
AC_CHECK_LIB([pthread],[pthread_create],[pthread_good="yes"],[pthread_good="no"])
AC_DEFINE(_THREAD_SAFE,[1], ["Enable thread safe"])
if test "x$pthread_good" = "xyes" ;then
    AC_DEFINE(USE_PTHREAD,[1], ["Enable thread"])
    if test "x$win32" != "xyes" ;then
    echo        OPTSprof="-pthread"
    fi
else
    AC_DEFINE(USE_PTHREAD,[0], ["Disable thread"])
fi
if test "x$use_duma" = "xyes" ;then	
    AC_CHECK_LIB([duma],[malloc],[],[AC_MSG_ERROR([duma not found])])
    LDSYS="$LDSYS -lduma"
fi 

if test x$use_oprof = "xyes" ;then
   OPTSprof="$OPTSprof -g -fno-omit-frame-pointer -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls"
fi
#Set optimize flags for different compilers
#-xO Enables SSE3, SSE2 and SSE instruction sets optimizations for non-Intel CPUs 
#-wd1292 : disables warning "error #1292: attribute "__alloc_size__" ignored"
COMMON_ICC="-wd981 -wd1292 -vec-report=0 -diag-disable remark"
OPTIM_ICC="-O3 -funroll-loops -DHAVE_SSE2=1 -static-intel"
OPTIM_GCC="-O3 -funroll-loops -DHAVE_SSE2=1 -msse2 -msse -fprefetch-loop-arrays -fno-math-errno -ffinite-math-only -static-libgcc"
OPTIM_OPENCC="-O3"
OPT_DEBUG="-g -rdynamic -ggdb -O0 -DDEBUG_TRACE"
AM_CONDITIONAL(USE_ICC,false)
AM_CONDITIONAL(USE_OPENCC,false)
use_ifort="no"
if test "x$use_mem" = "xyes" ;then
    AC_DEFINE(USE_MEM,[1],["Use memory management"])
    use_icc="no" #disable icc.
    debug="yes"
else
    AC_DEFINE(USE_MEM,[0],["Do not use memory management"])	
fi
if test "x$use_icc" = "xyes" ;then
    OPTIM="$OPTIM_ICC $COMMON_ICC"
    OPT_DEBUG="$OPT_DEBUG $COMMON_ICC"
    CC=icc
    AC_PROG_CC([icc])
    which icc || AC_MSG_ERROR([icc is not found in path.])
    AM_CONDITIONAL(USE_ICC,true)	
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
elif test "x$use_opencc" = "xyes" ;then
    OPTIM=$OPTIM_OPENCC
    CC=opencc
    AC_PROG_CC([opencc])
    AM_CONDITIONAL(USE_OPENCC,true)
else
    OPTIM=$OPTIM_GCC
    AC_PROG_CC

    #Test GCC version
    CC_VERSION=$($CC -dumpversion)
    if test -n $CC_VERSION ;then
	major=${CC_VERSION%%.*}
	minor=${CC_VERSION%.*}
	minor="${minor#4.}"
	echo major=$major minor=$minor
	if test $major -lt 4 ;then
	    AC_MSG_ERROR([At least version 4 is required when compiling with gcc.])
	fi
	if test $major -ge 4 -a $minor -ge 2;then
	    echo "Enabling -mtune=native"
	    OPTIM="$OPTIM -mtune=native" #introduced in GCC 4.2
	    OPTIM="$OPTIM -ftree-vectorize " #introduced in GCC 4
	fi
    fi
fi	
if test "x$debug" = "xyes" ;then
    AM_CONDITIONAL(DEBUG,true)
    AC_DEFINE(DEBUG,[1],["Debug the code"])
    OPTScc="$OPT_DEBUG"
else
    AM_CONDITIONAL(DEBUG,false)
    AC_DEFINE(DEBUG,[0],["Do not debug"])
    OPTScc="$OPTIM"
fi
AC_PROG_CC_STDC
AM_PROG_CC_C_O
AC_PROG_CC_C99 #enable C99.

#call convention of fortran libraries
#AC_F77_FUNC(dgemm)
#AC_F77_FUNC(dposv)
function download(){
   echo The $1 library is not available in the system
   echo Will download one from the website.
   echo This library is only for non-commercial use.
   fn=$2
   $wget ${BASEURL}/${fn} $wgetoutput /tmp/${fn} || exit 1
   if test ! -d ${LIBDIR} ;then
      mkdir -p ${LIBDIR} || exit 1
   fi
   tar -xjvf /tmp/${fn} -C ${LIBDIR}  || exit 1
   rm -rf /tmp/${fn}
}
if test ! "x$use_mkl" = "xyes" ;then
#Check for LAPACK library in the system
   if test -z "$LAPACK" ;then
       if test -n "$(pkg-config --libs blas lapack)" ;then
           LAPACK=$(pkg-config --libs blas lapack)
       else
           LAPACK="-lblas -llapack"
       fi
   fi
       echo Trying LAPACK=$LAPACK
       AC_CHECK_LIB([blas],  [dgemm_],[],[has_lapack=no],[$LAPACK])
       AC_CHECK_LIB([lapack],[dposv_],[],[has_lapack=no],[$LAPACK])
       if test  "x$has_lapack" = "xno" ;then
           if test "x$mkl_avail" = "xyes"; then
       	       use_mkl="yes"
	       AC_MSG_NOTICE([blas/lapack is not found. Download and use mkl])
           else
               AC_MSG_ERROR([Please install blas and lapack])
	   fi
       fi
fi
if test "x$use_mkl" = "xyes" ;then #Use shipped libmkl64.so. only for linux.
   if test ! -f "$mkl_lib"  ;then
      download MKL ${mklname}${mklextra}.tar.bz2
   fi
   LAPACK="-L${LIBDIR}/${mklname} -lmkl -pthread ${mkldep}"
   AC_DEFINE(USE_MKL,[1],["MKL is used"])
else
   AC_DEFINE(USE_MKL,[0],["MKL is not used"])
fi
#Check for FFTW3 header and library.
has_fftw='yes'
AC_CHECK_HEADERS([fftw3.h],[], [has_fftw='no'])
AC_CHECK_LIB([fftw3], [fftw_execute],[],[has_fftw='no'])
if test "$has_fftw" = "no" ;then #fftw does not exist in system.
   if test "$system" = "linux" ;then
      if test ! -d ${LIBDIR}/fftw${cpu};then 
         download FFTW fftw${cpu}.tar.bz2
      fi
      LDFFT="-L${LIBDIR}/fftw${cpu}"
      OPTScommon="$OPTScommon -I${LIBDIR}/fftw${cpu}"
   else
      AC_MSG_ERROR([Please install fftw library])
   fi
else
   LDFFT=""
fi

AC_CHECK_LIB([fftw3_threads], [fftw_init_threads], [has_fftw_threads="yes"],[has_fftw_threads="no"],[$LDFFT])
AC_CHECK_LIB([fftw3_threads], [omp_get_num_threads], [has_fftw_omp="yes"], [has_fftw_omp="no"],[$LDFFT])
if test "x$has_fftw_threads" = "xyes" ;then
    if test "x$has_fftw_omp" = "xyes" ;then
        AC_DEFINE(USE_FFTW_THREADS, [2], "We have FFTW3 threads support with openmp, conflict.")
	LDFFT="$LDFFT -lfftw3" 
    else
        AC_DEFINE(USE_FFTW_THREADS, [1], "We have FFTW3 threads support")
        LDFFT="$LDFFT -lfftw3_threads -lfftw3"
    fi
else
    LDFFT="$LDFFT -lfftw3"
    AC_DEFINE(USE_FFTW_THREADS, [0], "We do not have FFTW3 threads")
fi
#Check availability of gtk and libnotify. If available, we will build drawdaemon and the monitor.
AM_CONDITIONAL(DRAW,true) #for drawdaemon
AM_CONDITIONAL(MON,true)  #for parital monitor without notify
AM_CONDITIONAL(MONN,true) #for full monitor with notify
#make sure there are spaces before and after >=
PKG_CHECK_MODULES(DRAW, gtk+-2.0 >= 2.2 cairo >= 1.4,    [comp_draw="yes"],[comp_draw="no"])
PKG_CHECK_MODULES(MON,  gtk+-2.0 >= 2.2                 ,[comp_mon="yes"], [comp_mon="no"])
PKG_CHECK_MODULES(MONN, gtk+-2.0 >= 2.2 libnotify >= 0.1,[comp_monn="yes"],[comp_monn="no"])
#OLD GLIB has GTHREAD separated
PKG_CHECK_MODULES(GTHREAD, gthread-2.0 >= 2.0, [has_gthread="yes"], [has_gthread="no"])
#The above will define MON_CFLAGS and MON_LIBS. The following will substitute the variables name 
#with their values in Makefile.in 
#so that you can use $(DRAW_CFLAGS) in Makefile.am to represent the values of the shell variable here.
AM_CONDITIONAL(DRAW, test x$comp_draw = xyes)
AM_CONDITIONAL(MON,  test x$comp_mon  = xyes)
AM_CONDITIONAL(MONN, test x$comp_monn = xyes)
if test -n "$GTHREAD_LIBS" ;then
   DRAW_LIBS="$DRAW_LIBS $GTHREAD_LIBS"
   MON_LIBS="$MON_LIBS $GTHREAD_LIBS"
   MONN_LIBS="$MONN_LIBS $GTHREAD_LIBS"
fi
AC_SUBST(DRAW_CFLAGS)
AC_SUBST(DRAW_LIBS)
AC_SUBST(MON_CFLAGS)
AC_SUBST(MON_LIBS)
AC_SUBST(MONN_LIBS)
AC_SUBST(MONN_CFLAGS)

#Simply renaming the compiled cholmod to cholmod_long and cholmod_int does not help after linking
if test "x$use_long" = "xyes" -a "x$cpu" = "x64" ;then
   OPTScommon="$OPTScommon -DDLONG"
   SUFFIX="long"
   chol_fun=cholmod_l_factorize
   AM_CONDITIONAL(USE_LONG, true)
else
   SUFFIX="int"
   chol_fun=cholmod_factorize
   AM_CONDITIONAL(USE_LONG, false)
fi

AC_CHECK_LIB([cholmod_${SUFFIX}],[${chol_fun}], [have_chol="yes"], [have_chol="no"],[$LAPACK -pthread])
if test "x$have_chol" = "xyes" ;then
    LDCHOL="-lcholmod_${SUFFIX}"
    AM_CONDITIONAL(COMP_CHOL, false)
else
    LDCHOL="../external/cholmod/libcholmod_${SUFFIX}.la"
    AM_CONDITIONAL(COMP_CHOL, true)
fi
# Here cholmod_long and cholmod_int are custom made cholmod libraries. We do not
# use the cholmod in the system even if it is available because we want to
# customize the dependency on blas/lapack.

AC_CHECK_HEADERS([unistd.h])
#LDCHOL must be before LAPACK
LDOPTS="$LAPACK $LDSYS $LDFFT -lm"
LDEXE="${LDEXE} -no-install"
OPTS="$OPTScommon $OPTScc $OPTSprof $OPTSextra "
AC_SUBST(LDAOS)
AC_SUBST(OPTS)
AC_SUBST(LDOPTS)
AC_SUBST(LDEXE)
AC_SUBST(LDSYS)
AC_SUBST(LDCHOL)
AC_SUBST(LDFFT)
AC_SUBST(LAPACK)
AC_SUBST(OPTScommon)
AC_SUBST(SRC_DIR)
AC_SUBST(BUILD_DIR)
LIBS="" #somehow LIBS contain -fftw3
AC_SUBST(LIBS)
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

#The following must be after AC_PROG_CC
AC_PROG_MAKE_SET
AC_PROG_INSTALL
#if test "x${use_mex}" = "xyes" ; then
AC_ENABLE_SHARED
#else
#AC_DISABLE_SHARED
#fi
#if test "x${use_static}" = "xyes" ;then
AC_ENABLE_STATIC
#else
#AC_DISABLE_STATIC
#fi
LT_INIT([pic-only])
AC_PROG_LIBTOOL

AC_CONFIG_FILES( [Makefile
		 lib/Makefile
		 tools/Makefile
		 maos/Makefile 
		 skyc/Makefile
		 ])
if test -d ${srcdir}/test ;then
AC_CONFIG_FILES([test/Makefile])
fi

if test -d ${srcdir}/external/cholmod ;then
AC_CONFIG_FILES([external/cholmod/Makefile])
fi

if test "x${use_mex}" = "xyes" ; then
   AM_CONDITIONAL(USE_MEX, true)
   AC_CONFIG_FILES([mex/Makefile])
   AC_ARG_WITH(matlab, --with-matlab=DIR, [MEXDIR=${withval}], [MEXDIR=])
   if test -n "$MEXDIR" ;then
       if test -f ${MEXDIR}/mex ;then
           MEXDIR=${MEXDIR}
       elif test -f ${MEXDIR}/bin/mex ;then
           MEXDIR=${MEXDIR}/bin
       else
           AC_MSG_ERROR([mex not found in specified matlab path ${withval}])
       fi
   else
       MEX=$(which mex)
       if test -z "${MEX}" ;then
           AC_MSG_ERROR([mex not found. please reconfigure with --with-matlab=DIR to specify MATLAB location])
       elif test -L ${MEX} ;then
           MEX=$(readlink $MEX)
       fi
       MEXDIR=$(dirname $MEX)
       if test ! -d $MEXDIR ;then
           AC_MSG_ERROR([$MEXDIR does not exist. reconfigure with --disable-mex to bypass mex.])
       fi
   fi
   MEXEXT=$($MEXDIR/mexext)
   if test -z "$MEXEXT" ;then
       AC_MSG_ERROR([Unable to determine mex extension. reconfigure with --disable-mex to bypass mex.])
   fi
   case $MEXEXT in
   	mexglx)
		MEXLIB=glnx86
		MEXOPT="-O3 -DNDEBUG -fPIC -pthread -m32 -fexceptions -D_FILE_OFFSET_BITS=64"
		MEXOPT="${MEXOPT} -D_BSD_SOURCE -D_POSIX_C_SOURCE=200112L"
	;;
   	mexa64)
		MEXLIB=glnxa64
		MEXOPT="-O3 -DNDEBUG -fPIC -pthread -m64 -fexceptions -fno-omit-frame-pointer"
		MEXOPT="${MEXOPT} -DDLONG -D_BSD_SOURCE -D_POSIX_C_SOURCE=200112L"
	;;
	mexmaci)
		MEXLIB=maci
		MEXOPT="-O3 -DNDEBUG -fno-common -no-cpp-precomp -m32 -fexceptions"
	;;
	mexmaci64)
		MEXLIB=maci64
		MEXOPT="-O3 -DNDEBUG -fno-common -no-cpp-precomp -fexceptions -arch x86_64"
		MEXOPT="${MEXOPT} -DDLONG"
	;;
	mexw64)
		MEXLIB=win64
		MEXOPT=""
	;;
	mexw32)
		MEXLIB=win32
		MEXOPT=""
	;;
	*)
		AC_MSG_ERROR([Unknown mex extension "$MEXEXT"])
	;;
   esac
   LDMEX="-L${MEXDIR}/${MEXLIB} -Wl,-rpath,${MEXDIR}/${MEXLIB} -lmx -lmex -lmat -lstdc++ -lz"
   if test x$use_icc = xyes ;then
       LDMEX="${LDMEX} -static-intel"
   fi
   MEXOPT="${MEXOPT} -I${MEXDIR}/../extern/include -DMATLAB_MEX_FILE"
   AC_SUBST(MEXEXT)
   AC_SUBST(LDMEX)
   AC_SUBST(MEXOPT)
   CFLAGS="$CFLAGS"
else
   AM_CONDITIONAL(USE_MEX, false)
fi
if test "$use_cuda" = "yes" ;then
   if test "$use_static" = "yes"; then
      AC_MSG_ERROR("all-static conflicts with cuda support. Please disable all-static")
   fi
   AC_CONFIG_FILES([cuda/Makefile])
   AM_CONDITIONAL(USE_CUDA, true)
   AC_DEFINE(USE_CUDA, [1], "Has Cuda")
   AC_MSG_NOTICE([Enabling CUDA])
   AC_ARG_WITH(cuda, --with-cuda=DIR, [CUDADIR=${withval}], [CUDADIR=])
   LDCUDA=
   if test -n "$CUDADIR" ;then
       AC_MSG_NOTICE([User provided cuda directory $CUDADIR])
       NVCC=$CUDADIR/bin/nvcc
       if test ! -x $NVCC ;then
       	  NVCC=$CUDADIR/nvcc
	  CUDADIR=$(dirname $CUDADIR)
       fi
   else
       if test -n "$(which nvcc)"  ;then
           NVCC=$(which nvcc)
           if test -L $NVCC ;then
               NVCC=$(readlink $NVCC)
           fi
       elif test -f /usr/local/cuda/bin/nvcc ;then
           NVCC=/usr/local/cuda/bin/nvcc
       elif test -f /opt/cuda/bin/nvcc ;then
           NVCC=/opt/cuda/bin/nvcc
       fi
       CUDADIR=$(dirname $(dirname ${NVCC}))
   fi
   AC_MSG_NOTICE([CUDADIR=$CUDADIR])
   if test -f ${CUDADIR}/lib/libcudart.$ldsuffix ;then
       CUDA_L=${CUDADIR}/lib
   elif test -f ${CUDADIR}/lib${cpu}/libcudart.$ldsuffix ;then
       CUDA_L=${CUDADIR}/lib${cpu}
   fi
   if test -d "$CUDA_L" ;then
       LDCUDA="-L$CUDA_L -R$CUDA_L"
   else
       LDCUDA=
   fi

   CUDA_H=$(dirname "$(find $CUDADIR -name 'cuda.h' || find $CUDADIR/../ -name 'cuda.h')")
   if test -d $CUDA_H ;then
       CUDAOPT="-I$CUDA_H"
       if test ! -f $CUDA_H/cuda.h ;then
           AC_MSG_ERROR([Please use --with-cuda=DIR to specify CUDA location])
       fi
   else
       AC_CHECK_HEADERS(cuda.h, [] ,[AC_MSG_ERROR([Please use --with-cuda=DIR to specify CUDA location])])
       CUDAOPT=
   fi
   echo LDCUDA=$LDCUDA
   echo CUDAOPT=$CUDAOPT
   if test ! -x $NVCC ;then
        AC_MSG_ERROR([nvcc not found])
   fi
   AC_CHECK_LIB([cudart], [cudaSetDevice], [cuda_lib="yes"], [cuda_lib="no"], [$LDCUDA])
   if test x$cuda_lib = xno ;then
       AC_MSG_ERROR([libcudart not found. Please use --with-cuda=DIR to specify CUDA location])
   fi
   LDCUDA="$LDCUDA -lcudart -lcurand -lcusparse -lcufft -lcublas -lstdc++"
   AC_SUBST(CUDAOPT)
   AC_SUBST(LDCUDA)
   AC_SUBST(NVCC)
else
   AM_CONDITIONAL(USE_CUDA, false)
   AC_DEFINE(USE_CUDA, [0], "Don't have Cuda")
fi
AC_OUTPUT
echo SRC_DIR=$SRC_DIR
echo BUILD_DIR=$BUILD_DIR
echo host_os=$host_os
echo use_mex=$use_mex
echo use_icc=$use_icc
echo use_cuda=$use_cuda
echo use_static=$use_static
echo Compile Drawdaemon=$comp_draw
echo Compile Monitor=$comp_mon
echo Compile Monitor with Libnotify=$comp_monn
echo CC=$CC
echo CFLAGS=$CFLAGS
echo LDFLAGS=$LCFLAGS
echo LIBS=$LIBS
echo LDEXE=$LDEXE
echo OPTS=$OPTS
echo LDOPTS=$LDOPTS
echo LAPACK=$LAPACK
echo LDCHOL=$LDCHOL
echo LDFFT=$LDFFT
if test x$use_icc_default = xno -a x$system = xlinux ;then
echo Intel icc is not enabled. Free icc for noncommerical use can be obtain from Intel website
fi
echo "*********************************************************************************"
echo "*** MAOS Documentation can be found at http://lianqiw.github.com/maos        ****"
echo "*** MAOS Source can be obtained at http://github.com/lianqiw/maos/downloads  ****"
if test ! -d "${HOME}/.aos" ;then
echo "*** Please drop me a line to lianqiw@gmail.com if you find MAOS useful       ****"
fi
echo "*********************************************************************************"

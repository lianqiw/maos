AC_PREREQ(2.57)
AC_INIT(maos, 0.8.0, lianqiw@gmail.com)
AM_INIT_AUTOMAKE([1.9.6 gnu -Wall no-define])
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])]) #silent building #disable this if automake is too old.
BASEURL=http://github.com/downloads/lianqiw/maos
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
#Check availability of gtk and libnotify. If available, we will build drawdaemon and the monitor.
AM_CONDITIONAL(DRAW,true)
AM_CONDITIONAL(MON,true)  #for parital monitor without notify
AM_CONDITIONAL(MONN,true) #for full monitor with notify
#make sure there are spaces before and after >=
PKG_CHECK_MODULES(DRAW, gtk+-2.0 >= 2.2 cairo >= 1.4,    [],[AM_CONDITIONAL(DRAW,false)])
PKG_CHECK_MODULES(MON,  gtk+-2.0 >= 2.2                 ,[],[AM_CONDITIONAL(MON,false)])
PKG_CHECK_MODULES(MONN, gtk+-2.0 >= 2.2 libnotify >= 0.1,[],[AM_CONDITIONAL(MONN,false)])
#The above will define MON_CFLAGS and MON_LIBS. The following will substitute the variables name 
#with their values in Makefile.in 
#so that you can use $(DRAW_CFLAGS) in Makefile.am to represent the values of the shell variable here.
AC_SUBST(DRAW_CFLAGS)
AC_SUBST(DRAW_LIBS)
AC_SUBST(MON_CFLAGS)
AC_SUBST(MON_LIBS)
AC_SUBST(MONN_LIBS)
AC_SUBST(MONN_CFLAGS)
#Check operating system.
AM_CONDITIONAL(LINUX,  false)
AM_CONDITIONAL(BSD,    false)
AM_CONDITIONAL(DARWIN, false)
AM_CONDITIONAL(WIN,    false)
#top_srcdir and abs_top_srcdir is empty. no idea why. but srcdir sometimes is the full path.
SRC_DIR=$(cd $srcdir && pwd)
SVN_REV=$(svnversion -n $SRC_DIR)
BUILD_DIR=$(pwd)
AC_DEFINE_UNQUOTED(SVN_REV, ["$SVN_REV"], ["Define SVN revision number"])
AC_DEFINE_UNQUOTED(SRCDIR,  ["$SRC_DIR"], ["Define SRC folder"])
AC_DEFINE_UNQUOTED(BUILDDIR,["$BUILD_DIR"],["Define BUILD folder"])
AC_DEFINE_UNQUOTED(BASEURL, ["$BASEURL"], ["The location to download files"])
OPTScommon=" -Wshadow  -Wall " #-all-static will produce a static exe.
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
case "$host_cpu" in 
    *86)
        cpu="32";;
    *64)
        cpu="64";;
    *)
        AC_MSG_ERROR([Unrecognized host_cpu])
esac

mkl_lib="${SRC_DIR}/external/mkl${cpu}/libmkl.la"
win32=no
case "$host_os" in
*linux*)
    AM_CONDITIONAL(LINUX, true)
    LDSYS="$LDSYS -lrt" #Real time library
    #Enable System V Interface definition (IPC,SHM)
    OPTSextra="-Werror -g -D_XOPEN_SOURCE=600 -D_BSD_SOURCE"
    if test -f "$mkl_lib" ;then
       use_mkl_default="yes"
    else
       use_mkl_default="yes" #"no"
    fi
    use_static_default="no"
    if test -d /dev/shm ;then
       use_shm_default="yes"
    else
       use_shm_default="no"
    fi
    #We only have MKL redistributable for linux
    AC_ARG_ENABLE(mkl,   [  --enable-mkl   Enable MKL blas/lapack], [use_mkl="$enableval"], [use_mkl="$use_mkl_default"])
    AC_ARG_ENABLE(shm,   [  --disable-shm  Disable posix shared memory usage], [use_shm="$enableval"],[use_shm="$use_shm_default"])
    ;;
*darwin*)
 AM_CONDITIONAL(DARWIN, true)
 if ! test -e ${SRC_DIR}/scripts/drawdaemon.app ;then
     curl -O -L --progress-bar https://github.com/downloads/lianqiw/maos/drawdaemon.app.tar.gz &&\
     tar -xzf drawdaemon.app.tar.gz -C ${SRC_DIR}/scripts/ &&\
     rm  -rf drawdaemon.app.tar.gz
 fi
 if ! test -e ${SRC_DIR}/scripts/monitor.app ;then
     curl -O -L --progress-bar https://github.com/downloads/lianqiw/maos/monitor.app.tar.gz &&\ 
     tar -xzf monitor.app.tar.gz -C ${SRC_DIR}/scripts/ &&\
     rm	 -rf monitor.app.tar.gz
 fi
 echo "Please use the monitor downloaded to ${SRC_DIR}/scripts/monitor.app"
;;
*bsd*)
 AM_CONDITIONAL(BSD, true)
;;
*cygwin*)
 AM_CONDITIONAL(WIN, true)
 LAPACK="-lblas -llapack -lgfortran"
 win32=yes
 use_static_default=yes #blas/lapack are linked statically, but others are still dynamic. ship those
;;
*)
AC_MSG_ERROR([Unrecognized OS])
esac
AC_ARG_ENABLE(all-static,[  --diable-all-static Compile statically], [use_static="$enableval"],[use_static=$use_static_default])

#Common configure options.
AC_ARG_ENABLE(debug, [  --enable-debug  Turn on debug], [debug="$enableval"], [debug="no"])
AC_ARG_ENABLE(mem,   [  --enable-mem    Turn on memory management],[use_mem="$enableval"],[use_mem="no"])
AC_ARG_ENABLE(icc,   [  --enable-icc    Using intel compiler], [use_icc="$enableval"], [use_icc=$use_icc_default])
AC_ARG_ENABLE(oprof, [  --enable-oprof  Enable oprof for code profiling.],[use_oprof="$enableval"],[use_oprf="no"])
AC_ARG_ENABLE(opencc,[  --enable-opencc Enable opencc compiling], [use_opencc="$enableval"], [use_opencc="no"])
AC_ARG_ENABLE(long,  [  --enable-long   Use long for sparse matrix index], [use_long="$enableval"], [use_long="no"])
AC_ARG_ENABLE(duma,  [  --enable-duma   Use duma, the memory debugger],[use_duma="$enableval"],[use_duma="no"])
if test "x$use_long" = "xyes" ;then
   OPTScommon="$OPTScommon -DDLONG"
   AM_CONDITIONAL(USE_LONG, true)
else
   AM_CONDITIONAL(USE_LONG, false)
fi
if test "x$use_static" = "xyes" ;then
   LDEXE="-all-static" #Enable all static build.
   AM_CONDITIONAL(USE_STATIC, true)
else
   LDEXE=""
   AM_CONDITIONAL(USE_STATIC, false)
fi

if test "x$use_shm" = "xyes" ;then
   AC_DEFINE(USE_POSIX_SHM, [1], "Enable posix shared memory")
else
   AC_DEFINE(USE_POSIX_SHM, [0], "Disable posix shared memory")
fi 

AC_CHECK_LIB([pthread],[pthread_create],[pthread_good="yes"],[pthread_good="no"])
AC_DEFINE(_THREAD_SAFE,[1], ["Enable thread safe"])
if test "x$pthread_good" = "xyes" ;then
    AC_DEFINE(USE_PTHREAD,[1], ["Enable thread"])
    OPTSprof="-pthread"
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
OPTIM_ICC="-O3 -xO -funroll-loops -DHAVE_SSE2=1 -wd981 -vec-report=0 -static-intel"
OPTIM_GCC="-O3 -funroll-loops -DHAVE_SSE2=1 -msse2 -msse -fprefetch-loop-arrays -fno-math-errno -ffinite-math-only -static-libgcc"
OPTIM_OPENCC="-O3"
OPT_DEBUG="-g -rdynamic -ggdb -O0 -DDEBUG_TRACE"
AM_CONDITIONAL(USE_ICC,false)
AM_CONDITIONAL(USE_OPENCC,false)
use_ifort="no"
if test "x$use_mem" = "xyes" ;then
    AC_DEFINE(USE_MEM,[1],["Use memory management"])
else
    AC_DEFINE(USE_MEM,[0],["Do not use memory management"])	
fi
if test "x$use_icc" = "xyes" ;then
    OPTIM=$OPTIM_ICC
    OPT_DEBUG="$OPT_DEBUG -wd981"
    CC=icc
    AC_PROG_CC([icc])
    which icc || AC_MSG_ERROR([icc is not found in path.])
    AM_CONDITIONAL(USE_ICC,true)	
elif test "x$use_opencc" = "xyes" ;then
    OPTIM=$OPTIM_OPENCC
    CC=opencc
    AC_PROG_CC([opencc])
    AM_CONDITIONAL(USE_OPENCC,true)
else
    OPTIM=$OPTIM_GCC
    AC_PROG_CC
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


if test "x$use_icc" = "xno" ;then
    #Test GCC version
    GCC_VERSION=$($CC -dumpversion)
    if test -n $GCC_VERSION ;then
	major=${GCC_VERSION%%.*}
	minor=${GCC_VERSION%.*}
	minor="${minor#4.}"
	echo major=$major minor=$minor
	if test $major -lt 4 ;then
	    AC_MSG_ERROR([At least version 4 is required when compiling with gcc.])
	fi
	if test $major -ge 4 -a $minor -ge 2 ;then
	    echo "Enabling -mtune=native"
	    OPTScc="$OPTScc -mtune=native" #introduced in GCC 4.2
	    OPTScc="$OPTScc -ftree-vectorize " #introduced in GCC 4
	fi
    fi
fi
#call convention of fortran libraries
#AC_F77_FUNC(dgemm)
#AC_F77_FUNC(dposv)
if test "x$use_mkl" = "xyes" ;then #Use shipped libmkl64.so. only for linux.
   if test ! -f "$mkl_lib"  ;then
      echo The MKL library is not available at "$mkl_lib"
      echo Will download one form the website. 
      echo This library is only for non-commerical use only.
      fn=mkl${cpu}_redistributed.tar.bz2
      wget --no-check-certificate ${BASEURL}/${fn}  -O /tmp/${fn}  || exit 1
      if test ! -d ${SRC_DIR}/external ;then
      	  mkdir -p ${SRC_DIR}/external || exit 1
      fi
      tar -xjvf /tmp/${fn} -C ${SRC_DIR}/external  || exit 1
      rm -rf /tmp/${fn}
   fi
   LAPACK="-L${SRC_DIR}/external/mkl${cpu} -lmkl -lguide"
   AC_DEFINE(USE_MKL,[1],["MKL is used"])
else
   if test -z "$LAPACK" ;then
       if test -n "$(pkg-config --libs blas lapack)" ;then
           LAPACK=$(pkg-config --libs blas lapack)
       else
           LAPACK="-lblas -llapack"
       fi
   
   echo LAPACK=$LAPACK
   AC_CHECK_LIB([blas],  [dgemm_],[],[lapack_not_found=yes],[$LAPACK])
   AC_CHECK_LIB([lapack],[dposv_],[],[lapack_not_found=yes],[$LAPACK])
   if test  "x$lapack_not_found" = "xyes" ;then
       if test -n "$use_mkl_default" ;then
           AC_MSG_ERROR(Please install blas or use --enable-mkl to download MKL from MAOS website)
       fi
   fi
   fi
   AC_DEFINE(USE_MKL,[0],["MKL is not used"])
fi
#Check for FFTW3 header and library.
AC_CHECK_HEADERS([fftw3.h],[],AC_MSG_ERROR([Please install fftw3 library]))
AC_CHECK_LIB([fftw3], [fftw_execute],[],[AC_MSG_ERROR(Please install fftw3)])
LDFFT="-lfftw3"

# Here cholmod_long and cholmod_int are custom made cholmod libraries. We do not
# use the cholmod in the system even if it is available because we want to
# customize the dependency on blas/lapack.

AC_CHECK_LIB([cholmod_long],[cholmod_l_factorize], [have_chol_long="yes"], [have_chol_long="no"],[$LAPACK -lpthread])
AC_CHECK_LIB([cholmod_int],[cholmod_factorize], [have_chol_int="yes"], [have_chol_int="no"],[$LAPACK -lpthread])
if test "x$use_long" = "xyes" -a "x$have_chol_long" = "xyes" ;then
    echo "cholmod_long is installed in the system. will use it."
    LDCHOL="-lcholmod_long"
    AM_CONDITIONAL(COMP_CHOL, false)
elif test "x$use_long" = "xno" -a "x$have_chol_int" = "xyes"  ;then
    echo "cholmod_int is installed in the system. will use it."
    LDCHOL="-lcholmod_int"
    AM_CONDITIONAL(COMP_CHOL, false)
else
    echo "cholmod is not available, will use the shipped source."
    LDCHOL="../external/cholmod/libcholmod.la"
    AM_CONDITIONAL(COMP_CHOL, true)
fi
AC_CHECK_HEADERS([unistd.h])
#LDCHOL must be before LAPACK
LDOPTS="$LAPACK $LDSYS"
OPTS="$OPTScommon $OPTScc $OPTSprof $OPTSextra "
AC_SUBST(OPTS)
AC_SUBST(LDOPTS)
AC_SUBST(LDEXE)
AC_SUBST(LDSYS)
AC_SUBST(LDCHOL)
AC_SUBST(LDFFT)
AC_SUBST(LAPACK)
AC_SUBST(OPTScommon)
AC_SUBST(SRC_DIR)
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
AC_DISABLE_SHARED #disable shared object build on LTLIBARIRIES. should be before AC_PROG_LIBTOOL
AC_ENABLE_STATIC
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

AC_OUTPUT
echo host_os=$host_os
echo use_icc=$use_icc
echo use_static=$use_static

echo CC=$CC
echo LIBS=$LIBS
echo LDEXE=$LDEXE
echo SRC_DIR=$SRC_DIR
echo SVN_REV=$SVN_REV
echo BUILD_DIR=$BUILD_DIR
echo OPTS=$OPTS

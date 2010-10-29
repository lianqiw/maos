AC_PREREQ(2.57)
AC_INIT(maos, 0.6.2, lianqiw@gmail.com)
AM_INIT_AUTOMAKE
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])]) #silent building #disable this if automake is too old.
BASEURL=http://github.com/downloads/lianqiw/maos
####NOTES####
#USE F77 Instead of FC. FC is for fortran90/95 compiling.
#############
#If user does not supply any CFLAGS, set something so that AC_PROG_CC won't put a -g -O2 to CFLAGS
if test -z "$CFLAGS" ;then
   CFLAGS=" "
fi
#For CPP. not really used here. 
if test -z "$CXXFLAGS" ;then
   CXXFLAGS=" "
fi
#for F77.  don't want warning as it is external code
if test -z "$FFLAGS" ;then
   FFLAGS=" "
fi
#for f90 not really used here
if test -z "$FCFLAGS" ;then
   FCFLAGS=" "
fi
if test -z "$LIBS" ;then
   LIBS=" "
fi

AC_LANG([C])
#After enable HEADER, the package name will not appear in CFLAGS. Appear in config.h instead.
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADER([config.h])
AC_CANONICAL_HOST #defines host host_os host_cpu host_vendor etc

#echo host=$host
#echo host_alias=$host_alias
#echo host_cpu=$host_cpu
#echo host_vendor=$host_vendor
#echo host_os=$host_os


AC_TYPE_SIZE_T  
#Check availability of gtk and libnotify. If available, we will build drawdaemon and the monitor.
AM_CONDITIONAL(DRAW,true)
AM_CONDITIONAL(MON,true)
#make sure there are spaces before and after >=
PKG_CHECK_MODULES(DRAW, gtk+-2.0 >= 2.2 cairo >= 1.4,    [],[AM_CONDITIONAL(DRAW,false)])
PKG_CHECK_MODULES(MON,  gtk+-2.0 >= 2.2 libnotify >= 0.1,[],[AM_CONDITIONAL(MON,false)])
#The above will define MON_CFLAGS and MON_LIBS. The following will substitute the variables name 
#with their values in Makefile.in 
#so that you can use $(DRAW_CFLAGS) in Makefile.am to represent the values of the shell variable here.
AC_SUBST(DRAW_CFLAGS)
AC_SUBST(DRAW_LIBS)
AC_SUBST(MON_CFLAGS)
AC_SUBST(MON_LIBS)
#Check operating system.
AM_CONDITIONAL(LINUX,  false)
AM_CONDITIONAL(BSD,    false)
AM_CONDITIONAL(DARWIN, false)
#top_srcdir and abs_top_srcdir is empty. no idea why. but srcdir sometimes is the full path.
SRC_DIR=$(cd $srcdir && pwd)
SVN_REV=$(svnversion -n $SRC_DIR)
BUILD_DIR=$(pwd)
AC_DEFINE_UNQUOTED(SVN_REV, ["$SVN_REV"], ["Define SVN revision number"])
AC_DEFINE_UNQUOTED(SRCDIR,  ["$SRC_DIR"], ["Define SRC folder"])
AC_DEFINE_UNQUOTED(BUILDDIR,["$BUILD_DIR"],["Define BUILD folder"])
AC_DEFINE_UNQUOTED(BASEURL, ["$BASEURL"], ["The location to download files"])
OPTScommon=" -Wshadow  -Wall -fPIC" #-all-static will produce a static exe.
LDSYS=""
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
    AC_MSG_ERROR([Unrecognized host_cpu]);
esac
mkl_lib="${SRC_DIR}/external/mkl${cpu}/libmkl.la"
case "$host_os" in
 *linux*)
    AM_CONDITIONAL(LINUX, true)
    LDSYS="$LDSYS -lrt" #Real time library
    #Enable System V Interface definition (IPC,SHM)
    OPTSextra="-D_XOPEN_SOURCE=600 -D_BSD_SOURCE"
    if test -f "$mkl_lib" ;then
       use_mkl_default="yes"
       use_static_default="yes"
    else
       use_mkl_default="no"
       use_static_default="no"
    fi
    if test -d /dev/shm ;then
       use_shm_default="yes"
    else
       use_shm_default="no"
    fi
#We only have MKL redistributable for linux
AC_ARG_ENABLE(mkl,   [  --enable-mkl   Enable MKL blas/lapack], [use_mkl="$enableval"], [use_mkl="$use_mkl_default"])
AC_ARG_ENABLE(shm,   [  --disable-shm  Disable posix shared memory usage], [use_shm="$enableval"],[use_shm="$use_shm_default"])
AC_ARG_ENABLE(all-static,[  --diable-all-static Compile statically], [use_static="$enableval"],[use_static=$use_static_default])
;;
*darwin*)
 AM_CONDITIONAL(DARWIN, true)
;;
*bsd*)
 AM_CONDITIONAL(BSD, true)
;;
esac

#Common configure options.
AC_ARG_ENABLE(debug, [  --enable-debug  Turn on debug], [debug="$enableval"], [debug="no"])
AC_ARG_ENABLE(icc,   [  --enable-icc    Using intel compiler], [use_icc="$enableval"], [use_icc=$use_icc_default])
AC_ARG_ENABLE(oprof, [  --enable-oprof  Enable oprof for code profiling.],[use_oprof="$enableval"],[use_oprf="no"])
AC_ARG_ENABLE(opencc,[  --enable-opencc Enable opencc compiling], [use_opencc="$enableval"], [use_opencc="no"])

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
   AC_DEFINE(USE_POSIX_SHM, [0], "Enable posix shared memory")
fi 

AC_CHECK_LIB([pthread],[__pthread_unwind_next],[pthread_good="yes"],[pthread_good="no"])
OPTSprof="-pthread"
AC_DEFINE(_THREAD_SAFE,[1], ["Enable thread safe"])
if test "x$pthread_good" = "xyes" ;then
    AM_CONDITIONAL(DO_THRPOOL,true)
    AC_DEFINE(USE_PTHREAD,[2], ["Enable thread pool"])
else
    AM_CONDITIONAL(DO_THRPOOL,false)
    AC_DEFINE(USE_PTHREAD,[1], ["Disable thread pool"])
fi

if test x$use_oprof = "xyes" ;then
   OPTSprof="$OPTSprof -g -fno-omit-frame-pointer -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls"
fi
#Set optimize flags for different compilers
#-xO Enables SSE3, SSE2 and SSE instruction sets optimizations for non-Intel CPUs 
OPTIM_ICC="-Werror -g -O3 -xO -funroll-loops -DHAVE_SSE2=1 -wd981 -vec-report=0 " #don't set -fast globally -static-intel -static-libgcc
OPTIM_GCC="-Werror -g -O3 -funroll-loops -DHAVE_SSE2=1 -msse2 -msse -fprefetch-loop-arrays -fno-math-errno -ffinite-math-only" #-all-static panics gcc -static-libgfortran -static-libgcc
OPTIM_OPENCC="-O3"
OPT_DEBUG="-g -rdynamic -ggdb -Werror -O0 -DDEBUG_TRACE"
AM_CONDITIONAL(USE_ICC,false)
AM_CONDITIONAL(USE_OPENCC,false)
use_ifort="no"
if test "x$debug" = "xyes" ;then
    AM_CONDITIONAL(DEBUG,true)
    OPTScc="$OPT_DEBUG"
    AC_DEFINE(DEBUG,[1],["Debug the code"])
    AC_DEFINE(USE_MEM,[1],["Whether use memory management or not"])
    AC_PROG_CC
else
    AM_CONDITIONAL(DEBUG,false)
    AC_DEFINE(DEBUG,[0],["Do not debug"])
    AC_DEFINE(USE_MEM,[0],["Whether use memory management or not"])

    if test "x$use_icc" = "xyes" ;then
        OPTScc=" $OPTIM_ICC"
	CC=icc
	AC_PROG_CC([icc])
	which icc || AC_MSG_ERROR([icc is not found in path.])
	which ifort && use_ifort="yes"
	AM_CONDITIONAL(USE_ICC,true)	
	if test "x$use_opencc" = "xyes" ;then
            AC_MSG_ERROR([Shall not specify both icc and opencc])
	fi
    elif test "x$use_opencc" = "xyes" ;then
    	OPTScc="$OPTIM_OPENCC"
	CC=opencc
	AC_PROG_CC([opencc])
	AM_CONDITIONAL(USE_OPENCC,true)
    else
	OPTScc="$OPTIM_GCC"
	AC_PROG_CC
    fi	
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
AC_F77_FUNC(dgemm)
AC_F77_FUNC(dposv)
AC_F77_FUNC(dsaupd)
if test "x$use_mkl" = "xyes" ;then #Use shipped libmkl64.so. only for linux.
   if test ! -f "$mkl_lib"  ;then
      echo The MKL library is not available at "$mkl_lib"
      echo Will download one form the website. 
      echo This library is only for non-commerical use only.
      fn=mkl${cpu}_redistributed.tar.bz2
      wget ${BASEURL}/${fn}  -O /tmp/${fn}  || exit 1
      if test ! -d ${SRC_DIR}/external ;then
      	  mkdir -p ${SRC_DIR}/external || exit 1
      fi
      tar axvf /tmp/${fn} -C ${SRC_DIR}/external  || exit 1
      rm -rf /tmp/${fn}
   fi
   LAPACK="$mkl_lib ${SRC_DIR}/external/mkl${cpu}/libguide.la"
else
   if test -n "$(pkg-config --libs blas lapack)" ;then
       LAPACK=$(pkg-config --libs blas lapack)
   else
       LAPACK="-lblas -llapack"
   fi
   echo LAPACK=$LAPACK
   AC_CHECK_LIB([blas],[$dgemm],[],[AC_MSG_ERROR(Please install blas or use --enable-mkl to down MKL from MAOS website)],[$LAPACK])
   AC_CHECK_LIB([lapack],[$dposv],[],[AC_MSG_ERROR(Please install lapackor use --enable-mkl to down MKL from MAOS websit)],[$LAPACK])
fi

#Check for FFTW3 header and library.
AC_CHECK_HEADERS([fftw3.h],[],AC_MSG_ERROR([Please install fftw3 library]))
AC_CHECK_LIB([fftw3], [fftw_execute],[],[AC_MSG_ERROR(Please install fftw3)])
LDFFT="-lfftw3"
AC_CHECK_LIB([arpack], [$dsaupd],[have_arpack="yes"],[have_arpack="no"])
if test "x$have_arpack" = "xyes" ;then
   echo "Arpack is installed in the system. will use it."
   LDARPACK="-larpack"
   AM_CONDITIONAL(COMP_ARPACK, false)
else
   echo "Arpack is not available, will download the source and compile. We will require a fortran compiler."
   LDARPACK="../external/arpack/libarpack.la"
   AM_CONDITIONAL(COMP_ARPACK, true)
   if test "x$use_ifort" = "xyes";then
      AC_PROG_F77([ifort ifc]) #Use intel compiler
   else
      AC_PROG_F77 #use g77 or gfortran
   fi
   AC_PROG_F77_C_O #Check the compiler
   AC_F77_LIBRARY_LDFLAGS #create FLIBS FCLIBS
fi
AC_CHECK_LIB([cholmod],[cholmod_l_factorize], [have_chol="yes"], [have_chol="no"])
if test "x$have_chol" = "xyes" ;then
    echo "cholmod is installed in the system. will use it."
    LDCHOL="-lcholmod"
    AM_CONDITIONAL(COMP_CHOL, false)
else
    echo "cholmod is not available, will download the source and compile."
    LDCHOL="../external/cholmod/libcholmod.la"
    CFLAGS="$CFLAGS -I${srcdir}/external/cholmod/Include"
    AM_CONDITIONAL(COMP_CHOL, true)
fi
AC_CHECK_HEADERS([unistd.h])
#LDARPACK and LDCHOL must be before LAPACK and FLIBS
LDOPTS="$LDARPACK $LDCHOL $LAPACK $LDFFT -lz -lm $LDSYS $FLIBS"
#put all-static should be behind AC_F77_LIBRARY_LDFLAGS and/or AC_F77_FUNC, otherwise it will panics gcc linking with fortran.
OPTS="$OPTScommon $OPTScc $OPTSprof $OPTSextra "
AC_SUBST(OPTS)
AC_SUBST(LDOPTS)
AC_SUBST(LDEXE)
AC_SUBST(LDSYS)
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

#The following must be after AC_PROG_CC and AC_PROG_F77
AC_PROG_MAKE_SET
AC_PROG_INSTALL
AC_DISABLE_SHARED #disable shared object build on LTLIBARIRIES. should be before AC_PROG_LIBTOOL
AC_ENABLE_STATIC
AC_PROG_LIBTOOL

AC_CONFIG_FILES( [Makefile
		 lib/sys/Makefile 
		 lib/Makefile
		 maos/Makefile 
		 skyc/Makefile
		 ])
if test -d ${srcdir}/test ;then
AC_CONFIG_FILES([test/Makefile])
fi
if test -d ${srcdir}/external/arpack ;then
AC_CONFIG_FILES([external/arpack/Makefile])
fi
if test -d ${srcdir}/external/cholmod ;then
AC_CONFIG_FILES([external/cholmod/Makefile])
fi

echo use_icc=$use_icc
echo use_ifort=$use_ifort
echo use_static=$use_static
echo host_os=$host_os

echo CC=$CC
echo F77=$F77
echo FLIBS=$FLIBS
echo LIBS=$LIBS
echo LDEXE=$LDEXE
echo SRC_DIR=$SRC_DIR
echo SVN_REV=$SVN_REV
echo BUILD_DIR=$BUILD_DIR
AC_OUTPUT

#The folloing flags need to be tested before turned on.
#OPTS="$OPTS -funroll-loops"
#don't use the following.
#OPTS="$OPTS -funsafe-loop-optimizations -Wunsafe-loop-optimizations"

#-ffast-math changes the result, thus the atmosphere in gcc-4.3 (Tablet), but in polaris, the atm is the same as without. This flag is dangers. don't set it. rather, use individual flags. When something goes differently after changing compiler, check this part.
#OPTS="$OPTS -ffast-math"
#The following flags are dangerous. don't enable.
#OPTS="$OPTS -fassociative-math -fno-signed-zeros -fno-trapping-math -freciprocal-math "

#Make OPTS available in Makefile.am/in

#use FLIBS to link to fortran compiler library
#LAPACK="-lmkl_blas_lapack"
#LAPACK="-lblas -llapack"


#Notes for icc:
#-fast uses static binding, which is bad for drawdaemon
#-no-prec-div causes error
#-ipo causes dnew etc not found
#-axT 
#-xHOST causes illegal operations
#

#2009-12-28:
#Mac related: autoreconf -vfi in mac complains about the missing of AM_CONDITIONAL
#Copyed m4 folder in polaris to the mac. succeed. need to distrubite m4 directory.

#The grammar in configure.ac and configure is similar but not identical to shell
#scripts.  in if, don't use [, must use test because [ is treated as quote my
#autoconf.

#Maintainer mode causes configure to not regenerate Makefile when makefile.am is changed
#AM_MAINTAINER_MODE

#The bug I had in test_accphi is due to the fact that I didn't set  -DHAVE_SSE2=1 for icc.
#remove presented cflags
#We should not modify user variables CFLAGS LDFLAGS etc.
   #LAPACK="-Wl,-rpath,${SRC_DIR}/external/mkl -L${SRC_DIR}/external/mkl -lmkl64"

AC_DEFUN([MY_COMPILE],[
	echo Download and compile $1
	stop=0
	name=$1
	file=$2
	rm -rf ${TMP_DIR}/$name
	MY_DOWNLOAD([$name], [$file], [${TMP_DIR}/$name], [1])
	fnlog=${TMP_DIR}/compile.log
	fnlog2=${TMP_DIR}/compile_${name}.log
	(
		echo Compiling $name with $3
		if cd ${TMP_DIR}/$name/* ;then
			export CFLAGS="${CFLAGS/-std=c11/} -O3 -Wno-error" #-std=c11 causes error in compiling libwebsockets 4.4.1 on macos
			export CPPFLAGS="$CPPFLAGS" 
			ans=0
			if test -f configure ;then
				(./configure --prefix=${prefix} --libdir=${libdir} $3 && make -j4 install) > $fnlog2 2>&1|| ans=1
			elif test -f CMakeLists.txt ;then
				(mkdir build && cd build && cmake .. -DCMAKE_INSTALL_PREFIX:PATH=${prefix} -DCMAKE_INSTALL_LIBDIR=lib${libsuffix} -DLIB_SUFFIX=${libsuffix} $3 && make -j4 install) >$fnlog2 2>&1|| ans=1
			else
				echo "Unknown compile system." 
				ans=1
			fi
			if test $ans = 0 ;then
				echo "Compiled $name successfully."
			else
				echo "Failed to compile $name. Please check $fnlog and $fnlog2 for details." 
			fi
		fi
	) | tee -a $fnlog
])

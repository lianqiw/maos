AC_DEFUN([MY_DOWNLOAD],[
	stop=0
	fn="$2"
	OUTDIR="$3"
	if test "$enable_download" = "no" ;then
		AC_MSG_NOTICE([Download is disabled])
		stop=1
	fi
	if test -z "$OUTDIR" ;then
		AC_MSG_NOTICE([Download folder is not specified])
		stop=1
	else
		mkdir -p "${OUTDIR}" || stop=1
	fi
	mkdir -p "${DOWNLOAD_DIR}" || stop=1
	if test ${stop} = 0 -a "$enable_binary" = "no" -a -z "$4" ;then
		AC_MSG_NOTICE([Download pre-built binary is disabled])
		stop=1
	fi
	if test ${stop} = 0 ;then
		if ! test -f "${DOWNLOAD_DIR}/${fn}" ;then
			case "$fn" in
			http*)
				url="$fn"
				fn=`basename $fn`
				;;
			*)
				url=${BASEURL}/${fn}
				;;
			esac
			echo $wget "$url" "$wgetoutput" "${DOWNLOAD_DIR}/${fn}"	
			$wget "$url" "$wgetoutput" "${DOWNLOAD_DIR}/${fn}" || stop=1
		fi
		if test ${stop} = 0 -a -f "${DOWNLOAD_DIR}/${fn}" ;then
			if ! tar -xf "${DOWNLOAD_DIR}/${fn}" -C "${OUTDIR}" ;then
				AC_MSG_NOTICE([Failed to extract the file from ${BASEURL}/${fn}])
				rm -rf ${DOWNLOAD_DIR}/${fn}
				stop=1
			fi
		fi
	fi
])

#!/bin/sh
#Clean maos results folder, leaving only Res_${seed}.bin, Res_${seed}.done, *.conf, *.log

for fd in "$@" ;do
    (
	cd ${fd}
	if ls Res_*.bin >/dev/null 2>&1 ;then
	    echo cleaning ${fd}
	    rm -rfv setup Res[ofzC]*_*.bin Timing_*.bin psd*.bin extra
	else
	    for sub in * ;do
		if [ -d $sub ];then
		    eval $0 $sub
		fi
	    done
	fi
    )
done
    

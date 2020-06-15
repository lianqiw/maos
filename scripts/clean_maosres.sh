#!/bin/sh
#Clean maos results folder, leaving only Res_${seed}.bin, Res_${seed}.done, *.conf, *.log

for fd in "$@" ;do
    (
	cd ${fd}
	if ls Res_*.bin >/dev/null;then
	    echo cleaning ${fd}
	    rm -rf setup Res[a-zA-Z]*.bin Timing_*.bin psd*.bin
	else
	    for sub in * ;do
		if [ -d $sub ];then
		    eval $0 $sub
		fi
	    done
	fi
    )
done
    

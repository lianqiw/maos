#!/bin/sh
#Clean maos results folder, leaving only Res_${seed}.bin, Res_${seed}.done, *.conf, *.log
pwd=`pwd`
function fixlink(){
	#echo "$1 $2 $3"
	if [ -L "$1" ];then
		mv -f $(readlink "$1") "$1"
	fi
	if [ ! -f "$1" ];then
		fn=`ls -t "$2" 2>/dev/null |head -n1` 
		if [ -f "$fn" ];then
			mv -f $fn "$1"
		fi
	fi
	if [ -f "$1" ] ;then
		rm -rf "$2" "$3"
	fi
}

for fd in "$@" ;do
	cd $pwd/${fd}
	if ls Res_[0-9]*.bin >/dev/null 2>&1 ;then #maos results
	    echo cleaning maos $pwd/${fd}
		if ls extra/Resp_[0-9]*.bin >/dev/null 2>&1;then
			mv -n extra/Resp_[0-9]*.bin . 
		fi
		#keep Res_*.bin and Resp_*.bin
		if [ -f setup/aloc ];then
			rm -rf setup
		fi
		if ls extra/Resclemp_*.bin > /dev/null 2>&1 ; then
			rm -rf extra
		fi
	    rm -rf Res[ofzcCu]*_[0-9]*.bin Res[u]*_[0-9]*.txt Timing_*.bin  powfs[0-9]*_*.bin otf*.bin range_*.bin wfs[0-9]_ints*.bin wfs[0-9]_grad*.bin dm_ncpa.bin surfevl_*.bin
		for fn in run_*.log ;do
			if grep 'There are no seed to run.' $fn >/dev/null 2>&1 ;then
				rm -rf $fn
			fi
		done
		if ls Res_[0-9].done ;then
			fixlink maos_done.conf 'maos*_[0-9]*.conf'
			fixlink run_done.log 'run_[0-9]*.done' 'run_[0-9]*.log'
		fi
		rm -rf run*_[0-9]*.killed
	elif ls Res[0-9]*_[0-9]*.bin >/dev/null 2>&1 ;then #skyc
 		echo cleaning skyc $pwd/${fd}
		rm -rf setup
		if ls Res[0-9]*_[0-9]*.conf ;then
			fixlink skyc_done.conf 'skyc*_[0-9]*.conf'
			fixlink run_done.log 'run_[0-9]*.done' 'run_[0-9]*.log'
		fi
		rm -rf run*_[0-9]*.killed
	fi

	for sub in * ;do
		if [ -d "$sub" ];then
			eval $0 "$sub"
		fi
	done
done
    

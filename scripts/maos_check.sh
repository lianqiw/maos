#!/bin/bash
if [ ! -f "../bin/maos" ];then
	echo "Please run in maos build directory"
	exit
fi
unset MAOS_PYWFS_DEBUG 
export MAOS_LOG_LEVEL=-1 #reduce log output
D=30
#Sanity check maos
#NFIRAOS default
shopt -s extglob
case "$1" in
	(+([0-9.]))
	D=$1
	shift;;
esac
if [ "$D" != "30" ];then
    args="aper.d=$D "
fi
args+=" $@"

ARCH=CPU
case "$@" in
  *-g-1* | *-G0*)
	;;
  *)
	if nvidia-smi -i 0 >/dev/null 2>&1 ;then
		ARCH=GPU
	fi
		;;
esac

case $D in
    2)
		REF_GPU=(595.58 94.55 103.43 91.58 91.27 129.20 118.57 136.70 120.19 100.21 98.93 322.83 322.79 98.98 94.81 119.55 113.42 114.60 112.48 121.30 121.39 126.07 124.68 123.30 121.81 117.79 130.26 0.00) #2025-01-13 kepler D=2
		REF_CPU=(595.58 90.39 105.95 90.46 90.65 129.43 117.37 136.45 119.85 96.17 98.42 329.94 324.01 99.60 95.10 120.98 113.52 143.84 115.39 123.43 123.78 123.58 126.62 124.33 123.21 118.97 133.20 0.00) #2025-01-13 kepler D=2
	;;
    5)
		REF_GPU=(1092.86 111.69 107.94 100.62 100.45 115.40 115.41 109.80 109.23 115.92 130.11 538.77 549.30 113.25 108.71 143.20 122.27 121.47 126.00 153.06 152.35 153.63 152.14 154.36 153.32 130.37 168.12 0.00) #2025-01-13 kepler D=5
		REF_CPU=(1092.86 106.82 108.46 100.59 100.02 115.85 114.75 110.59 109.27 116.00 128.95 537.98 549.50 111.81 109.02 143.25 121.90 154.93 129.90 154.03 153.27 153.67 153.93 159.34 154.52 128.53 177.02 0.00) #2025-01-13 kepler D=5
	;;
    10)
		REF_GPU=(1338.49 110.78 102.43 102.81 103.97 111.84 111.04 110.75 108.60 166.95 182.54 837.98 838.58 113.11 111.80 111.05 120.90 121.29 124.32 128.59 127.79 129.53 129.71 127.37 128.60 127.97 151.12 0.01) #2025-01-13 kepler D=10
		REF_CPU=(1338.49 109.63 103.12 103.03 102.48 111.57 111.06 110.78 108.53 167.39 181.36 838.03 838.54 108.09 110.60 111.75 120.86 140.25 123.98 128.17 128.47 129.25 129.20 134.70 128.16 127.91 146.99 0.01) #2025-01-13 kepler D=10
	;;
    30)
		REF_GPU=(1700.11 113.54 121.63 112.56 113.78 118.27 116.09 129.27 117.33 345.00 365.09 867.34 898.31 119.65 124.73 122.60 127.68 132.64 138.15 141.25 138.22 143.45 140.56 141.35 138.78 275.42 0.86 188.81 148.71) #2025-01-13 kepler D=30
		REF_CPU=(1700.11 113.43 120.17 112.25 112.39 118.33 116.11 129.20 117.47 345.40 363.34 866.73 898.58 118.06 125.94 122.29 127.57 131.92 138.06 141.40 138.41 141.63 0 141.26 137.37 269.45 0.92 185.90 0) #2025-01-13 kepler D=30
	;;
    *)
		REF_GPU=()
		REF_CPU=()
	;;
esac
echo "Using $ARCH"
if [ "$ARCH" = GPU ];then
	REF=(${REF_GPU[@]})
	REF2=(${REF_CPU[@]})
else
	REF=(${REF_CPU[@]})
	REF2=(${REF_GPU[@]})
fi
fnlog=maos_check_${D}.log #log of all
fntmp=maos_check_${D}.tmp #log of current simulation
fnerr=maos_check_${D}.err #err of current simulation
fnres=maos_check_${D}.res #result summary
#fnerr=maos_check_${D}.err

echo $(date) > $fnlog
#echo $(date) > $fnerr
ans=0 #result code
ii=0
printf "%-20s    Res    Ref  Diff Diff2\n" "D=${D}m" | tee $fnres
function run_maos(){
	aotype=$1
	shift
    eval "../bin/maos sim.end=100 $* $args >$fntmp 2>$fnerr"
    if [ $? -eq 0 ];then
		RMS[ii]=$(grep 'Mean:' $fntmp |tail -n1 |cut -d ' ' -f 2)
		a=${RMS[$ii]}
	else
		a=
	fi
	if [ x$a = x ];then
		cat $fntmp >> $fnlog
		RMS[ii]=000.00
		a=0
		ans=$((ans+1)) #failed to run
	fi
    
	echo $aotype $* >> $fnlog
	cat $fntmp >> $fnlog
    b=${REF[$ii]:-0}
    b2=${REF2[$ii]:-0}
	diff=$(echo "200*($a-$b)/($a+$b+1)" | bc)
	diff2=$(echo "200*($a-$b2)/($a+$b2+1)" | bc)
	if [ $diff -gt 10 -o $diff -lt -10 ];then
		ans=$((ans+1)) #mark failure
	fi
	printf "%-20s %6.1f %6.1f %4d%% %4d%%\n" "$aotype" "$a" "$b" "$diff" "$diff2"| tee -a $fnres
	ii=$((ii+1)) 
}
function run_maos_gpu(){
	if [ "$ARCH" = GPU ];then
		run_maos "$@"
	else
		printf "%-20s   skipped in CPU mode.\n" "$1" | tee -a $fnres
		RMS[ii]=0
		ii=$((ii+1)) 
	fi
}

run_maos "Openloop:        " -cmcao_lgs.conf sim.evlol=1

run_maos "Ideal tomo:      " -cmcao_lgs.conf sim.idealtomo=1 

run_maos "Evl tomo:        " -cmcao_lgs.conf evl.tomo=1 

run_maos "NGS SCAO (inte): " -cscao_ngs.conf recon.split=0

run_maos "NGS SCAO (ahst):" -cscao_ngs.conf 

run_maos "NGS SCAO (lsq):  " -cscao_ngs.conf recon.alg=1

run_maos "NGS SCAO (lsq,modal)" -cscao_ngs.conf recon.modal=1 recon.alg=1

run_maos "NGS PYWFS (zonal):" -cscao_pywfs.conf recon.modal=0 sim.end=500

run_maos "NGS PYWFS (modal):" -cscao_pywfs.conf recon.modal=1 sim.end=500

run_maos "LGS SCAO (inte): " -cscao_lgs.conf recon.split=0

run_maos "LGS SCAO (ahst):" -cscao_lgs.conf 

run_maos "LGS GLAO (inte): " -cglao.conf recon.split=0 evl.psfmean=0

run_maos "LGS GLAO (ahst):" -cglao.conf  evl.psfmean=0

run_maos "LGS LTAO (inte): " -cmcao_lgs.conf dm_single.conf fov_oa.conf recon.split=0

run_maos "LGS LTAO (ahst):" -cmcao_lgs.conf dm_single.conf fov_oa.conf powfs.astscale=[1 0.1 0.1] #ahst with wide NGS asterism suffers (from aliasing error?)

run_maos "LGS MOAO:        " -cmcao_lgs.conf evl.moao=0 moao.dx=[1/2]

run_maos "NGS MCAO (inte): " -cmcao_ngs.conf recon.split=0

run_maos "NGS MCAO (ahst):" -cmcao_ngs.conf 

run_maos "LGS MCAO (inte): " -cmcao_lgs.conf recon.split=0 tomo.precond=0

run_maos "LGS MCAO (CG):   " -cmcao_lgs.conf tomo.precond=0 cn2.pair=[0 1 2 5] recon.psd=1 tomo.assemble=0 fit.assemble=1

run_maos "LGS MCAO (FDPCG):" -cmcao_lgs.conf tomo.precond=1 tomo.assemble=1 fit.assemble=0

run_maos "LGS MCAO (CBS):  " -cmcao_lgs.conf tomo.alg=0 fit.alg=0 atmr.os=[2 2 1 1 1 1 1]

if [ ${D%.*} -le 10 ];then
run_maos "LGS MCAO (SVD):  " -cmcao_lgs.conf tomo.alg=2 fit.alg=2 atmr.os=[2 2 1 1 1 1 1]
run_maos "LGS MCAO (MVM):  " -cmcao_lgs.conf atmr.os=[2] tomo.precond=1 tomo.maxit=100 fit.alg=0 recon.mvm=1
else
run_maos_gpu "LGS MCAO (MVM):  " -cmcao_lgs.conf atmr.os=[2] tomo.precond=1 tomo.maxit=100 fit.alg=0 recon.mvm=1
fi
run_maos "LGS MCAO PCCD:  " -cmcao_lgs.conf tomo.precond=0 cn2.pair=[0 1 2 5] recon.psd=1 powfs.radpix=[16,0,0] powfs.pixpsa=[6,0,0]

run_maos "SLGS MCAO (inte): " -cmcao_lgs.conf cn2.pair=[0 1 2 5] recon.psd=1 powfs.fnllt=['llt_SL.conf',,] powfs.pixpsa=[16,0,0] recon.split=0

run_maos "SLGS MCAO (ahst):  " -cmcao_lgs.conf cn2.pair=[0 1 2 5] recon.psd=1 powfs.fnllt=['llt_SL.conf',,] powfs.pixpsa=[16,0,0]

run_maos "LGS NCPA noatm: " -cmcao_lgs.conf sim.noatm=1 ncpa.surf=["'rms=150;L0=20;D=60;SURFEVL=0;'"] sim.wspsd= powfs.noisy=[0] atm.r0z=10 powfs0_llt.fnsurf="'rms=150;mode=5;D=0.5;dx=1/64'" #should be close to 0.

if [ "$D" = 30 ];then
run_maos "NFIRAOS LGS:    " -cnfiraos_lgs.conf
run_maos_gpu "NFIRAOS PYWFS:" -cnfiraos_ngs.conf
fi
echo "REF_${ARCH}=(${RMS[*]}) #$(date +%Y-%m-%d) $HOSTNAME D=$D"
if [ $ans -ne 0 ];then
	echo "$ans tests failed"
fi
exit $ans

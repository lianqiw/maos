#!/bin/bash
if [ ! -f "./maos" ];then
	echo "Please run in maos bin directory"
	exit
fi
#Sanity check maos
#NFIRAOS default
if [ -n "$1" ];then
    D=$1
    shift
    args="aper.d=$D $@"
else
    D=30
fi
has_gpu=0
if nvidia-smi -i 0 >/dev/null 2>&1 ;then
	case "$@" in
		"-g-1" | "-G0")
		has_gpu=0
		;;
	*)
		has_gpu=1
		;;
	esac
fi
case $D in
    2)
		REF_GPU=(624.74 87.53 86.32 116.18 124.88 126.96 137.45 133.29 154.33 125.97 333.07 333.78 96.20 96.05 112.54 114.53 111.13 108.47 182.17 215.01 102.20 116.27 124.01 134.59) #2022-03-08 Cassiopeia
		REF_CPU=(624.74 85.64 86.33 111.53 125.08 125.07 130.34 132.90      0 123.06 331.32 333.89 96.26 97.20 113.42 142.91 104.08 107.89 143.15 215.34 101.22 113.99 125.71 135.00) #2022-03-08 Cephei
	;;
    5)
		REF_GPU=(1137.59 107.22 106.03 125.50 151.39 153.36 155.71 153.64 153.89 142.64 573.55 567.56 101.01 104.74 121.33 125.36 117.62 129.99 132.54 203.98 104.87 110.65 151.66 164.96) #2022-03-08 Cassiopeia
		REF_CPU=(1137.59 105.71 106.08 126.67 151.92 152.54 155.20 153.08      0 143.19 573.86 567.80 102.00 103.79 121.59 161.44 117.04 129.37 118.53 202.99 104.99 110.05 152.71 170.46) #2022-03-08 Cephei
	;;
    10)
		REF_GPU=(1353.41 109.88 109.72 123.89 126.71 126.02 128.62 128.16 127.90 110.34 720.58 714.38 102.74 103.28 120.79 127.73 167.76 180.87 125.08 191.02 106.56 107.49 126.50 145.43) #2022-03-08 Cassiopeia
		REF_CPU=(1353.41 109.42 109.74 124.20 126.05 126.09 128.72 128.05      0 109.93 720.61 714.33 102.88 103.21 120.50 140.38 166.79 182.87 112.89 190.40 106.03 107.67 126.25 143.75) #2022-03-08 Cephei
	;;
    30)
		REF_GPU=(1712.96 113.71 113.81 137.65 144.37 137.99 143.05 140.07 122.93 870.52 883.52 112.75 113.77 128.27 141.01 346.93 366.19 122.84 140.14 117.64 118.05 144.32 197.92 199.76 91.97) #2022-09-26 Maxwell
		REF_CPU=(1712.96 113.42 113.71 137.54 144.39 137.80 140.94      0 122.64 870.77 883.79 112.41 113.42 127.60 130.98 346.06 365.48 124.62 140.81 116.24 118.47 144.18 193.61 196.27 0) #2022-03-08 Cephei
	;;
    *)
		REF_GPU=()
		REF_CPU=()
	;;
esac
if [ $has_gpu -eq 1 ];then
	echo "Using GPU" 
	REF=(${REF_GPU[@]})
else
	echo "Using CPU"
	REF=(${REF_CPU[@]})
fi
fnlog=maos_check_${D}.log #log of all
fntmp=maos_check_${D}.tmp #log of current simulation
fnres=maos_check_${D}.res #result summary
#fnerr=maos_check_${D}.err

echo $(date) > $fnlog
#echo $(date) > $fnerr

ii=0
printf "%-20s   Res   Ref     %%\n" "D=${D}m DM is $((D*2))x$((D*2))" | tee $fnres
function run_maos(){
	aotype=$1
	shift
    cmd="./maos sim.end=100 $* $args > $fntmp 2>&1"
	eval $cmd
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
		ans=1
	fi
    
	echo $aotype $* >> $fnlog
	cat $fntmp >> $fnlog
    b=${REF[$ii]}
	if [ x$b = 'xerror' -o x$b = x ];then
		b=0
	fi
	as=${a/./} #%.*
	bs=${b/./}
    if [ "$as" -gt 0 -a "$bs" -gt 0 ];then
		diff=$(((as-bs)*100/as))
	else
		diff=100
    fi
	printf "%-20s %6.1f %6.1f %4d%%\n" "$aotype" "$a" "$b" "$diff" | tee -a $fnres
	ii=$((ii+1)) 
}
function run_maos_gpu(){
	if [ $has_gpu -eq 1 ];then
		run_maos "$@"
	else
		printf "%-20s   skipped in CPU mode.\n" "$1" | tee -a $fnres
		RMS[ii]=0
		ii=$((ii+1)) 
	fi
}
export MAOS_LOG_LEVEL=-1

run_maos "Openloop:        " -cmcao_lgs.conf sim.evlol=1

run_maos "Ideal fit:       " -cmcao_lgs.conf sim.idealfit=1 

run_maos "Ideal tomo:      " -cmcao_lgs.conf sim.idealtomo=1 

run_maos "LGS MCAO (inte): " -cmcao_lgs.conf recon.split=0 tomo.precond=0

run_maos "LGS MCAO (CG):   " -cmcao_lgs.conf tomo.precond=0 cn2.pair=[0 1 2 5] recon.psd=1 tomo.assemble=0 fit.assemble=1

run_maos "LGS MCAO (FDPCG):" -cmcao_lgs.conf tomo.precond=1 tomo.assemble=1 fit.assemble=0

run_maos "LGS MCAO (CBS):  " -cmcao_lgs.conf tomo.alg=0 fit.alg=0 atmr.os=[2 2 1 1 1 1 1]

if [ $D -le 10 ];then
run_maos "LGS MCAO (SVD):  " -cmcao_lgs.conf tomo.alg=2 fit.alg=2 atmr.os=[2 2 1 1 1 1 1] gpu.tomo=0
fi
run_maos_gpu "LGS MCAO (MVM):  " -cmcao_lgs.conf atmr.os=[2] tomo.precond=1 tomo.maxit=100 fit.alg=0 recon.mvm=1

run_maos "LGS MOAO:        " -cmcao_lgs.conf evl.moao=0 moao.dx=[1/2]

run_maos "LGS GLAO (inte): " -cmcao_lgs.conf glao.conf recon.split=0 evl.psfmean=0

run_maos "LGS GLAO (split):" -cmcao_lgs.conf glao.conf recon.split=1 evl.psfmean=0

run_maos "NGS SCAO (inte): " -cscao_ngs.conf recon.split=0

run_maos "NGS SCAO (split):" -cscao_ngs.conf recon.split=1

run_maos "NGS MCAO (inte): " -cmcao_ngs.conf recon.split=0

run_maos "NGS MCAO (split):" -cmcao_ngs.conf recon.split=1

run_maos "SCAO LGS (inte): " -cscao_lgs.conf recon.split=0

run_maos "SCAO LGS (split):" -cscao_lgs.conf recon.split=1

run_maos "LGS LTAO (inte): " -cmcao_lgs.conf dm_single.conf fov_oa.conf recon.split=0

run_maos "LGS LTAO (split):" -cmcao_lgs.conf dm_single.conf fov_oa.conf recon.split=1 

run_maos "NGS SCAO (lsq,inte)" -cscao_ngs.conf recon.split=0 recon.alg=1

run_maos "NGS SCAO (lsq,split)" -cscao_ngs.conf recon.split=1 recon.alg=1

run_maos "LGS MCAO PCCD:  " -cmcao_lgs.conf tomo.precond=0 cn2.pair=[0 1 2 5] recon.psd=1 powfs.radpix=[16,0,0] powfs.pixpsa=[6,0,0]

run_maos "LGS MCAO SL:    " -cmcao_lgs.conf tomo.precond=0 cn2.pair=[0 1 2 5] recon.psd=1 powfs.fnllt=['llt_SL.conf',,] powfs.pixpsa=[16,0,0]

if [ $D -eq 30 ];then
run_maos "NFIRAOS LGS: "	  -c nfiraos_lgs_full.conf
run_maos_gpu "NFIRAOS PYWFS:" -c nfiraos_ngs.conf
fi
echo "REF=(${RMS[*]})"

exit $ans

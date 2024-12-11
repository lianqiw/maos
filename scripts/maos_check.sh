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

has_gpu=0

case "$@" in
  *-g-1* | *-G0*)
	has_gpu=0
	ARCH=CPU
	;;
  *)
	if nvidia-smi -i 0 >/dev/null 2>&1 ;then
		has_gpu=1
		ARCH=GPU
	fi
		;;
esac

case $D in
    2)
		#REF_GPU=(624.74 94.57 182.51 95.92 95.14 132.22 120.51 150.92 119.60 110.50 106.36 356.67 345.91 186.97 215.25 124.40 114.22 113.37 116.76 124.62 127.07 129.90 127.24 128.57 124.48 133.21 0 0.00) #2024-06-04 maxwell 2
		REF_GPU=(595.58 87.44 104.23 92.27 93.51 127.71 117.29 142.71 123.74 99.62 126.08 313.00 309.61 102.40 143.76 144.74 111.72 119.35 110.69 152.05 147.65 162.55 157.73 164.60 158.69 115.99 175.62 0.02) #2024-12-13 kepler D=2
		REF_CPU=(624.74 90.45 160.47 95.97 96.98 133.59 121.36 151.82 118.62 102.58 107.25 362.80 351.38 151.03 215.30 123.58 113.82 146.38 111.57 125.24 125.46 128.21 127.38 129.47 124.94 136.33 0 0.00) #2024-06-04 kepler D=2
	;;
    5)
		#REF_GPU=(1137.59 111.83 131.19 101.13 102.72 117.49 116.35 114.50 113.45 116.48 129.23 585.33 579.60 130.44 203.05 142.52 121.86 125.00 125.55 153.03 152.77 153.19 153.82 154.91 152.86 164.90 0 0.00) #2024-06-04 maxwell 5
		REF_GPU=(1092.86 105.96 108.17 101.55 103.67 115.40 114.76 111.68 110.19 117.65 131.98 544.46 554.66 112.72 112.47 145.50 121.88 125.83 120.96 154.78 153.71 155.25 157.41 155.50 155.57 133.00 175.48 0.00) #2024-12-17 kepler D=5
		REF_CPU=(1137.59 106.92 119.81 100.64 103.08 116.96 115.86 114.51 113.84 116.69 129.66 583.40 579.29 117.97 203.56 144.09 121.69 161.44 126.82 152.88 153.17 153.51 155.28 156.23 153.78 169.79 0 0.00) #2024-06-04 maxwell
	;;
    10)
		#REF_GPU=(1353.41 110.89 125.10 102.36 103.50 111.80 111.01 112.45 107.69 166.37 181.06 734.84 732.12 125.45 191.27 111.03 121.09 126.70 124.02 127.63 126.99 128.98 129.18 127.79 127.64 144.37 0 0.01) #2024-06-04 maxwell
		REF_GPU=(1338.49 109.47 102.26 103.43 103.63 111.96 111.11 112.91 110.32 166.95 182.43 599.23 592.16 108.24 111.49 110.84 121.01 121.93 124.18 127.31 126.92 129.57 128.56 128.31 128.20 127.96 154.55 0.01) #2024-12-13 maxwell D=10
		REF_CPU=(1353.41 109.74 110.27 102.53 102.67 111.00 111.46 112.43 107.62 166.71 180.60 735.92 731.93 112.67 189.96 111.22 120.62 140.85 124.66 127.12 127.08 128.61 128.51 145.28 127.51 142.65 0 0.01) #2024-06-04 kepler
	;;
    30)
		#REF_GPU=(1711.79 113.72 128.88 112.47 112.73 119.55 117.22 380.73 451.25 345.98 364.96 872.14 883.45 124.15 141.63 122.42 128.13 140.03 137.25 144.70 137.57 142.94 139.96 144.23 197.38 0 0.86 190.49 159.72) #2024-06-04 maxwell
		REF_GPU=(1700.11 113.31 122.14 112.68 113.50 118.27 117.64 129.41 117.14 344.59 365.15 810.84 837.82 118.12 125.26 124.07 127.67 133.97 137.98 141.85 139.29 143.84 141.32 141.54 138.26 287.03 0.86 192.26 150.94) #2024-12-13 kepler D=30
		REF_CPU=(1711.79 113.62 126.09 113.09 112.95 119.37 117.29 381.48 451.14 343.83 366.23 871.85 883.12 120.87 141.41 122.23 127.49 132.50 137.47 144.18 138.34 141.65      0 144.65 195.68 0 0.92 189.62 0) #2024-06-04 kepler
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
fnerr=maos_check_${D}.err #err of current simulation
fnres=maos_check_${D}.res #result summary
#fnerr=maos_check_${D}.err

echo $(date) > $fnlog
#echo $(date) > $fnerr
ans=0 #result code
ii=0
printf "%-20s    Res    Ref  Diff\n" "D=${D}m" | tee $fnres
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
    b=${REF[$ii]}
	if [ x$b = 'xerror' -o x$b = x ];then
		b=0
	fi
	as=${a#0.}
	bs=${b#0.}
	as=${as%.*} 
	bs=${bs%.*}
    if [ "$as" -eq "$bs" ];then #in case both are zero
		diff=0
	else
		diff=$(((as-bs)*200/(as+bs)))
    fi
	if [ $diff -gt 5 -o $diff -lt -5 ];then
		ans=$((ans+1)) #mark failure
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

run_maos "LGS GLAO (inte): " -cmcao_lgs.conf glao.conf recon.split=0 evl.psfmean=0

run_maos "LGS GLAO (ahst):" -cmcao_lgs.conf glao.conf  evl.psfmean=0

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

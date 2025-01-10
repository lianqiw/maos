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
		REF_GPU=(595.58 87.44 104.23 92.27 93.51 127.71 117.29 142.71 123.74 99.62 126.08 313.00 309.61 102.40 143.76 144.74 111.72 119.35 110.69 152.05 147.65 162.55 157.73 164.60 158.69 115.99 175.62 0.02) #2024-12-13 kepler D=2
		REF_CPU=(595.58 86.97 105.54 92.79 92.28 128.36 117.82 142.83 122.95 97.10 128.16 322.78 341.68 99.39 140.22 150.00 112.63 144.85 114.63 160.61 155.88 172.26 160.66 182.36 148.58 117.51 190.69 0.00) #2024-12-20 Capri D=2
	;;
    5)
		REF_GPU=(1092.86 105.96 108.17 101.55 103.67 115.40 114.76 111.68 110.19 117.65 131.98 544.46 554.66 112.72 112.47 145.50 121.88 125.83 120.96 154.78 153.71 155.25 157.41 155.50 155.57 133.00 175.48 0.00) #2024-12-17 kepler D=5
		REF_CPU=(1092.86 105.59 107.85 100.98 102.30 115.85 115.19 111.83 110.24 117.01 132.14 540.29 550.97 110.52 111.38 146.76 122.14 161.00 122.61 155.61 153.78 157.31 158.33 156.79 154.99 131.08 186.09 0.00) #2024-12-20 Capri D=5
	;;
    10)
		REF_GPU=(1338.49 109.47 102.26 103.43 103.63 111.96 111.11 112.91 110.32 166.95 182.43 599.23 592.16 108.24 111.49 110.84 121.01 121.93 124.18 127.31 126.92 129.57 128.56 128.31 128.20 127.96 154.55 0.01) #2024-12-13 maxwell D=10
		REF_CPU=(1338.49 109.32 103.01 102.99 103.24 111.57 111.32 112.68 110.31 167.47 181.89 837.67 838.39 108.22 110.91 110.75 120.78 141.39 124.06 127.51 127.41 129.43 128.57 138.60 127.65 127.88 147.77 0.01) #2024-12-20 Capri D=10
	;;
    30)
		REF_GPU=(1700.11 113.31 122.14 112.68 113.50 118.27 117.64 129.41 117.14 344.59 365.15 810.84 837.82 118.12 125.26 124.07 127.67 133.97 137.98 141.85 139.29 143.84 141.32 141.54 138.26 287.03 0.86 192.26 150.94) #2024-12-13 kepler D=30
		REF_CPU=(1700.11 113.26 121.78 114.54 112.15 118.68 116.67 129.33 117.33 343.56 363.98 881.55 902.45 118.14 124.91 123.00 127.44 132.65 137.34 142.03 139.03 142.20 0 141.75 137.66 282.60 0.92 185.82 0) #2024-12-19 kepler D=30
	;;
    *)
		REF_GPU=()
		REF_CPU=()
	;;
esac
echo "Using $ARCH"
if [ "$ARCH" = GPU ];then
	REF=(${REF_GPU[@]})
else
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
	diff=$(echo "200*($a-$b)/($a+$b)" | bc)
	if [ $diff -gt 10 -o $diff -lt -10 ];then
		ans=$((ans+1)) #mark failure
	fi
	printf "%-20s %6.1f %6.1f %4d%%\n" "$aotype" "$a" "$b" "$diff" | tee -a $fnres
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

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
    5)
	REF_GPU=(1153.3 104.2 107.0 120.2 121.8 120.2 117.3 124.7 138.4 615.7 617.3 127.8 114.4 138.3 141.4 146.1 128.6 129.6 149.7 154.3 154.6 156.9 160.3 158.9 154.3 152.6 160.7 0.0) #2025-04-24 kepler D=5
	REF_CPU=(1153.3 104.3 104.1 121.4 122.1 120.8 117.3 126.2 138.4 615.1 616.5 130.2 115.4 144.2 141.7 145.7 123.3 137.7 161.0 157.4 160.6 161.4 161.3 156.5 161.7 163.0 162.9 0.0) #2025-04-24 kepler D=5
	#REF_GPU=(1153.27 105.26 104.49 120.00 122.03 120.60 117.21 125.56 139.31 613.02 616.78 126.97 116.08 164.02 141.30 146.58 128.64 127.86 149.53 173.27 175.16 173.93 174.33 175.88 175.59 155.10 194.71 0.01) #2025-01-27 kepler D=5
	#REF_CPU=(1153.27 105.08 105.42 120.30 121.94 121.11 117.10 125.73 138.63 611.52 616.12 131.02 115.92 165.82 141.53 146.26 123.33 134.69 155.19 177.18 176.81 178.04 178.07 178.92 178.23 157.58 202.46 0.01) #2025-01-27 maxwell D=5
	;;
    10)
	REF_GPU=(1431.3 109.4 110.2 119.3 118.7 119.8 116.3 178.8 193.1 934.5 942.6 120.7 117.2 125.6 138.8 139.5 128.3 115.8 144.3 150.2 149.7 151.5 153.1 155.5 150.1 149.8 156.2 0.0) #2025-04-24 kepler D=10
	REF_CPU=(1431.3 108.8 109.5 118.8 119.2 119.4 116.3 178.5 194.1 933.9 941.8 120.7 116.4 124.1 139.0 139.7 127.0 117.5 144.7 150.0 149.4 150.1 152.0 155.1 151.0 149.3 151.6 0.0) #2025-04-24 kepler D=10
	#REF_GPU=(1431.26 109.31 110.02 118.17 118.48 118.32 115.89 179.06 193.30 932.46 942.26 117.90 122.58 125.13 139.32 139.65 128.33 114.60 143.01 149.91 148.68 150.68 151.21 149.23 149.12 148.23 164.41 0.02) #2025-01-27 kepler D=10
	#REF_CPU=(1431.26 110.44 109.16 118.80 119.05 118.35 115.86 178.90 193.08 932.60 942.34 116.69 121.43 125.43 139.10 139.88 127.04 115.26 143.27 149.65 149.05 150.54 150.29 149.33 149.50 148.27 161.88 0.02) #2025-01-27 maxwell D=10
	;;
    30)	
	REF_GPU=(1819.6 119.2 120.9 125.5 123.1 144.1 127.9 367.5 388.1 814.0 823.5 128.6 125.8 131.5 145.5 151.8 132.2 129.6 155.5 156.8 156.5 160.4 162.8 156.8 156.8 158.2 0.6 200.1 159.7) #2025-04-24 kepler D=30
	REF_CPU=(1819.6 119.0 118.1 125.2 123.4 144.2 127.8 367.6 389.4 814.0 823.4 128.8 126.9 132.7 145.6 147.9 132.1 130.5 155.3 156.6 157.2 160.3 000.0 157.6 156.9 159.0 0.7 203.4 000.0) #2025-04-24 kepler D=30
	#REF_GPU=(1819.55 119.35 118.51 125.09 123.15 143.65 127.16 367.56 390.37 815.74 822.66 129.86 137.22 132.22 145.66 150.15 132.19 131.32 155.11 160.93 156.48 159.38 158.79 156.22 156.45 251.88 0.62 208.72 159.99) #2025-01-27 kepler D=30
	#REF_CPU=(1819.55 120.74 119.19 125.05 123.74 143.48 127.19 368.28 390.32 815.76 822.34 128.82 138.50 132.23 145.43 146.91 132.09 131.42 155.28 160.92 156.67 159.18 000.00 156.36 156.69 282.88 0.66 207.01 000.00) #2025-01-27 kepler D=30
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
fnpre=maos_check_${ARCH}_${D}
fnlog=${fnpre}.all #log of all simulations
fnerr=${fnpre}.err #log of all failed simulations
fntmp=${fnpre}.log #log of current simulation
fnres=${fnpre}.res #result summary
fnref=maos_check.ref #all ference values
if [ -f $fnres ];then
	mv ${fnres} ${fnres}_$(date +%Y%m%d-%H%M%S)
fi
echo $(date) > $fnlog
echo $(date) > $fnerr
ans=0 #result code
ii=0
s_start0=`date +%s`
printf "%-20s    Res    Ref   Ref2     %%     %%  time  total\n" "D=${D}m" | tee $fnres
function run_maos(){
	aotype=$1
	shift
	s_start=`date +%s`
    eval "../bin/maos sim.end=100 $* $args" >$fntmp 2>&1
	a=
    if [ $? -eq 0 ];then
		a=`printf %.1f $(grep 'Mean:' $fntmp |tail -n1 |cut -d ' ' -f 2)`
	fi
	if [ x$a = x ];then
		a=000.0
		ans=$((ans+1)) #failed to run
		echo $aotype $* >> $fnerr
		cat $fntmp >> $fnerr
	fi
	RMS[ii]=$a
    s_end=`date +%s`
	s_diff=$((s_end-s_start)) #time per simulation
	s_diff2=$((s_end-s_start0)) #cumulative time
	echo $aotype $* >> $fnlog
	cat $fntmp >> $fnlog
	rm $fntmp
    b=${REF[$ii]:-0}
    b2=${REF2[$ii]:-0}
	diff=$(echo "200*($a-$b)/($a+$b+1)" | bc)
	diff2=$(echo "200*($a-$b2)/($a+$b2+1)" | bc)
	if [ $diff -gt 10 -o $diff -lt -10 ];then
		ans=$((ans+1)) #mark failure
	fi
	printf "%-20s %6.1f %6.1f %6.1f %5.1f %5.1f %5d %6d\n" "$aotype" "$a" "$b" "$b2" "$diff" "$diff2" "$s_diff" "$s_diff2"| tee -a $fnres
	ii=$((ii+1)) 
}
function run_maos_gpu(){
	if [ "$ARCH" = GPU ];then
		run_maos "$@"
	else
		printf "%-20s   skipped in CPU mode.\n" "$1" | tee -a $fnres
		RMS[ii]=000.0
		ii=$((ii+1)) 
	fi
}
{
run_maos "Openloop:        " -cmcao_lgs.conf sim.evlol=1

run_maos "NGS SCAO (inte): " -cscao_ngs.conf recon.split=0

run_maos "NGS SCAO (ahst): " -cscao_ngs.conf recon.split=1

run_maos "NGS SCAO (lsq):  " -cscao_ngs.conf recon.alg=1

run_maos "NGS SCAO (modal):" -cscao_ngs.conf recon.modal=1 recon.alg=1

run_maos "NGS PWFS (zonal):" -cscao_pywfs.conf recon.modal=0 powfs.dx=[1/16] sim.end=500

run_maos "NGS PWFS (modal):" -cscao_pywfs.conf recon.modal=1 powfs.dx=[1/16] sim.end=500

run_maos "LGS SCAO (inte): " -cscao_lgs.conf recon.split=0

run_maos "LGS SCAO (ahst): " -cscao_lgs.conf recon.split=1

run_maos "LGS GLAO (inte): " -cglao.conf recon.split=0 evl.psfmean=0

run_maos "LGS GLAO (ahst): " -cglao.conf recon.split=1 evl.psfmean=0

run_maos "LGS LTAO (inte): " -cmcao_lgs.conf dm_single.conf fov_oa.conf recon.split=0

run_maos "LGS LTAO (ahst) :" -cmcao_lgs.conf dm_single.conf fov_oa.conf recon.split=1 powfs.astscale=[1 0.1 0.1] #ahst with wide NGS asterism suffers (from aliasing error?)

run_maos "LGS MOAO (ahst): " -cmcao_lgs.conf recon.split=1 evl.moao=0 moao.dx=[1/2]

run_maos "NGS MCAO (inte): " -cmcao_ngs.conf recon.split=0

run_maos "NGS MCAO (ahst): " -cmcao_ngs.conf recon.split=1 

run_maos "LGS MCAO (fit):  " -cmcao_lgs.conf sim.idealtomo=1

run_maos "LGS MCAO (tomo): " -cmcao_lgs.conf evl.tomo=1 

run_maos "LGS MCAO (inte): " -cmcao_lgs.conf recon.split=0 tomo.assemble=1 fit.assemble=1 #also test assembled matrices

run_maos "LGS MCAO (CG):   " -cmcao_lgs.conf tomo.precond=0 tomo.assemble=1 fit.assemble=1 #also test assembled matrices

run_maos "LGS MCAO (FDPCG):" -cmcao_lgs.conf tomo.precond=1 

run_maos "LGS MCAO (CBS):  " -cmcao_lgs.conf tomo.alg=0 fit.alg=0 atmr.os=[2 2 1 1 1 1 1] recon.split=1 #split=0 CBS does not work well in GPU

if [ ${D%.*} -le 10 ];then
run_maos "LGS MCAO (SVD):  " -cmcao_lgs.conf tomo.alg=2 fit.alg=2 atmr.os=[2 2 1 1 1 1 1]
run_maos "LGS MCAO (MVM):  " -cmcao_lgs.conf fit.alg=0 recon.mvm=1
else
run_maos_gpu "LGS MCAO (MVM):  " -cmcao_lgs.conf fit.alg=0 recon.mvm=1
fi
run_maos "LGS MCAO PCCD:  " -cmcao_lgs.conf powfs.radpix=[16,0,0] powfs.pixpsa=[6,0,0]

run_maos "SLGS MCAO (inte): " -cmcao_lgs.conf powfs.fnllt=['llt_SL.conf',,] powfs.pixpsa=[16,0,0] recon.split=0

run_maos "SLGS MCAO (ahst):  " -cmcao_lgs.conf powfs.fnllt=['llt_SL.conf',,] powfs.pixpsa=[16,0,0] recon.split=1 tomo.splitlrt=2

run_maos "LGS NCPA noatm: " -cmcao_lgs.conf sim.noatm=1 ncpa.surf=["'rms=150;L0=20;D=60;SURFEVL=0;'"] sim.wspsd= powfs.noisy=[0] atm.r0z=10 powfs0_llt.fnsurf="'rms=150;mode=5;D=0.5;dx=1/64'" #should be close to 0.

if [ "$D" = 30 ];then
run_maos "NFIRAOS LGS:    " -cnfiraos_lgs.conf
run_maos_gpu "NFIRAOS PYWFS:" -cnfiraos_ngs.conf
fi
echo "REF_${ARCH}=(${RMS[*]}) #$(date +%Y-%m-%d) $HOSTNAME D=$D" | tee -a $fnres >> $fnref
if [ $ans -ne 0 ];then
	echo "$ans tests failed"
fi
};exit $ans 
